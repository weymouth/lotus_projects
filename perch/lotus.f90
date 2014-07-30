!-------------------------------------------------------!
!---------------- 2D perching test case ----------------!
!-------------------------------------------------------!
program perch
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use geom_shape  ! to create geometry
  use gridMod,    only: xg
  implicit none
  integer,parameter  :: f=3*2**5           ! resolution  
  real,parameter     :: Re = 2e3           ! Reynolds number
  real,parameter     :: Xi = 0.125         ! Shape change number
  integer,parameter  :: b(3) = (/4,4,1/)   ! blocks
  integer,parameter  :: d(3) = (/6,6,1/)   ! domain size
!
  integer,parameter  :: ndims = 2   ! dimensions
  real,parameter     :: L = f       ! length
  integer,parameter  :: m(3) = f*d  ! points
  real,parameter     :: nu = L/Re   ! viscosity
  real,parameter     :: T = L/Xi    ! motion period
  integer            :: n(3)
  real               :: t0,t1,dt,dtPrint=0.1
  real               :: area,u,u0,a,force(3)
!
  type(fluid)        :: flow
  type(body)         :: foil
  type(model_info)   :: info
!
! -- Initialize MPI (if MPI is ON)
#if MPION
  call init_mympi(ndims,set_blocks=b(:ndims))
!
! -- Get grid size
  n = m/b
#else
  n = m
#endif
  if(ndims==2) n(3) = 1
  xg(1:2)%left = -m(1:2)/2.
  if(mympi_rank()==0) call xg(1)%write
  if(mympi_rank()==0) call xg(2)%write
  if(mympi_rank()==0) call xg(3)%write
!
! -- Init
  if(mympi_rank()==0) print *, '-- Perch Test --'
  if(mympi_rank()==0) print '("   L=",i0," nu=",f0.4)', f,nu
!
! -- Initialize the foil geometry
  foil = naca(L,0.12,(/0.16,0.,0./))
  area = L
  if(ndims==3) area = L*m(3)
!
! -- Initialize fluid
  call flow%init(n,foil,V=(/1.,0.,0./),nu=nu)
  flow%velocity%e%exit = .true.
  call flow%resume
  if(flow%time==0.) flow%time = -5*L
  if(mympi_rank()==0) print *, '-- init complete --'
!
! -- Time update loop
  do while (flow%time/T<2)
     dt = flow%dt/T
     t0 = flow%time/T
     t1 = t0+dt
!
! -- Deccelerate
     u = 1.; u0 = 1
     if(t0>0) u0 = 1.-t0
     if(t0>1) u0 = 0.
     if(t1>0) u = 1.-t1
     if(t1>1) u = 0.
     flow%velocity%e(1)%bound_val = u
     a = (u-u0)/flow%dt
     flow%g(1) = a

     if(mympi_rank()==0) print 1,t1,a,u
1    format("   t=",f0.4," g=",f0.4," u=",f0.4)
!
!-- update and write fluid
     call flow%update
     if(mod(abs(t1),dtPrint)<dt) call flow%write
!
! -- print force
     force = foil%pforce(flow%pressure)
     write(9,'(f10.4,f8.4,3e16.8)') t1,dt*T,2.*force/area
     flush(9)
  end do
  if(mympi_rank()==0) write(6,*) '--- complete --- '
!
! -- Finalize MPI
#if MPION
  call mympi_end
#endif

contains
  type(set) function naca(c,t,o)
    real,intent(in)   :: c,t,o(3)
    integer :: m = 20,i
    real    :: p(3)=0,q(3)=0,n(3)

    p = 0; 
    q(1) = (1./m)**2; q(2) = t*oord(q(1))
    print *,q
    n = (/p(2)-q(2),q(1)-p(1),0./)
    n = n/sqrt(sum(n*n))
    naca = plane(4,1,n,c*(p-o),0,0).and.plane(4,1,mirror(n),mirror(c*(p-o)),0,0)

    do i=2,m
       p = q
       q(1) = (real(i)/m)**2; q(2) = t*oord(q(1))
       print *,q
       n = (/p(2)-q(2),q(1)-p(1),0./)
       n = n/sqrt(sum(n*n)) 
       naca = naca.and.plane(4,1,n,c*(p-o),0,0).and.plane(4,1,mirror(n),mirror(c*(p-o)),0,0)
    end do

  end function naca

  function mirror(a) result(b)
    real,intent(in) :: a(3)
    real :: b(3)
    b = a; b(2) = -b(2)
  end function mirror

  real function oord(x) result(y)
    real,intent(in) :: x
    y = 5*(0.2969*sqrt(x)-0.1260*x-0.3516*x**2+0.2843*x**3-0.1015*x**4);
  end function oord
end program perch
