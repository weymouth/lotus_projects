!-------------------------------------------------------!
!---------------- 2D perching test case ----------------!
!-------------------------------------------------------!
program perch
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg
  use geom_shape  ! to define geom (set,eps,plane, etc)
  implicit none
  integer,parameter  :: f=3*2**5           ! resolution  
  real,parameter     :: Re = 2e3           ! Reynolds number
  real,parameter     :: Xi = 0.5           ! Shape change number
  integer,parameter  :: b(3) = (/4,4,1/)   ! blocks
  integer,parameter  :: d(3) = (/4,4,1/)   ! domain size
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
  call xg(1)%init(m(1),0.75*L,1.25*L,1.0,r=1.018)
  call xg(2)%init(m(2),0.25*L,1.75*L,1.0,r=1.018)
  if(mympi_rank()==0) call xg(1)%write
  if(mympi_rank()==0) call xg(2)%write
  if(mympi_rank()==0) call xg(3)%write
!
! -- Init
  if(mympi_rank()==0) print *, '-- Perch Test --'
  if(mympi_rank()==0) print '("   L=",i0," nu=",f0.4)', f,nu
!
! -- Initialize the foil geometry
  foil = naca(L,0.15).map.init_rigid(6,alpha,omega)
  area = L
  if(ndims==3) area = L*m(3)
!
  ! -- Initialize fluid
  call flow%init(n,foil,V=(/1.,0.,0./),nu=nu)

  flow%velocity%e%exit = .true.
!  if(ndims==3) call flow%velocity%e(1)%perturb(0.05)
  if(flow%time==0.) flow%time = -L*5
  if(mympi_rank()==0) print *, '-- init complete --'
  !
! -- Time update loop
  do while (flow%time/T<2)
     dt = flow%dt/T
     t0 = flow%time/T
     t1 = t0+dt
!
! -- Deccelerate
     u  = min(1.,max(0.,1.-t1))
     u0 = min(1.,max(0.,1.-t0))
     a = (u-u0)/flow%dt
     flow%velocity%e(1)%bound_val = u
     flow%g(1) = a
!
! -- rotate
     if(t1>0 .and. t0<1) then
        call foil%update(t1)
     end if
     if(mympi_rank()==0) print 1,t1,a,u,alpha(DBLE(t1)),omega(DBLE(t1))
1    format("   t=",f0.4," g=",f0.4," u=",f0.4," alpha=",f0.4," omega=",f0.4)
!
!-- update and write fluid
     call flow%update(foil)
     if(mod(t1,dtPrint)<dt .and. t1>0) call flow%write(foil)
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
!
! -- geometry definitions
  type(set) function naca(c,o)
    real,intent(in) :: c,o
    type(model_info) :: info
    type(affn) :: shift
    info%file = 'naca_square.IGS'
    info%x = (/-4.219,-10.271,-18.876/)
    info%s = 0.36626*c*(/1,1,-1/)
    info%n = (/c,c,1./)
    shift = init_affn()+DBLE(c*(/0.5-o,0.,0./)) ! shift pivot to origin
    surface_debug = .true.
    eps = 2.0
    naca = model_init(info).map.shift
  end function naca
!
! -- motion definitions
  real(8) pure function alpha(ts)  ! rotation angle
    real(8),intent(in) :: ts
    if(ts>1) then
       alpha = pi/2.
    else if(ts>0) then
       alpha = pi/2.*(ts-sin(2.*pi*ts)/(2.*pi))
    end if
  end function alpha
  real(8) pure function omega(ts)  ! rotation velocity
    real(8),intent(in) :: ts
    omega = 0
    if(ts>0 .and. ts<1) then
       omega = pi/2.*(1-cos(2.*pi*ts))/T
    end if
  end function omega
end program perch
