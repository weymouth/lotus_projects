!-------------------------------------------------------!
!---------------- 2D perching test case ----------------!
!-------------------------------------------------------!
program perch
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  use geom_shape  ! to define geom (set,eps,plane, etc)
  implicit none
  real,parameter     :: f = 1              ! scaling factor
  real,parameter     :: L = 200/f          ! length scale
  real,parameter     :: Re = 2e3           ! Reynolds number
  real,parameter     :: Xi = 1./2.         ! Shape change number
!
  integer,parameter  :: ndims = 2          ! dimensions
  integer,parameter  :: d(3) = (/5,5,1/)   ! domain size
  real,parameter     :: nu = L/Re          ! viscosity
  real,parameter     :: T = L/Xi           ! motion period
  logical,parameter  :: p(3) = (/.false.,.false.,.true./)  ! periodic BCs
  integer            :: b(3) = (/4,4,1/)   ! blocks
  integer            :: n(3),i
  real               :: t0,t1,dt,dtPrint=0.1
  real               :: area,u,u0,a,force(3)
  real               :: Ufric,yp
!
  type(fluid)        :: flow
  type(body)         :: foil
!
! -- Initialize MPI (if MPI is ON)
#if MPION
  call init_mympi(ndims,set_blocks=b(:ndims),set_periodic=p(:ndims))
#else
  b=1
#endif
!
! -- Get array size
  n = composite(L*d/b)
  if(ndims==2) n(3) = 1
!
! -- Get grid
  call xg(1)%init(n(1)*b(1),1.0*L,1.75*L,1.0,c=2.6,f=f)
  call xg(2)%init(n(2)*b(2),0.75*L,2.0*L,1.0,c=2.6,f=f)
!!$  call xg(1)%init(n(1)*b(1),1.0*L,1.75*L,1.0,c=2.8,f=f)
!!$  call xg(2)%init(n(2)*b(2),0.75*L,2.0*L,1.0,c=2.8,f=f)
  if(mympi_rank()==0) call xg(1)%write
  if(mympi_rank()==0) call xg(2)%write
  if(mympi_rank()==0) call xg(3)%write
!
! -- Init
  if(mympi_rank()==0) print *, '-- Perch Test --'
  if(mympi_rank()==0) print '("   L=",f0.4,", points=",i0)', L,product(n*b)
  
  Ufric = sqrt(0.026/Re**(1./7.)/2.)
  yp = Ufric/nu
  if(mympi_rank()==0) print '("   nu=",f0.4,", y+=",f0.4)', nu,yp

!
! -- Initialize the foil geometry
  foil = naca(L,0.16).map.init_rigid(6,alpha,omega)
  area = L
  if(ndims==3) area = L*n(3)*b(3)*xg(3)%h
!
! -- Initialize fluid
  call flow%init(n,foil,V=(/1.,0.,0./),nu=nu)
  if(flow%time==0.) flow%time = -5*L
  if(mympi_rank()==0) print *, '-- init complete --'
!
! -- Time update loop
  do while (flow%time/T<1.2)
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
     if(mympi_rank()==0 .and. mod(abs(t1),dtPrint)<dt) &
          print 1,t1,a,u,alpha(REAL(t1,8)),omega(REAL(t1,8))
1    format("   t=",f0.4," g=",f0.4," u=",f0.4," alpha=",f0.4," omega=",f0.4)
!
!-- update and write fluid
     call flow%update(foil)
     if(mod(t1,dtPrint)<dt .and. t1>-dt) call flow%write(foil)
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
    info%xmax(1) = c
    info%n = (/c,c,1./)
    shift = init_affn()+REAL(c*(/0.5-o,0.,0./),8) ! shift pivot to origin
!    surface_debug = .true.
    eps = 2.0
    naca = model_init(info).map.shift
  end function naca
!
! -- motion definitions
  real(8) pure function alpha(ts)  ! rotation angle
    real(8),intent(in) :: ts
    if(ts<0) then
       alpha = 0.
    else if(ts>1) then
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
