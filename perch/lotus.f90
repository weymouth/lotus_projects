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
  integer,parameter  :: d(3) = (/8,8,1/)   ! domain size
!
  integer,parameter  :: ndims = 2   ! dimensions
  real(8),parameter  :: L = f       ! length
  integer,parameter  :: m(3) = f*d  ! points
  real,parameter     :: nu = L/Re   ! viscosity
  real,parameter     :: T = L/Xi    ! motion period
  integer            :: n(3)
  real               :: t0,t1,dt,dtPrint=0.1
  real               :: area,u,u0,force(3)
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
  xg(1:2)%left = -m(1:2)
!
! -- Init
  if(mympi_rank()==0) print *, '-- Perch Test --'
  if(mympi_rank()==0) print '("   L=",i0," nu=",f0.4)', f,nu
!
! -- Initialize the foil geometry
  info%file = 'input.IGS'
  info%x = (/-4.219,-10.271,-18.876/)
  info%s = 0.36626*L*(/1,1,-1/)
  foil = model_init(info).map.jerk(5,6,0,T,0,0)
  area = L
  if(ndims==3) area = L*m(3)
!
! -- Initialize fluid
  call flow%init(n,foil,nu=nu)
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
     u = 1; u0 = 1
     if(t0>0) u0 = 1-t0
     if(t1>0) u  = 1-t1
     if(t0>1) u0 = 0
     if(t1>1) u  = 0
     flow%velocity%e(1)%bound_val = u
     flow%g(1) = (u-u0)/flow%dt
     if(mympi_rank()==0) print 1,t1,flow%g(1),u
1    format("   t=",f0.4," g=",f0.4," u=",f0.4)
!
!-- update and write fluid
     call flow%update
     if(mod(t1,dtPrint)<dt) call flow%write
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
end program perch
