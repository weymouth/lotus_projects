!-------------------------------------------------------!
!--------------- Foil Impulse Test case ----------------!
!-------------------------------------------------------!
program foil_impulse
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use geom_shape  ! to create geometry
  use gridMod,    only: xg
  implicit none
  integer,parameter  :: L0 = 2**6          ! base resolution  
  real,parameter     :: f = 3.0            ! refinement level
  real,parameter     :: Re = 1.4e4         ! Reynolds number
  real,parameter     :: alpha = 10         ! AOA
  integer,parameter  :: b(3) = (/2,4,2/)   ! blocks
  integer,parameter  :: d(3) = (/4,4,6/)   ! domain size
  logical,parameter  :: yank=.false.       ! ramp or yank?
!
  real,parameter     :: L = L0*f    ! resolution
  integer,parameter  :: ndims = 2   ! dimensions
  integer,parameter  :: m(3) = L*d  ! points
  real,parameter     :: nu = L/Re   ! viscosity
  integer            :: n(3)
  real               :: force(3),area,u,t0,t1,dt
!
  type(fluid)        :: flow
  type(body)         :: foil
  type(set)          :: geom
  type(model_info)   :: info
!
! -- Set up run parameters
  real    :: V(3)=0, tStop=30, dtPrint=3   ! ramp parameters
  real    :: T=10, Umax=1.0, Tend=5, tShift=0
  integer :: dim=1
  if(yank) then                            ! yank parameters
     V = (/1,0,0/); tShift = tStop; tStop = 2; dtPrint = 0.02
     Umax = -6; T = -pi*0.9/Umax; Tend = T
     dim=ndims
  end if
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
  call xg(1)%init(m(1),0.55*L,2.00*L,0.50,r=1.0235,f=1/f,h=1.0,c=2.)
  call xg(2)%init(m(2),0.25*L,0.25*L,1.00,r=1.0235,f=1/f,h=0.5,c=8.)
  if(ndims==3) call xg(3)%init(m(3),2.0*L,1.0*L,1.0,f=1/f)
  if(mympi_rank()==0) print *, '-- Foil Impulse --'
  if(mympi_rank()==0) print '("   L=",i0," nu=",f0.4)', int(L),nu
  if(mympi_rank()==0) call xg(1)%write
  if(mympi_rank()==0) call xg(2)%write
  if(mympi_rank()==0) call xg(3)%write
!
! -- Initialize the foil geometry
  surface_debug = .true.
  info%file = 'naca_square.IGS'
  info%x = (/-4.219,-10.271,-18.876/)
  info%s = 0.36626*L*(/1,1,-1/)
  if(ndims==2) info%n = (/L,L,1./) 
  geom = (model_init(info) &
       .map.(init_affn()**real((/alpha,0.,0./),8))) ! rotate by alpha
  if(ndims==3) geom = geom.and.plane(4,1,(/0,0,-1/),0,0,0)
  foil = geom
  area = L
  if(ndims==3) area = L*xg(3)%right
!
! -- Initialize fluid
  call flow%init(n,foil,V=V,nu=nu)
  flow%velocity%e%exit = .true.
  call flow%resume
  if(mympi_rank()==0) print *, '-- init complete --'
!
! -- Time update loop
  do while (flow%time/L<tStop+tShift)
     dt = flow%dt/L
     t0 = flow%time/L-tShift ! for motion
     t1 = t0+dt
!
! -- Accelerate the reference frame
     u = uref(t1)
     if(t1>Tend) u = uref(Tend)
     flow%velocity%e(dim)%bound_val = u
     flow%g(dim) = (u-uref(t0))/(dt*L)
     if(t0>Tend) flow%g(dim) = 0
     t1 = t1+tShift ! for printing
     if(mympi_rank()==0 .and. mod(t1,dtPrint)<dt) &
          print 1,t1,flow%g(dim),u
1    format("   t=",f0.4," g=",f0.4," u=",f0.4)
!
! -- shift grid
     if(yank) then
        call xg(1)%shift(-V(1)*flow%dt)
        call xg(dim)%shift(-u*flow%dt)
     end if
!
!-- update and write fluid
     call flow%update
     if(mod(t1,dtPrint)<dt) call flow%write
!
! -- print force
     force = foil%pforce(flow%pressure)
     write(9,'(f10.4,f8.4,3e16.8)') t1,dt*L,2.*force/area
     flush(9)
  end do
  if(mympi_rank()==0) write(6,*) '--- complete --- '
!
! -- Finalize MPI
#if MPION
  call mympi_end
#endif
contains
  real function uref(time)
    real,intent(in) :: time
    if(time/T*10<0.5) then
       uref = sin(time/T*pi*10)
    else if((T-time)/T*10<0.5) then
       uref = sin((T-time)/T*pi*10)
    else
       uref = 1.
    end if
    uref = Umax*uref*sin(time/T*pi)
  end function uref
end program foil_impulse
!
