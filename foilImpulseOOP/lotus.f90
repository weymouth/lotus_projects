!-------------------------------------------------------!
!--------------- Foil Impulse Test case ----------------!
!-------------------------------------------------------!
program foil_impulse
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use geom_shape  ! to create geometry
  use gridMod,    only: xg,composite
  implicit none
  real,parameter     :: L = 100            ! length
  real,parameter     :: Re =  1000         ! Reynolds number
  integer,parameter  :: b(3) = (/2,2,4/)   ! blocks
  logical,parameter  :: yank=.true.        ! ramp or yank?
  logical,parameter  :: hollow=.false.     ! internal slug flow?
  integer,parameter  :: ndims = 3          ! dimensions
  real,parameter     :: nu = L/Re          ! viscosity
!
  type(fluid)        :: flow
  type(body)         :: foil
  type(set)          :: core
  integer            :: n(3)
  real               :: force(3),area,t0,t1,dt,u=0,dprnt
  logical            :: root
!
! -- Set up run parameters
  real    :: V(3)=0, tStop=8.5, dtPrint=0.5 ! ramp parameters
  real    :: T=10, Umax=1.0, Tend=5, tShift=0
  integer :: dim=1
  if(yank) then                            ! yank parameters
     V = (/1,0,0/); tShift = tStop+0.3; tStop = 0.6; dtPrint = 0.015
     Umax = -6; T = -pi*0.9/Umax; Tend = T
     dim=ndims
  end if
!
! -- Initialize
  call init_mympi(ndims,set_blocks=b(:ndims))
  root = mympi_rank()==0

  n = composite(L*(/3.1,2.6,4./), prnt=root)
  if(ndims==2) n(3) = 1

  call xg(1)%stretch(n(1),-2*L, -0.6*L, L, 5*L, h_max=5., prnt=root)
  call xg(2)%stretch(n(2),-2*L, -0.5*L, 0.5*L, 2*L, h_min=0.5, prnt=root)
  if(ndims==3) call xg(3)%stretch(n(3),-4*L, -2*L, L, 5*L, prnt=root)

  foil = geom('naca_square.IGS')
  area = L
  if(ndims==3) area = L*xg(3)%right
!
! -- Initialize fluid
  call flow%init(n/b, foil, V=V, nu=nu)
  if(root) print *, '-- init complete --'
!
! -- Time update loop
  do while (flow%time/L<tStop+tShift)
     dt = flow%dt/L
     t0 = flow%time/L-tShift ! for motion
     t1 = t0+dt
!
! -- Accelerate the reference frame
     u = uref(t1)
     if(t1<0)    u = uref(0.)
     if(t1>Tend) u = uref(Tend)
     flow%velocity%e(dim)%bound_val = u
     flow%g(dim) = (u-uref(t0))/(dt*L)
     if(t0<0. .or. t0>Tend) flow%g(dim) = 0
     t1 = t1+tShift ! for printing
     if(root) print 1,t1,flow%g(dim),u
1    format("   t=",f0.4," g=",f0.4," u=",f0.4)
!
!-- update and write fluid
     if(hollow.and.yank.and.u.ne.0) call foil%update(t1)
     call flow%update(foil)

     dprnt = merge(0.003, dtPrint, yank.and.t0>0.and.t0<0.3)
     if(mod(t1,dprnt)<dt) call flow%write(lambda=mod(t1,dtPrint)<dt)
!
! -- print force
     force = -foil%pforce(flow%pressure)
     write(9,'(f10.4,f8.4,2e16.8)') t1,dt*L,2.*force(:2)/area
     flush(9)
  end do

  if(root) write(6,*) '--- complete --- '
  call mympi_end
contains
!
! -- Initialize the foil geometry
  type(set) function geom(name)
    character(*)       :: name
    type(model_info)   :: info
    real(8),parameter  :: alpha = 10         ! AOA
!  surface_debug = .true.
    info%file = name
    info%x = (/-4.219,-10.271,-18.876/)
    info%s = 0.36626*L*(/1,1,-1/)
    eps = 2
    info%xmax(1) = L
    info%n = 50
    geom = (model_init(info) &
         .map.(init_affn()**(/alpha,0.D0,0.D0/))) ! rotate by alpha
    core = geom
    if(ndims==3 .and. info%file == 'naca_square.IGS' ) &
         geom = geom.and.plane(4,1,(/0,0,-1/),0,0,0)
    if(hollow) geom = geom.map.init_velocity(slug)
!  call shape_write(100,geom)
  end function geom

  pure function slug(x) result(v)
    real(8),intent(in) :: x(3)
    real(8) :: v(3)
    type(prop) :: p
    p = core%at(x)
    v = 0; v(dim) = u*max(0.,min(-p%distance,1.))
  end function slug

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
