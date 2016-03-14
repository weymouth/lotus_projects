program squeeze
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  use geom_shape, only: sphere,operator(.map.),init_scale,init_rigid,pi
  use bodyMod,    only: body
  use fluidMod,   only: fluid
  use gridMod,    only: xg
  use imageMod,   only: display
  implicit none
!
! -- Define parameters, declare variables
  real,parameter    :: D = 64                      ! sphere size in cells
  real,parameter    :: Re_D = 3200                 ! Reynolds number based on D
  integer           :: n(3) = D*(/8.,4.,1.5/)      ! number of cells in ijk
  integer           :: b(3) = (/4,4,1/)            ! MPI domain cuts in ijk

  real,parameter    :: T = 4*D, amp = AMP          ! size-change length/time scale
  logical,parameter :: sharp = .FLAG.              ! sharp growth flag
  real,parameter    :: m = pi/6.*D**3, ma = m/2.   ! mass and added mass guess
  real,parameter    :: kn = (2.*pi/T)**2*(m+ma)    ! spring constant for T
  real,parameter    :: k = kn                      ! spring constant
  real              :: pos=-POS0*D/2.,vel=0        ! force and motion variables
  real              :: force(3)=0,a0 = 0

  real,parameter    :: lim = 10./D                 ! vorticity level in image
  real,parameter    :: dprnt = 0.04                ! how often to print
  logical     :: root,there = .FALSE.              ! flag for stopping
  real(8)     :: ts
  type(fluid) :: flow
  type(body)  :: geom
!
! -- Initialize mpi, grid, body and fluid
  call init_mympi(ndims=3,set_blocks=b)
  root = mympi_rank()==0

  call xg(1)%stretch(n(1),-5.*D,-3.5*D,3.5*D,5.*D, prnt=root)  ! x (wider for motion)
  call xg(2)%stretch(n(2),-4.*D,-1.5*D,1.5*D,4.*D, prnt=root)      ! y
  call xg(3)%stretch(n(3),0.,0.,D,4.*D, prnt=root)      ! z

  geom = sphere(2, 1, radius=D/2, center=0) &   ! sphere
          .map.init_scale(0,length,rate) &      ! scale in all dimensions
          .map.init_rigid(1,position,velocity)  ! move in x

  call flow%init(n/b, geom, nu=D/Re_D)          ! Initialize the fluid
  flow%dt = min(flow%dt,1.)                     ! limit the time step
!
! -- Time update loop
  do while(flow%time<20*T .and. .not.there)
    ts = (flow%time+flow%dt)/T                  ! get non-dimensional time
    call rigid_update(force(1))                 ! update rigid motion
    call geom%update(real(ts))                  ! apply mapping to geom
    call flow%update(geom)                      ! update the flow
    flow%dt = min(flow%dt,2.)                   ! limit the time step
    force = -2.*geom%pforce(flow%pressure)      ! compute the force
!
! -- write to file
    if(root) write(9,'(f10.4,f8.4,3e16.8)') ts,flow%dt,force
    if(mod(ts,dprnt)<flow%dt/T) then
      call display(flow%velocity%vorticity_Z(),'out_vort', lim=lim)
      call flow%write(geom,lambda=.TRUE.)
      if(root) print '(f10.4,5f9.4)',ts,flow%dt, &
            length(ts),rate(ts),position(ts),velocity(ts)
    end if
    inquire(file='.kill', exist=there)
  end do
contains
 real(8) pure function length(ts)  ! length scale
    real(8),intent(in) :: ts
    real(8) :: t,c
    t = 2.*ts-0.25                         ! pump at double the frequency
    if(sharp) then
      t = mod(ts,0.5)                      ! include shrink phase only
      c = 40                               ! growth rate
      t = t+(tanh((0.5-t)*c)-tanh(t*c))/4. ! add sharp growth phase
    end if
    length = 1+amp*cos(2.*pi*t)
  end function length
  real(8) pure function rate(ts)  ! rate of change of length
    real(8),intent(in) :: ts
    rate = (length(ts+1e-6)-length(ts-1e-6))/(T*2e-6)
  end function rate
  subroutine rigid_update(force)
    real :: force,accel,ma_t
    ma_t = ma*length(ts)**3
    accel = (force-k*pos+ma_t*a0)/(m+ma_t)
    pos = pos+flow%dt*(vel+flow%dt*accel/2.)  ! second order
    vel = vel+flow%dt*accel                   ! first order :(
    a0 = accel
  end subroutine rigid_update
  real(8) pure function position(ts)
    real(8),intent(in) :: ts
    position = pos
  end function position
  real(8) pure function velocity(ts)
    real(8),intent(in) :: ts
    velocity = vel
  end function velocity
end program squeeze
