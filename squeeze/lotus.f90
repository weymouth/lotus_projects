program squeeze
  use geom_shape, only: sphere,operator(.map.),init_scale,init_rigid,pi
  use bodyMod,    only: body
  use fluidMod,   only: fluid
  use gridMod,    only: xg
  use imageMod,   only: display
  implicit none
!
! -- Define parameters, declare variables
  real,parameter    :: D = 32, T = 4*D, amp = 0.25 ! length/time scale
  logical,parameter :: sharp = .TRUE.              ! sharp growth
  integer,parameter :: n(3) = D*(/5,3,3/)          ! number of cells in ijk
  real,parameter    :: m = pi/6.*D**3, ma = m/2.   ! mass and added mass guess
  real,parameter    :: k = (2.*pi/T)**2*(m+ma)     ! spring constant for T
  real        :: force(3)=0,pos=-D,vel=0
  real(8)     :: ts
  type(fluid) :: flow
  type(body)  :: geom
!
! -- Initialize grid, body and fluid
  call xg(1)%stretch(n(1),-3.*D,-2*D,2*D,3.*D)  ! center the domain
  call xg(2)%stretch(n(2),-2.*D,-D,D,2.*D)      ! center the domain
  call xg(3)%stretch(n(3),-2.*D,-D,D,2.*D)      ! center the domain
  geom = sphere(2, 1, radius=D/2, center=0) &   ! sphere
          .map.init_scale(0,length,rate) &      ! scale in all dimensions
          .map.init_rigid(1,position,velocity)

  call flow%init(n, geom, nu=0.01)
  flow%dt = min(flow%dt,1.)                     ! limit the time step
!
! -- Time update loop
  do while(flow%time<15*T)
    ts = (flow%time+flow%dt)/T                  ! get non-dimensional time
    call rigid_update(force(1))                 ! update rigid motion
    call geom%update(real(ts))                  ! apply mapping to geom
    call flow%update(geom)                      ! update the flow
    flow%dt = min(flow%dt,1.)                   ! limit the time step
    force = -geom%pforce(flow%pressure)         ! compute the force
!
! -- write to file
    write(9,'(f10.4,f8.4,3e16.8)') ts,flow%dt,force
    if(mod(ts,0.04)<flow%dt/T) then
      call display(flow%velocity%vorticity_Z(),'out_vort')
      write(6,'(f10.4,5f9.4)') ts,flow%dt,length(ts),rate(ts),position(ts),velocity(ts)
    end if
  end do
  call flow%write()
contains
 real(8) pure function length(ts)  ! length scale
    real(8),intent(in) :: ts
    real(8) :: t,c
    t = 2*ts-0.5
    if(sharp) then
      c = 40
      t = mod(ts,0.5)
      t = t+(tanh((0.5-t)*c)-tanh(t*c))/4.
    end if
    length = 1+amp*cos(2.*pi*t)
  end function length
  real(8) pure function rate(ts)  ! rate of change of length
    real(8),intent(in) :: ts
    rate = (length(ts+1e-6)-length(ts-1e-6))/(T*2e-6)
  end function rate
  subroutine rigid_update(force)
    real :: force,accel
    accel = (force-k*pos)/m
    vel = vel+flow%dt*accel
    pos = pos+flow%dt*vel
  end subroutine rigid_update
  real(8) pure function position(ts)
    real(8),intent(in) :: ts
!    position = -D*cos(2.*pi*ts)
    position = pos
  end function position
  real(8) pure function velocity(ts)
    real(8),intent(in) :: ts
!    velocity = (position(ts+1e-6)-position(ts-1e-6))/(T*2e-6)
    velocity = vel
  end function velocity
end program squeeze
