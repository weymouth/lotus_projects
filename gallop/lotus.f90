program gallop
  use mympiMod,   only: mympi_end,mympi_rank,init_mympi
  use gridMod,    only: xg
  use geom_shape, only: cylinder,operator(.map.),init_rigid,pi
  use bodyMod,    only: body
  use fluidMod,   only: fluid
  use imageMod,   only: display
  implicit none
!
! -- Define parameters, declare variables
  real,parameter    :: L = 32                      ! cylinder length in cells
  real,parameter    :: Re_L = 3200                 ! Reynolds number based on D
  integer           :: n(3) = (/6*L,4*L,1./)        ! number of cells in ijk
  integer           :: b(3) = (/1,1,1/)            ! MPI domain cuts in ijk

  real,parameter    :: Ur = 6                      ! reduced velocity
  real,parameter    :: m = pi/4.*L**2, ma = m      ! mass and added mass guess
  real,parameter    :: T = Ur*L, k=(2.*pi/T)**2*m  ! natural period and spring constant
  real              :: pos=-L/2,vel=0              ! force and motion variables
  real              :: force(3)=0,a0=0
  real(8)           :: ts

  real,parameter    :: lim = 10./L                 ! vorticity level in image
  real,parameter    :: dprnt = 0.04                ! how often to print
  logical           :: root,there = .FALSE.        ! flag for root, stopping
  type(fluid)       :: flow
  type(body)        :: geom
!
! -- Initialize mpi, grid, body and fluid
  call init_mympi(ndims=2,set_blocks=b)
  root = mympi_rank()==0

  call xg(1)%stretch(n(1),-5*L,-0.5*L,2*L,10*L,h_max=6.,prnt=root)  ! x
  call xg(2)%stretch(n(2),-5*L,-L,L,5*L, prnt=root)                  ! y

  geom = cylinder(axis=3, radius=L/2, center=0) & ! circular cylinder
          .map.init_rigid(2,position,velocity)    ! move in y

  call flow%init(n/b, geom, V=(/1.,0.,0./), nu=L/Re_L) ! Initialize the fluid
!
! -- Time update loop
  do while(flow%time<10*T .and. .not.there)
    ts = (flow%time+flow%dt)/T                  ! get non-dimensional time
    call rigid_update(force(2))                 ! update rigid motion
    call geom%update(real(ts))                  ! apply mapping to geom
    call flow%update(geom)                      ! update the flow
    force = -geom%pforce(flow%pressure)         ! compute the force
!
! -- write to file
    if(root) write(9,'(f10.4,f8.4,4e16.8)') ts,flow%dt,force(1:2)/(0.5*L),pos/L,vel
    if(mod(ts,dprnt)<flow%dt/T) then
      call display(flow%velocity%vorticity_Z(),'out_vort', lim=lim)
      if(root) write(6,'(f10.4,f8.4,4e16.8)') ts,flow%dt,force(1:2)/(0.5*L),pos/L,vel
    end if
    inquire(file='.kill', exist=there)
  end do
  call flow%write(geom)
  call mympi_end
contains
  subroutine rigid_update(force)
    real :: force,accel
    accel = (force-k*pos+ma*a0)/(m+ma)        ! update acceleration
    pos = pos+flow%dt*(vel+flow%dt*accel/2.)  ! second order position update
    vel = vel+flow%dt*accel                   ! first order velocity update :(
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
end program gallop
