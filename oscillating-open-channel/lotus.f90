!-------------------------------------------------------!
!---------------- Turbulent Channel Flow ---------------!
!-------------------------------------------------------!
program channel
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use gridMod,    only: xg,composite
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use imageMod,   only: display
  use geom_shape
  implicit none
! - Domain
  integer            :: b(3) = [4,32,4], n(3) ! blocks
  real               :: m(3) = [1/512., 1/64., 1/512.]
! - Length scales
  integer,parameter  :: a = 2**16
  real,parameter     :: nu = a/102000
  real, parameter    :: delta=a/226
  real               :: h_roughness=delta,thickness,lambda=delta
! - Time scales
  real,parameter     :: St=0.318, f=St/(2.*a)
  real,parameter     :: finish=20., print_int=0.05, init_time=finish-1.
  real               :: t,dt
! - Init sim
  type(fluid)        :: flow
  type(body)         :: wall
  logical            :: root,there=.false.
! - io
  real               :: enstrophy, tke
! -- Initialize MPI
  call init_mympi(3,set_blocks=b,set_periodic=[.true.,.false.,.true.])
  root = mympi_rank()==0
  !
  ! -- Init grid
  if(root) print *,'Setting up the grid, body and fluid'
  if(root) print *,'-----------------------------------'
  n = composite(a*m,prnt=root)
  xg(1)%h = (4.*lambda)/n(1)
  call xg(2)%stretch(n(2),0.,0.,1.*delta,10.*delta,prnt=root)
  xg(3)%h = (4.*lambda)/n(3)
  if(root) print *,'x,z aspect ratio', xg(1)%h 
!
! -- Init channel
  wall = upper(2.*h_roughness+2).map.init_rigid(1,position,velocity)
  eps=1.0
!
! -- Initialize fluid
  call flow%init(n/b,wall,nu=nu,exit=.true.)
  flow%dt = min(flow%dt,1.)
!
! -- Time update loop
  if(root) print *,'Starting time update loop'
  if(root) print *,'-----------------------------------'
!
  do while (flow%time<finish/f.and..not.there)
    t = flow%time
    dt = flow%dt
    call wall%update(t+dt)
    call flow%update(wall)
    if(mod(t,print_int/f)<dt) then
      if(root) print "('Time:',f15.3)",t*f
      if(t>init_time/f) call flow%write(wall, write_vtr=.true.)
    end if
  ! - Write some stuff
    tke = flow%velocity%tke(mean=.true.)
    enstrophy = flow%velocity%enstrophy(mean=.true.)
    if(root) write(9,'(f10.4,4e16.8)') t*f, &
            enstrophy, tke
  end do
  !
  if(root) print *,'Loop complete: writing restart files and exiting'
  if(root) print *,'-----------------------------------'
  call mympi_end()
contains
!
! -- Kinematics
real(8) pure function position(t)
  real(8),intent(in) :: t
  position = a*sin(2*pi*f*t)
end function position
real(8) pure function velocity(t)
  real(8),intent(in) :: t
  velocity = 2*pi*f*a*cos(2*pi*f*t)
end function velocity
!
type(set) function upper(thickness) result(geom)
  real,intent(in) :: thickness
  geom = plane([0.,1.,0.],[0.,thickness,0.]).map.init_warp(2,egg_top,dotegg_top,degg_top)
end function
!
real pure function egg_top(x)
  real,intent(in) :: x(3)
  egg_top = (h_roughness)*sin((2*pi*x(1)/lambda)-pi/2)*cos(2*pi*x(3)/lambda)
end function egg_top
pure function degg_top(x)
  real,intent(in) :: x(3)
  real            :: degg_top(3)
  degg_top = 0
  degg_top(1) = (h_roughness)*cos((2*pi*x(1)/lambda)-pi/2)&
                                *sin((2*pi*x(3)/lambda))*(2*pi/lambda)
  degg_top(3) = -(h_roughness)*cos((2*pi*x(1)/lambda-pi/2))&
                                *sin((2*pi*x(3)/lambda))*(2*pi/lambda)
end function degg_top
real pure function dotegg_top(x)
  real,intent(in) :: x(3)
  dotegg_top = 0
end function dotegg_top
end program channel
