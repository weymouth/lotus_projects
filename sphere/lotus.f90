program sphere_flow
  use mympiMod,   only: init_mympi
  use gridMod,    only: xg
  use geom_shape, only: sphere
  use bodyMod,    only: body
  use fluidMod,   only: fluid
  use imageMod,   only: display
  implicit none
!
! -- Define parameter, declare variables
  real,parameter :: D = 4*32, Re = 30e3, dprnt = 0.5, lim = 20./D
  integer :: n(3) = 32*(/4,2,2/), box(4) = (/-100,-320,2080,640/), b(3) = (/2,1,1/)
  logical :: there = .false.
  type(fluid) :: flow
  type(body) :: geom
  call init_mympi(ndims=3,set_blocks=b)
!
! -- Initialize grid, body and fluid
  call xg(1)%stretch(n(1), -2*D, -0.5*D, 0.5*D, 5*D, h_min=4., h_max=12.)
  call xg(2)%stretch(n(2), -2*D, -0.5*D, 0.5*D, 2*D, h_min=4.)
  call xg(3)%stretch(n(3), -2*D, -0.5*D, 0.5*D, 2*D, h_min=4.)
  geom = sphere(2, 1, radius=0.5*D, center=0)
  call flow%init(n/b,geom,V=(/1.,0.,0./),nu=D/Re)
!
! -- Time update loop
  do while(.not.there)
    call flow%update
    if(mod(flow%time,D*dprnt)<flow%dt) &
      call display(flow%velocity%vorticity_Z(),'vort', lim=lim, box=box)
    inquire(file='.kill', exist=there)
  end do
  call flow%write()
end program sphere_flow
