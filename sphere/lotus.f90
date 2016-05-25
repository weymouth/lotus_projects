program sphere_flow
  use gridMod,    only: xg
  use geom_shape, only: sphere,pi
  use bodyMod,    only: body
  use fluidMod,   only: fluid
  use imageMod,   only: display
  implicit none
!
! -- Define parameter, declare variables
  real,parameter :: D = 2*32, Re = 30e3 ! diameter and Reynolds number
  integer :: n(3) = 32*(/4,2,2/)        ! numer of points
  type(fluid) :: flow
  type(body) :: geom
!
  print *,'Setting up the grid, body and fluid'
  print *,'-----------------------------------'
  call xg(1)%stretch(n(1), -2*D, -0.5*D, 0.5*D, 5*D, h_min=2., h_max=12.,prnt=.true.)
  call xg(2)%stretch(n(2), -2*D, -0.5*D, 0.5*D, 2*D, h_min=2.,prnt=.true.)
  call xg(3)%stretch(n(3), -2*D, -0.5*D, 0.5*D, 2*D, h_min=2.,prnt=.true.)
  geom = sphere(2, 1, radius=0.5*D, center=0)
  call flow%init(n,geom,V=(/1.,0.,0./),nu=D/Re)
!
  print *,'Starting time update loop'
  print *,'-----------------------------------'
  print *,' -t- , -dt- '
  do while(flow%time<10*D)
    call flow%update
    write(9,'(f10.4,f8.4,3e16.8)') flow%time/D,flow%dt, &
          2.*geom%pforce(flow%pressure)/(pi*D**2/4.)
    if(mod(flow%time,0.1*D)<flow%dt) then
      print '(f6.1,",",f6.3)',flow%time/D,flow%dt
      call display(flow%velocity%vorticity_Z(), 'out_vort', lim = 20./D)
    end if
  end do
  print *,'Loop complete: writing restart files and exiting'
  print *,'-----------------------------------'
  call flow%write()
end program sphere_flow
