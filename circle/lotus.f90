program circle_flow
  use bodyMod,    only: body
  use fluidMod,   only: fluid
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  use imageMod,   only: display
  use geom_shape
  implicit none
!
! -- Define parameter, declare variables
  real,parameter :: D = 128, Re = 100, nu = D/Re ! diameter and Reynolds number
  integer :: n(3), b(3) = (/4,4,1/), box(4) = (/-2*D,-2*D,7.5*D,4*D/)   ! MPI blocks (product must equal n_procs)
  logical :: root, there = .false.     ! root processor
  type(fluid) :: flow
  type(body) :: geom
  integer :: i
!
! -- Initialize MPI (if MPI is OFF, b is set to 1)
  call init_mympi(2,set_blocks=b)
  root = mympi_rank()==0
!
  if(root) print *,'Setting up the grid, body and fluid'
  if(root) print *,'-----------------------------------'
  n = composite((/8.*D,4.*D,0./),prnt=root)
  call xg(1)%stretch(n(1), -10*D, -0.6*D, 3*D, 20*D, h_max=8.,prnt=root)
  call xg(2)%stretch(n(2), -15*D, -0.9*D, 0.9*D, 15*D, prnt=root)
  geom = cylinder(axis=3, radius=0.5*D, center=0)!.and.plane(norm=(/1,0,0/),center=0.)
  call flow%init(n/b,geom,V=(/1.,0.,0./),nu=nu,exit=.true.)
  call display(flow%velocity%vorticity_Z(), 'out_vort', lim = 0.25, box=box)
!
  if(root) print *,'Starting time update loop'
  if(root) print *,'-----------------------------------'
  if(root) print *,' -t- , -dt- '
  do while(flow%time<300*D.and..not.there)
    call flow%update()
    write(9,'(f10.4,f8.4,6e16.8)') flow%time/D,flow%dt, &
       2./D*nu*geom%vforce(flow%velocity), &
      -2./D*geom%pforce(flow%pressure)
      ! (flow%velocity%e(i)%max(),flow%velocity%e(i)%max(neg=.TRUE.),0.,i=1,2)
    flush(9)
    if(mod(flow%time,1*D)<flow%dt) then
      if(root) print '(f6.1,",",f6.3)',flow%time/D,flow%dt
      call display(flow%velocity%vorticity_Z(), 'out_vort', lim = 0.25, box=box)
    end if
    inquire(file='.kill', exist=there)
  end do
  if(root) print *,'Loop complete: writing restart files and exiting'
  if(root) print *,'-----------------------------------'
  call flow%write()
  call mympi_end
end program circle_flow
