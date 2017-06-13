!-------------------------------------------------------!
!---------------- Static square cylinder ---------------!
!-------------------------------------------------------!
program square_abreast
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  use imageMod,   only: display
  use geom_shape
  implicit none
  real,parameter :: L = 32                   ! length scale
  real,parameter :: Re = 100                 ! Reynolds number
  real,parameter :: g = 3                    ! gap
  real,parameter :: f = 0.15                 ! pertubation frequency
!
  real,parameter :: nu = L/Re                ! viscosity
  real,parameter :: m(3) = (/24*L,12*L,1./)  ! points (approx)
  integer        :: b(3) = (/4,4,1/)         ! blocks
!
  integer        :: n(3)
  logical        :: root,there=.false.
  type(fluid)    :: flow
  type(body)     :: square
!
! -- Initialize MPI (if MPI is OFF, b is set to 1)
  call init_mympi(2,set_blocks=b(:2))
  root = mympi_rank()==0
!
! -- Print run info
  if(root) print *, '-- Square Abreast Cylinders --'
!
! -- Initialize array size
  n = composite(m,prnt=root)
!
! -- Initialize and print grid
  call xg(1)%stretch(n(1),-15.5*L,-1.5*L,10*L,45.5*L,h_max=5.,prnt=root)
  call xg(2)%stretch(n(2),-22.5*L,-(g+5)*L/2.,(g+5)*L/2.,22.5*L,prnt=root)
!
! -- Initialize the square geometry
  square = (plane(4,1,(/-1,0,0/),(/-0.5*L,0.,0./),0,0) &
       .and.plane(4,1,(/1,0,0/),(/0.5*L,0.,0./),0,0) &
       .and.plane(4,1,(/0,1,0/),(/0.,(g/2.+1)*L,0./),0,0) &
       .and.plane(4,1,(/0,-1,0/),(/0.,-(g/2.+1)*L,0./),0,0) &
          -(plane(4,1,(/0,-1,0/),(/0.,-g/2.*L,0./),0,0) &
       .and.plane(4,1,(/0,1,0/),(/0.,g/2.*L,0./),0,0))) &
       .map.init_rigid(2,position,velocity)  ! move in y
!
! -- Initialize fluid
  call flow%init(n/b,square,V=(/1.,0.,0./),nu=nu, exit=.true.)
  if(flow%time==0) flow%time = -10*L/f
  if(root) print *, '-- init complete --'
!
! -- Time update loop
  do while (flow%time/L<600.and..not.there)
     if(flow%time<0.) then
       call square%update(flow%time+flow%dt)
       call flow%update(square)
     else
       call flow%update()
     end if
!
! -- Print fields and force
     write(9,'(f10.4,f8.4,3e14.6)') flow%time/L,flow%dt, &
        -2./L*square%pforce(flow%pressure)
     if(mod(abs(flow%time),L)<flow%dt) then
        flush(9)
        if(root) print *,flow%time/L
        call display(flow%velocity%vorticity_Z(),'01vort', &
                     lim=4./L,box=(/-100,-320,6400,640/))
     end if
     inquire(file='.kill', exist=there)
  end do
  call flow%write()
  if(root) print *, '--- complete ---'
  call mympi_end()
contains
  real(8) pure function position(t) ! location
     real(8),intent(in) :: t
     position = 0.1*L*(1-cos(2.*pi*min(t/L*f,0.)))
   end function position
   real(8) pure function velocity(t)
     real(8),intent(in) :: t
     velocity = (position(t+1e-6)-position(t-1e-6))/2e-6
   end function velocity
end program square_abreast
