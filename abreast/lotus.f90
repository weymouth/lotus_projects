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
  real,parameter :: L = 50                   ! length scale
  real,parameter :: Re = 100                 ! Reynolds number
  real,parameter :: g = 3                    ! gap
!
  real,parameter :: nu = L/Re                ! viscosity
  real,parameter :: m(3) = (/25*L,10*L,1./)  ! points (approx)
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
  if(root) print *, '-- Square Cylinder --'
!
! -- Initialize array size
  n = composite(m,prnt=root)
!
! -- Initialize and print grid
  call xg(1)%stretch(n(1),-15*L,-L,L,100*L,h_max=5.,prnt=root)
  call xg(2)%stretch(n(2),-15*L,-(g/2.+1.5)*L,(g/2.+1.5)*L,15*L,prnt=root)
!
! -- Initialize the square geometry
  square = plane(4,1,(/-1,0,0/),(/-0.5*L,0.,0./),0,0) &
       .and.plane(4,1,(/1,0,0/),(/0.5*L,0.,0./),0,0) &
       .and.plane(4,1,(/0,1,0/),(/0.,(g/2.+1)*L,0./),0,0) &
       .and.plane(4,1,(/0,-1,0/),(/0.,-(g/2.+1)*L,0./),0,0) &
          -(plane(4,1,(/0,-1,0/),(/0.,-g/2.*L,0./),0,0) &
       .and.plane(4,1,(/0,1,0/),(/0.,g/2.*L,0./),0,0))
!
! -- Initialize fluid
  call flow%init(n/b,square,V=(/1.,0.,0./),nu=nu)
  if(root) print *, '-- init complete --'
!
! -- Time update loop
  do while (.not.there)
     call flow%update()
!
! -- Print fields and force
     write(9,'(f10.4,f8.4,3e14.6)') flow%time/L,flow%dt, &
        -2./L*square%pforce(flow%pressure)
     flush(9)
     if(mod(flow%time,0.5*L)<flow%dt) call display( &
        flow%velocity%vorticity_Z(),'vort',lim=4./L,box=(/-100,-320,6400,640/))
     inquire(file='.kill', exist=there)
  end do
  call flow%write()
  if(root) print *, '--- complete ---'
  call mympi_end()
end program square_abreast
