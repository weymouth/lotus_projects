!-------------------------------------------------------!
!------------------- Cylinder Array  -------------------!
!-------------------------------------------------------!
program array
  use fluidMod, only: fluid
  use bodyMod,  only: body
  use mympiMod, only: init_mympi,mympi_end,mympi_rank
  use gridMod,  only: xg,composite
  use geom_shape
  implicit none

  integer,parameter :: rows = 4, finish = 1021, dtPrint = 100
  real,parameter    :: D = 16, DG = 21*D, Re_D = 100, nu = D/Re_D
  integer,parameter :: b(3) = (/4,4,1/)

  type(fluid) :: flow
  type(body)  :: bodies
  integer     :: n(3)
  real        :: dt, t1, pforce(2), vforce(2)
  logical     :: root
!
! -- Initialize simulation
  call init_mympi(2,set_blocks=b(:2))
  root = mympi_rank()==0

  n = composite(DG*(/4,2,0/), root)
  call xg(1)%stretch(n(1), -10*DG, -0.6*DG, 0.6*DG, 25*DG, h_max=15., prnt=root)
  call xg(2)%stretch(n(2), -10*DG, -0.6*DG, 0.6*DG, 10*DG, prnt=root)

  bodies = make_array()
  call flow%init(n/b, bodies, nu=nu)
  if(flow%time==0.) call flow%write(bodies)

  if(root) print *, '-- init complete --'
!
! -- Time update loop
  do while (flow%time/D<finish)
     dt = flow%dt
     t1 = flow%time+dt
!
! -- smoothly increase flow speed
     flow%velocity%e(1)%bound_val = tanh(t1/(50*D))
     flow%g(1) = sech(t1/(50*D))**2/(50*D)
!
! -- update flow
     call flow%update()
!
! -- print force coefficients
     pforce = -bodies%pforce(flow%pressure)
     vforce = nu*bodies%vforce(flow%velocity)
     write(9,'(f10.4,f8.4,4e16.8)') t1/D,dt,2.*pforce/DG,2.*vforce/DG
     flush(9)
!
! -- print flow
     if(mod(t1,dtPrint*D)<dt) call flow%write()
  end do

  call mympi_end

contains
!
! -- hyperbolic secant that is stable for large x
  real pure function sech(x)
    real,intent(in) :: x
    sech = 2*exp(-x)/(1+exp(-2*x))
  end function sech
!
! -- make an array of cylinders with group diameter DG
  type(set) function make_array()
    integer,parameter :: n(4) = (/6,13,19,25/)
    real,parameter    :: R = 0.5*DG-0.5*D
    real    :: theta,xc,yc
    integer :: i,j
    make_array = place_cyl(0.,0.)
    do j=1,rows
       do i=1,n(j)
          theta = 2.*3.14159*(i-1.)/real(n(j))
          xc = R*sin(theta)*real(j)/rows
          yc = R*cos(theta)*real(j)/rows
          make_array = make_array.or.place_cyl(xc,yc)
       end do
    end do
  end function make_array
  type(cylinder) function place_cyl(xc,yc)
    real,intent(in) :: xc,yc
    place_cyl = cylinder(1,1,3,0.5*D,center=(/xc,yc,0./))
  end function place_cyl
end program array
