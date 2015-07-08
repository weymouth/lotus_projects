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

  integer,parameter :: rows = 2, ndims = 3, finish = 600
  real,parameter    :: D = 16, DG = 7*D, Re_D = 100, nu = D/Re_D
  integer,parameter :: b(3) = (/4,4,1/), dtPrint = 10
  logical,parameter :: p(3) = (/.false.,.false.,.true./), single=.false.

  type(fluid) :: flow
  type(body)  :: bodies
  integer     :: n(3)
  real        :: dt, t1, pforce(ndims), vforce(ndims), A
  logical     :: root
!
! -- Initialize simulation
  call init_mympi(ndims,set_blocks=b(:ndims),set_periodic=p(:ndims))
  root = mympi_rank()==0

  if(ndims==3) then
     n = composite((/6.0*DG,2.4*DG,16*D/), root)
     xg(3)%h = 0.5*DG/D
  else
     n = composite((/6.0*DG,2.4*DG,0./), root)
  end if
  call xg(1)%stretch(n(1), -10*DG, -0.6*DG, 0.6*DG, 20*DG, h_max=DG/D, prnt=root)
  call xg(2)%stretch(n(2), -10*DG, -0.6*DG, 0.6*DG, 10*DG, prnt=root)

  if(single) then
     bodies = cylinder(1,1,3,0.5*DG,center=0.)
  else
     bodies = make_array()
  end if

  A = DG*n(3)*xg(3)%h
  call flow%init(n/b, bodies, nu=nu)

  if(root) print *, '-- init complete --'
!
! -- Time update loop
  do while (flow%time/DG<finish)
     dt = flow%dt
     t1 = flow%time+dt
!
! -- smoothly increase flow speed
     flow%velocity%e(1)%bound_val = tanh(t1/(50*DG))
     flow%g(1) = sech(t1/(50*DG))**2/(50*DG)
!
! -- update flow
     call flow%update()
!
! -- print force coefficients
     pforce = -bodies%pforce(flow%pressure)
     vforce = nu*bodies%vforce(flow%velocity)
     write(9,'(f10.4,f8.4,4e16.8)') t1/DG,dt,2.*pforce(:2)/A,2.*vforce(:2)/A
     flush(9)
!
! -- print flow
     if(mod(t1,dtPrint*DG)<dt) call flow%write()
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
!    integer,parameter :: n(4) = (/6,13,19,25/)
    integer,parameter :: n(4) = (/6,12,18,26/) ! edit for symmetry
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
