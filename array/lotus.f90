!-------------------------------------------------------!
!------------------- Cylinder Array  -------------------!
!-------------------------------------------------------!
program array
  use fluidMod, only: fluid
  use bodyMod,  only: body
  use mympiMod, only: init_mympi,mympi_end,mympi_rank
  use gridMod,  only: xg,composite
  use geom_shape  ! to define geom (set,eps,plane, etc)
  implicit none

  type(fluid) :: flow
  type(body)  :: bodies

  integer,parameter :: rows = 2, b(3) = (/1,1,1/), finish = 1000, dtPrint = 100
  real,parameter    :: D = 16, Re = 100, DG = 11*D, nu = D/Re
  integer :: n(3), uni
  real    :: dt, t1, pforce(2), vforce(2)
!
! -- Initialize simulation
  call init_mympi(2,set_blocks=b(:2))

  bodies = make_array()

  n = composite(DG*(/5,2,0/), mympi_rank()==0)
  uni = 0.55*DG
  call xg(1)%stretch(n(1), -12*DG, -uni, 3*uni, 30*DG, h_max=15., prnt=mympi_rank()==0)
  call xg(2)%stretch(n(2), -22*DG, -uni, uni, 22*DG, prnt=mympi_rank()==0)

  call flow%init(n/b, bodies, nu=nu)

  call flow%write(bodies)
  if(mympi_rank()==0) print *, '-- init complete --'
!
! -- Time update loop
  do while (flow%time/D<finish)
!
! -- update flow
     dt = flow%dt
     t1 = flow%time+dt
     flow%velocity%e(1)%bound_val = tanh(t1/(50*D))
     flow%g(1) = sech(t1/(50*D))**2/(50*D)

     call flow%update()
!
! -- print force
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
  real pure function sech(x) ! stable for large x
    real,intent(in) :: x
    sech = 2*exp(-x)/(1+exp(-2*x))
  end function sech

  type(set) function make_array()
    integer :: i,j,n(4) = (/6,13,19,25/)
    real :: theta,xc,yc,R=0.5*DG-0.5*D
    make_array = place_cyl(0.,0.).map.init_affn()
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
