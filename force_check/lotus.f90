program force_check
  use geom_shape, only: sphere,operator(.map.),init_rigid
  use gridMod,    only: xg
  use bodyMod,    only: body
  use fluidMod,   only: fluid
  implicit none
!
! -- Define parameters, declare variables
  real,parameter    :: D = 64 ! diameter
  integer           :: n(3) = D*(/2,2,2/),time,i
  type(fluid) :: flow
  type(body)  :: geom
!
! -- Initialize
  xg%left = -D  ! uniform grid
  ! do i=1,3 ! stretched grid
  !   call xg(i)%stretch(n(i),-2*D,-0.5*D-6,0.5*D+6,2*D,prnt=.true.)
  ! end do
  geom = sphere(radius=D/2, center=0) &  ! sphere
          .map.init_rigid(1,x).map.init_rigid(2,y).map.init_rigid(3,z)

  call flow%init(n,geom)  ! fluid
  call flow%pressure%eval(linear)
  write(9,'(i4,3e16.8)') 0,geom%pforce(flow%pressure)
!
! -- Time update loop
  do time = 1,50
    call geom%update(real(time))  ! apply mapping to geom
    call geom%measure(flow%mu0,flow%mu1,flow%ub)  ! remeasure
    write(9,'(i4,3e16.8)') time,geom%pforce(flow%pressure)
  end do

contains

  real(8) pure function x(t)
    real(8),intent(in) :: t
    x = 2*sin(t/3)
  end function x
  real(8) pure function y(t)
    real(8),intent(in) :: t
    y = 2*sin(t/5)
  end function y
  real(8) pure function z(t)
    real(8),intent(in) :: t
    z = 2*cos(t/7)
  end function z
  real pure function linear(x)
    real,intent(in) :: x(3)
    linear = sum(x)
  end function linear
end program force_check
