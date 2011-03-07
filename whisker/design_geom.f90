!-------------------------------------------------------!
!----------------- Dummy Geom File ---------------------!
!-------------------------------------------------------!
program design_geom
  use analytic
  implicit none
  type(set) :: dummy
  dummy = .set.cylinder(1,3,1,0,0,0)
!
! -- write it up
  open(7,file='inp.ana')
  call set_write(7,dummy)  ! body
  call set_write(7,dummy)  ! free surface
  call set_write(7,dummy)  ! velocity
  call set_write(7,dummy)  ! domain geometry

end program design_geom
