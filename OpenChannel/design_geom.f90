!-------------------------------------------------------!
!-------------------- Stern test case ------------------!
!-------------------------------------------------------!
program design_geom
  use analytic
  implicit none
  character(20) :: string
  integer :: ndims=3
  real(8) :: corner(3),BoD,norm(4,3),corner2(3)
  type(set) :: stern,water
!
! -- init
  call getarg(1,string)
  read(string,*) BoD
  open(7,file='inp.geom')
!
! -- geom_body 
  write(7,'(a4)') 'body'
  corner    = (/0.,BoD,-1./)
  corner2    = (/0.,-BoD,-1./)
  norm(1,:) = (/ 1.,0.,0./)  ! face forward
  norm(2,:) = (/0., 1.,0./)  ! face left
  norm(3,:) = (/0.,-1.,0./)  ! face right
  norm(4,:) = (/0.,0.,-1./)  ! face down
  stern = .set.plane(-1,norm(1,:),corner ,0,0).and.&
          .set.plane(-1,norm(2,:),corner ,0,0).and.&
          .set.plane(-1,norm(3,:),corner2,0,0).and.&
          .set.plane(-1,norm(4,:),corner ,0,0)
  call set_write(7,stern)
!
! -- geom_fint
  write(7,'(a4)') 'fint'
  water = .set.plane(1,(/0,0,1/),0,0,0)
  call set_write(7,water)

end program design_geom


