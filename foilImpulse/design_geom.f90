!-------------------------------------------------------!
!----------------- Jet Impact Test case ----------------!
!-------------------------------------------------------!
program design_geom
  use analytic
  implicit none
  integer,parameter :: ndims=3
  real(8),parameter :: r0 = 0.5, v0 = 1.0, x1 = 2.0
  type(set) :: pipe,water,falling_water,drain
  real (8)  :: velo(3),tip(3),norm(3)
  velo = 0; velo(ndims) = -v0
  tip  = 0; tip (ndims) =  x1
  norm = 0; norm (ndims) = 1.
!
! -- geom_body 
  pipe = .set.cylinder(-1,ndims,r0,0,0,velo)-.set.plane(1,norm,tip,0,velo)
!
! -- geom_fint
  water = .set.cylinder(1,ndims,r0,0,0,0).or..set.plane(1,norm,0,0,0)
!
! -- geom_velo
  falling_water = .set.cylinder(1,ndims,r0,0,0,velo)
!
! -- geom_domain
  drain = .set.cylinder(1,ndims,r0,0,0,velo)
!
! -- write it up
  open(7,file='geom.txt')
  call set_write(7,pipe)
  call set_write(7,water)
  call set_write(7,falling_water)
  call set_write(7,drain)

end program design_geom
