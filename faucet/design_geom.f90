!-------------------------------------------------------!
!----------------- Jet Impact Test case ----------------!
!-------------------------------------------------------!
program design_geom
  use analytic
  implicit none
  integer,parameter :: ndims=3
  real(8),parameter :: r0 = 0.5, v0 = 1.0, x1 = 0.8
  type(set) :: pipe,water
  real (8)  :: velo(3),tip(3),norm(3)

  open(7,file='inp.geom')
  velo = 0; velo(ndims) = -v0
  tip  = 0; tip (ndims) =  x1
  norm = 0; norm (ndims) = 1.
!
! -- geom_body 
  write(7,'(a4)') 'body'
  pipe = .set.cylinder(1,ndims,r0,0,0,velo)-.set.plane(1,norm,tip,0,velo) &
     .or.(.set.cylinder(1,ndims,r0*1.1,0,0,0)-.set.cylinder(1,ndims,r0,0,0,0)-.set.plane(1,norm,tip,0,0))
  call set_write(7,pipe)

! -- geom_fint
  write(7,'(a4)') 'fint'
  water = .set.cylinder(1,ndims,r0,0,0,0)-.set.plane(1,norm,tip,0,velo) &
       .or..set.plane(1,norm,0,0,0)
  call set_write(7,water)

end program design_geom
