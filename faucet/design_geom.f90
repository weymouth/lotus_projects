!-------------------------------------------------------!
!----------------- Jet Impact Test case ----------------!
!-------------------------------------------------------!
program design_geom
  use analytic
  implicit none
  integer,parameter :: ndims=3
  real(8),parameter :: r0 = 0.5, v0 = 1.0, x1 = 4.0
  type(set) :: pipe,water,falling_water,drain,scrn,totalpipe,dummy
  character(8) :: string
  real (8)  :: velo(3),tip(3),norm(3),h,off(3)

  open(7,file='inp.geom')
  velo = 0; velo(ndims) = -v0
  tip  = 0; tip (ndims) =  x1
  norm = 0; norm (ndims) = 1.
  off = 0; off(ndims) = 0.25

  call getarg(1,string)
  read(string,*) h

!
! -- geom_body 
  write(7,'(a4)') 'body'
  pipe = .set.cylinder(1,ndims,r0,0,0,velo)-.set.plane(1,norm,tip+0.5,0,velo) &
     .or.(.set.cylinder(1,ndims,r0*1.1,0,0,0)-.set.cylinder(1,ndims,r0,0,0,0)-.set.plane(1,norm,tip,0,0))! &
     !.or..set.sphere(1,0.05,tip+(/0.1,0.4,0.25/),0,0) &
     !.or..set.sphere(1,0.05,tip+(/0.1,0.1,0.25/),0,0) &
     !.or..set.sphere(1,0.05,tip+(/0.4,0.1,0.25/),0,0)
  scrn = .set.sphere(1,h,tip+off+(/0.0D0,0.0D0,0.D0/),0,0).or. &
         .set.sphere(1,h,tip+off+(/0.125D0,0.125D0,0.D0/),0,0).or. &
         .set.sphere(1,h,tip+off+(/0.25D0,0.0D0,0.D0/),0,0).or. &
         .set.sphere(1,h,tip+off+(/0.25D0,0.25D0,0.D0/),0,0).or. &
         .set.sphere(1,h,tip+off+(/0.0D0,0.25D0,0.D0/),0,0).or. &
         .set.sphere(1,h,tip+off+(/0.3750D0,0.1250D0,0.D0/),0,0) .or. &
         .set.sphere(1,h,tip+off+(/0.125D0,0.375D0,0.D0/),0,0) 
  totalpipe=pipe.or.scrn
  !totalpipe=scrn

  call set_write(7,totalpipe)

! -- geom_fint
  write(7,'(a4)') 'fint'
  water = .set.cylinder(1,ndims,r0,0,0,0).or..set.plane(1,norm,0,0,0)
  call set_write(7,water)
!
! -- geom_velo
  write(7,'(a4)') 'velo'
  falling_water = .set.cylinder(1,ndims,r0,0,0,velo)-.set.plane(1,norm,-tip,0,velo)
  call set_write(7,falling_water)

end program design_geom
