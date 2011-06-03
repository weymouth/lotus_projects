!-------------------------------------------------------!
!---------- JUST a screen with a few bubbles -----------!
!-------------------------------------------------------!
program design_geom
  use analytic
  implicit none
  type(set) :: dummy,scrn
  character(8) :: string
  real(8) :: h

  call getarg(1,string)
  read(string,*) h

  h = h*0.5 ! this is half a grid space. it is needed to get the
                       !   SGS screen to line up with the center of the x-faces
!
! -- write it up
  scrn = .set.cylinder(1,2,0.025D0,(/0.5D0,0.25D0,0.25D0+h/),0,0).or. &
         .set.cylinder(1,3,0.025D0,(/0.5D0,0.25D0+h,0.25D0/),0,0).or. &
         .set.cylinder(1,2,0.025D0,(/0.5D0,0.50D0,0.50D0+h/),0,0).or. &
         .set.cylinder(1,3,0.025D0,(/0.5D0,0.50D0+h,0.50D0/),0,0).or. &
         .set.cylinder(1,2,0.025D0,(/0.5D0,0.75D0,0.75D0+h/),0,0).or. &
         .set.cylinder(1,3,0.025D0,(/0.5D0,0.75D0+h,0.75D0/),0,0)
  open(7,file='inp.geom')
  write(7,'(a4)') 'body'
  call set_write(7,scrn)  ! body
  dummy = .set.sphere(1,0.1D0,(/1.0D0,0.5D0,0.5D0/),0,0).or. &
          .set.sphere(1,0.1D0,(/1.5D0,0.5D0,0.5D0/),0,0).or. &
          .set.sphere(1,0.1D0,(/2.0D0,0.5D0,0.5D0/),0,0).or. &
          .set.sphere(1,0.1D0,(/2.5D0,0.5D0,0.5D0/),0,0).or. &
          .set.sphere(1,0.1D0,(/3.0D0,0.5D0,0.5D0/),0,0)
  dummy = .not.dummy
  write(7,'(a4)') 'fint'
  call set_write(7,dummy)  ! free surface
  write(7,'(a4)') 'velo'
  call set_write(7,dummy)  ! velocity

end program design_geom
