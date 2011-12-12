!-------------------------------------------------------!
!------------------ Whisker motion ---------------------!
!-------------------------------------------------------!
program design_geom
  use analytic
  use path
  implicit none
  type(pathSet) :: motion(3)
  character(20) :: string
  real(8) :: a,f
  integer :: n
!
! -- read command line input
  call getarg(1,string)

  read(string,*,IOSTAT=n) a
  if(n<0.or.a.eq.0.D0) then
     print *, 'not making path file'
  else
     call getarg(2,string)
     read(string,*,IOSTAT=n) f
     if(n<0) stop 'missing f'
     print *,'a,f = ',a,f
!
! -- write the motion
     motion(1) = .set.line(0,0)
     motion(2) = .set.trig(a,f*2*acos(-1.D0),0.D0)
     motion(3) = .set.line(0,0)
     open(7,file='inp.mot')
     call pathSet_write(7,motion(1))  ! x
     call pathSet_write(7,motion(2))  ! y
     call pathSet_write(7,motion(3))  ! z
     close(7)
  end if
!
! -- write the geom
  open(7,file='inp.geom')
  write(7,'(a4)') 'body'
  call set_write(7,.set.(/'inp.IGS ','inp.IGS2'/))
  close(7)

end program design_geom
