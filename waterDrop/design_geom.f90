!-------------------------------------------------------!
!-------------- Water Body Free Drop -------------------!
!-------------------------------------------------------!
program design_geom
  use analytic
  implicit none
  character(255) :: cmdln
  character(20) :: string
  type(set) :: water
  real(8)   :: dropHeight,waterHeight
  integer   :: i
!
! -- open file and set defaults
  open(7,file='inp.geom')
  dropHeight = 1.; waterHeight = 0.5
!
! -- read command line input
  do i=1,iargc()
     call getarg(i,string)
     select case(string(2:2))
     case('d')
        read(string(4:),*) dropHeight
        print *,'dropHeight',dropHeight
     case('w')
        read(string(4:),*) waterHeight
        print *,'waterHeight',waterHeight
     case default
        print *,string
        stop 'unknown argument'
     end select
  end do
!
! -- geom_fint
  write(7,'(a4)') 'fint'
  water = (.set.cylinder(1,3,0.5,0,0,0) &
       .and..set.plane(1,(/0,0,-1/),(/0.,0.,dropHeight/),0,0) &
       .and..set.plane(1,(/0,0, 1/),(/0.,0.,dropHeight+waterHeight/),0,0)) &
       .or..set.plane(1,(/0,0,1/),0,0,0)
  call set_write(7,water)
!
! -- print commandline
  call get_command(cmdln)
  write(7,*) 'waterDrop geom: generated using command set:'
  write(7,*) trim(cmdln)
end program design_geom
