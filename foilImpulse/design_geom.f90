!-------------------------------------------------------!
!--------------- Foil Impulse Test case ----------------!
!-------------------------------------------------------!
program design_geom
  use analytic
  use transform
  implicit none
  type(set) :: body
  character(20) :: string
  character(16) :: name
  character(255) :: cmdln
  real(8) :: d,a
  integer :: i
!
! -- files
  open(7,file='inp.geom')
  open(7+1,file='srf.dat',STATUS='REPLACE')
  open(7+2,file='grd.dat',STATUS='REPLACE')
  open(7+3,file='xcp.dat',STATUS='REPLACE')
!
! -- defaults
  d = 0.5; a = 10 
!
! -- read command line input
  do i=1,iargc()
     call getarg(i,string)
     select case(string(2:2))
     case('n')
        read(string(4:),*) name
        print *,'file name',name
     case('a')
        read(string(4:),*) a
        print *,'angle of attack',a
     case('d')
        read(string(4:),*) d
        print *,'displacement time',d
     case default
        print *,string
        stop 'unknown argument'
     end select
  end do
!
! -- construct the stingray
  write(7,'(a4)') 'body'
  body = (.set.name.map.(init_affn()**(/a,0.,0./))).map.init_jerk(3,7.,7.+d,10./7.)
  call set_write(7,body)
!
! -- print commandline
  call get_command(cmdln)
  write(7,*) 'stingray geom: generated using command set:'
  write(7,*) trim(cmdln)

end program design_geom
