!-------------------------------------------------------!
!--------------- Foil Impulse Test case ----------------!
!-------------------------------------------------------!
program design_geom
  use geom_shape
  implicit none
  class(shape),pointer :: body
  type(model_info) :: info
  character(30) :: string
  character(26) :: name
  character(255) :: cmdln
  real(8) :: a=10,z=0,r=0
  integer :: i
!
! -- default
  name = 'naca_square'
  info%file = 'naca_square.IGS'
  info%x = (/-4.219,-10.271,0./)
  info%s = 0.36626
  info%n = (/31,31,141/)
!
! -- read command line input
  do i=1,iargc()
     call getarg(i,string)
     select case(string(2:2))
     case('n')
        read(string(4:),*) name
        print *,'file name:',name
        info%file = trim(name)//'.IGS'
     case('a')
        read(string(4:),*) a
        print *,'angle of attack: ',a
     case('z')
        read(string(4:),*) z
        print *,'z offset: ',z
        info%x(3) = z
     case('r')
        read(string(4:),*) r
        print *,'initial rotation: ',r
        info%r(2) = r
     case default
        print *,string
        stop 'unknown argument'
     end select
  end do
!
! -- files
  surface_debug = .true.
  model_fill = .false.
  open(7,file=trim(name)//'.geom')
  open(7+1,file='srf.dat',STATUS='REPLACE')
  open(7+2,file='xcp.dat',STATUS='REPLACE')
!
! -- construct the foil
  write(7,'(a4)') 'body'
  call shape_write(7,model_init(info) &
       .map.(init_affn()**(/a,0.,0./)) &
       .map.jerk(axis=3,t0=7.12,dt=0.5,dis=10./7.))
!
! -- print commandline
  call get_command(cmdln)
  write(7,*) 'foilImpluse geom: generated using command set:'
  write(7,*) trim(cmdln)

end program design_geom
