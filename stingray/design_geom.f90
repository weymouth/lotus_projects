!-------------------------------------------------------!
!---------- Stingray with Rajiform Motion --------------!
!-------------------------------------------------------!
program design_geom
  use analytic
  implicit none
  type(set) :: body
  character(20) :: string
  character(255) :: cmdln
  real(8) :: amp,f,k,off,peak,tail
  integer :: i
!
! -- files
  open(7,file='inp.geom')
  open(7+1,file='srf.dat',STATUS='REPLACE')
  open(7+2,file='grd.dat',STATUS='REPLACE')
  open(7+3,file='xcp.dat',STATUS='REPLACE')
!
! -- defaults
  amp=0.2; f=1.33; k=2.2; off=0; peak=90; tail=0
!
! -- read command line input
  do i=1,iargc()
     call getarg(i,string)
     select case(string(2:2))
     case('a')
        read(string(4:),*) amp
        print *,'a',amp
     case('f')
        read(string(4:),*) f
        print *,'f',f
     case('k')
        read(string(4:),*) k
        print *,'k',k
     case('o')
        read(string(4:),*) off
        print *,'o',off
     case('p')
        read(string(4:),*) peak
        print *,'p',peak
     case('t')
        read(string(4:),*) tail
        print *,'t',tail
     case default
        print *,string
        stop 'unknown argument'
     end select
  end do
!
! -- construct the stingray
  write(7,'(a4)') 'body'
  body = .set.'sim_body_part1.IGS'.or.&
       (.set.'body2_part2.IGS'.and..set.plane(1,(/0, 1,0/),(/0., .633,0./),0,0)).or.&
       (.set.'body2_part3.IGS'.and..set.plane(1,(/0,-1,0/),(/0.,-.633,0./),0,0))
  body = body.map.init_raji(amp,f,k,peak,tail,off)
  body = body.or.(.set.'body2_tail.IGS'.and..set.plane(1,(/1,0,0/),(/3,0,0/),0,0))
  call set_write(7,body)
!
! -- print commandline
  call get_command(cmdln)
  write(7,*) 'stingray geom: generated using command set:'
  write(7,*) trim(cmdln)
end program design_geom
