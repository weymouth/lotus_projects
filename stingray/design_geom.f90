!-------------------------------------------------------!
!---------- Stingray with Rajiform Motion --------------!
!-------------------------------------------------------!
program design_geom
  use geom_shape
  implicit none
  type(model_info) :: tail,pos,neg,main
  character(20) :: string
  character(255) :: cmdln
  real(8) :: amp=0.2,f=1.33,k=1.475
  logical :: bio=.true.
  integer :: i
!
! -- files
  open(7,file='inp.geom')
  open(7+1,file='srf.dat',STATUS='REPLACE')
  open(7+2,file='xcp.dat',STATUS='REPLACE')
  model_fill = .false.
  surface_debug = .true.
!
! -- read command line input
  do i=1,iargc()
     call getarg(i,string)
     select case(string(2:2))
     case('a')
        read(string(4:),*) amp
        print *,'amplitude',amp
     case('f')
        read(string(4:),*) f
        print *,'frequency',f
     case('k')
        read(string(4:),*) k
        print *,'wavenumber',k
     case('b')
        read(string(4:),*) bio
        print *,'bio',bio
     case default
        print *,string
        stop 'unknown argument'
     end select
  end do
!
! -- model info
  tail%file = 'body2_tail.IGS'
  tail%r(2) = 90
  tail%s = 8.14E-3
  tail%buff = 0.05
!
  pos = tail
  pos%file = 'body2_part2.IGS'
!
  neg = tail
  neg%file = 'body2_part3.IGS'
!
  main = tail
  main%file = 'sim_body_part1.IGS'
  main%n = (/41,41,61/)
!
! -- construct the stingray
  write(7,'(a4)') 'body'
  call shape_write(7,((model_init(main) &
       .or.(model_init(pos).and.plane(norm=(/0, 1,0/),center=(/0., .633,0./))) &
       .or.(model_init(neg).and.plane(norm=(/0,-1,0/),center=(/0.,-.633,0./)))) &
       .map.raji(amp=amp,omega=2*acos(-1.)*f,k=2*k,bio=bio)) &
       .or.(model_init(tail).and.plane(norm=(/1,0,0/),center=(/3,0,0/))))
!
! -- print commandline
  call get_command(cmdln)
  write(7,*) 'stingray geom: generated using command set:'
  write(7,*) trim(cmdln)
end program design_geom
