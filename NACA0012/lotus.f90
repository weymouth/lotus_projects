!-------------------------------------------------------!
!--------------- Canonical foil tests ------------------!
!-------------------------------------------------------!
program foilTest
  use fluidMod, only: fluid
  use bodyMod,  only: body
  use gridMod,  only: xg
  use mympiMod
  implicit none
  type(fluid)        :: flow
  type(body)         :: foil
  integer,parameter  :: n(3) = (/1024,512,1/)
  real,parameter     :: c = 325
  real,parameter     :: alpha = 10
  real,parameter     :: Re = 5.3e3
  real,parameter     :: tStop=30,dtPrint=1
  real,parameter     :: nu = c/Re
  integer            :: b(3) = (/4,4,1/)
  logical            :: root
#if MPION
  call init_mympi(2,set_blocks=b(1:2))
#else
  b=1
#endif
  root = mympi_rank()==0
  if(root) then
     print *, '-- Foil test case --'
     print '("    c/h=",i0,", nu=",f0.4,", y+=",f0.4)', int(c),nu,sqrt(0.026/Re**(1./7.)/2.)/nu
  end if

  call xg(1)%stretch(n(1),-2.5*c,-0.6*c,1.4*c,5*c,h_max=4.,prnt=root)
  call xg(2)%stretch(n(2),-5*c,-0.6*c,0.3*c,5*c,prnt=root)

  foil = naca(c,alpha)
  call flow%init(n/b,foil,V=(/1.,0./),nu=nu)
  flow%velocity%e%exit = .true.
  call flow%write(foil)
!  flow%body = .false. ! turns off O(2)

  if(root) print *, '-- init complete --'

!
! -- Time update loop
  do while (flow%time/c<tStop)
     call flow%update
     write(9,'(f10.4,f8.4,3e16.8)') flow%time/c,flow%dt,2.*foil%pforce(flow%pressure)/c
     flush(9)
     if(mod(flow%time/c,dtPrint)<flow%dt/c) then
        call flow%write
        if(root) print *,flow%time/c,flow%dt
     end if
  end do

  if(root) print *,'--- complete --- '
  call mympi_end
contains
  type(set) function naca(chord,alpha)
    use geom_shape
    real,intent(in) :: chord,alpha
    type(model_info) :: info
    type(affn) :: shift
    info%file = 'naca_square.IGS'
    info%x = (/-4.219,-10.271,-18.87/)
    info%r = (/alpha,0.,0./)
    info%s = 0.36626*chord*(/1,1,-1/)
    info%xmax(1) = chord
    info%n = (/chord,chord,1./)
!    surface_debug = .true.
    eps = 2.0
    naca = model_init(info)
  end function naca
end program foilTest
