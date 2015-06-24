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
  integer,parameter  :: ndims = 3
  real,parameter     :: c = 164
  real,parameter     :: alpha = -10
  real,parameter     :: Re = 5.3e3
  real,parameter     :: tStop = 80, dtPrint = 1
  real,parameter     :: nu = c/Re
  logical,parameter  :: per(3) = (/.false.,.false.,.true./)
  integer            :: n(3) = (/384,256,512/)
  integer            :: b(3) = (/4,4,1/)
  logical            :: root
#if MPION
  call init_mympi(ndims,set_blocks=b(1:ndims),set_periodic=per(1:ndims))
#else
  b=1
#endif
  root = mympi_rank()==0

  if(root) then
     print *, '-- Foil test case --'
     print '("    c/h=",i0,", nu=",f0.4,", del=",f0.4)', &
          int(c),nu,6.*sqrt(c/Re)
  end if

  call xg(1)%stretch(n(1),-2.5*c,-0.6*c,3*c,5*c,h_min=2.,h_max=10.,prnt=root)
  call xg(2)%stretch(n(2),-5*c,-0.3*c,0.6*c,5*c,prnt=root)
  if(ndims==3) then
     xg(3)%h = 2
  else
     n(3) = 1
  end if
  foil = naca(c,alpha)
  call flow%init(n/b,foil,V=(/1.,0.,0./),nu=nu)

  if(root) print *, '-- init complete --'
!
! -- Time update loop
  do while (flow%time/c<tStop)
     call flow%update
     write(9,'(f10.4,f8.4,3e16.8)') flow%time/c,flow%dt, &
          -2.*foil%pforce(flow%pressure)/(c*xg(3)%h*n(3))
     flush(9)
     if(mod(flow%time/c,dtPrint)<flow%dt/c) then
        call flow%write(lambda=.true.,average=.true.)
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
