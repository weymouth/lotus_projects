program foilTest
  use fluidMod, only: fluid
  use bodyMod,  only: body
  use gridMod,  only: xg
  use imageMod, only: display
  use mympiMod
  implicit none
  type(fluid)        :: flow
  type(body)         :: foil
  integer,parameter  :: ndims = 3
  real,parameter     :: c = 160
  real,parameter     :: Re = 50e3
  real,parameter     :: amp=0.1/c, freq = 2*3.14159*2
  real,parameter     :: tStop = 45, dtPrint = 1
  real,parameter     :: nu = c/Re
  logical,parameter  :: per(3) = (/.false.,.false.,.true./)
  integer            :: n(3) = (/384,256,64/)
  integer            :: b(3) = (/4,4,1/)
  real               :: mag = 0
  logical            :: root

  call init_mympi(ndims,set_blocks=b(1:ndims),set_periodic=per(1:ndims))
  root = mympi_rank()==0
  if(root) print *, '-- Foil test case --'

  call xg(1)%stretch(n(1),-2.5*c,-0.6*c,3*c,5*c,h_min=2.,h_max=10.,prnt=root)
  call xg(2)%stretch(n(2),-5*c,-0.4*c,0.4*c,5*c,prnt=root)
  if(ndims==3) then
     xg(3)%h = 4
  else
     n(3) = 1
  end if
  if(root) print *,'n(3),xg(3)%h',n(3),xg(3)%h

  foil = naca(c,-4.)
  call flow%init(n/b,foil,V=(/1.,0.,0./),nu=nu)
  if(root) print *, '-- init complete --'

  do while (flow%time/c<tStop)
     mag = amp*sin(freq*(flow%time+flow%dt)/c)
     call flow%ub%eval(velo)
     call flow%update
     write(9,'(f10.4,f8.4,3e16.8)') flow%time/c,flow%dt, &
          -2.*foil%pforce(flow%pressure)/(c*xg(3)%h*n(3))
     flush(9)
     if(mod(flow%time/c,dtPrint)<flow%dt/c) then
        call display(flow%velocity%vorticity_Z(),'vort',&
            lim=15./c, box=(/-160,-320,1280,640/))
        call flow%write(average=.true.)
        if(root) print *,flow%time/c,flow%dt
     end if
  end do
  if(root) print *,'--- complete --- '
contains
  type(set) function naca(chord,aoa)
    use geom_shape
    real,intent(in) :: chord,aoa
    type(model_info) :: info
    type(affn) :: shift
    info%file = 'naca_square.IGS'
    info%x = (/-4.219,-10.271,-18.87/)
    info%r = (/aoa,0.,0./)
    info%s = 0.36626*chord*(/1,1,-1/)
    info%xmax(1) = chord*1.5
    info%xmin(2) = -chord/2
    info%n = (/chord,chord,1./)
!    surface_debug = .true.
    eps = 2.0
    naca = model_init(info)
  end function naca

  pure function velo(x)  ! rotation velocity
    real(8),intent(in) :: x(3)
    real(8) :: velo(3)
    velo = (/-mag*x(2),mag*x(1),real(0.,8)/)
  end function velo
end program foilTest
