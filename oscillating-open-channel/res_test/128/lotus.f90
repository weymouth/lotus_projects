!-------------------------------------------------------!
!---------------- Turbulent Channel Flow ---------------!
!-------------------------------------------------------!
program channel
    use fluidMod,   only: fluid
    use bodyMod,    only: body
    use gridMod,    only: xg,composite
    use mympiMod,   only: init_mympi,mympi_end,mympi_rank
    use imageMod,   only: display
    use geom_shape
    implicit none
  ! - Domain
    integer            :: b(3) = [2,2,2], n(3) ! blocks
    real               :: m(3) = [1*1./8, 2., 1*1./8]
  ! - Length scales
    integer,parameter  :: S = 32  ! Stokes length
    real,parameter     :: nu = S/100000
    real               :: finish=2*pi*10*S, print_int=2*pi*0.02*S, init_time=2*pi*3*S
    real               :: t,dt
    real               :: h_roughness=0.1  ! Set roughness height to same as stokes length
    real               :: thickness
    real               :: lambda=S*1./2   ! Roughness wavelength
  ! - Init sim
    type(fluid)        :: flow
    type(body)         :: wall
    logical            :: root,there=.false.
  !
  ! - io
    real               :: enstrophy, tke
  !
  ! -- Initialize MPI
    call init_mympi(3,set_blocks=b,set_periodic=[.true.,.false.,.true.])
    root = mympi_rank()==0
    !
    ! -- Init grid
    if(root) print *,'Setting up the grid, body and fluid'
    if(root) print *,'-----------------------------------'
    n = composite(S*m,prnt=root)
    xg(1)%h = 8.0
    call xg(2)%stretch(n(2),0.,0.,1.*S,3.0*S,prnt=root)
    xg(3)%h = 8.0
  !
  ! -- Init channel
    wall = upper(2.*h_roughness*S+2).map.init_rigid(1,position,velocity)
    eps=1.0
  !
  ! -- Initialize fluid
    call flow%init(n/b,wall,nu=nu)
    flow%dt = min(flow%dt,1.)
  !
  ! -- Time update loop
    if(root) print *,'Starting time update loop'
    if(root) print *,'-----------------------------------'
  !
    do while (flow%time<finish.and..not.there)
      t = flow%time
      dt = flow%dt
      call wall%update(t+dt)
      call flow%update(wall)
      if(mod(t,print_int)<dt) then
        if(root) print "('Time:',f15.3)",t/(2*pi*S)
        if(t>init_time) call flow%write(wall, write_vtr=.true.)
      end if
    ! - Write some stuff
      tke = flow%velocity%tke(mean=.true.)
      enstrophy = flow%velocity%enstrophy(mean=.true.)
      write(9,'(f10.4,f8.4,4e16.8)')&
        t/(2*pi*S),dt/t/(2*pi*S), enstrophy, tke
    end do
    !
    if(root) print *,'Loop complete: writing restart files and exiting'
    if(root) print *,'-----------------------------------'
    call mympi_end()
  contains
  !
  ! -- Kinematics
  real(8) pure function position(t)
    real(8),intent(in) :: t
    position = S*sin(t/S)
  end function position
  real(8) pure function velocity(t)
    real(8),intent(in) :: t
    velocity = cos(t/S)
  end function velocity
  !
  type(set) function upper(thickness) result(geom)
    real,intent(in) :: thickness
    geom = plane([0.,1.,0.],[0.,thickness,0.]).map.init_warp(2,egg_top,dotegg_top,degg_top)
  end function
  !
  real pure function egg_top(x)
    real,intent(in) :: x(3)
    egg_top = (h_roughness*S)*sin((2*pi*x(1)/lambda)-pi/2)*cos(2*pi*x(3)/lambda)
  end function egg_top
  pure function degg_top(x)
    real,intent(in) :: x(3)
    real            :: degg_top(3)
    degg_top = 0
    degg_top(1) = (h_roughness*S)*cos((2*pi*x(1)/lambda)-pi/2)&
                                  *sin((2*pi*x(3)/lambda))*(2*pi/lambda)
    degg_top(3) = -(h_roughness*S)*cos((2*pi*x(1)/lambda-pi/2))&
                                  *sin((2*pi*x(3)/lambda))*(2*pi/lambda)
  end function degg_top
  real pure function dotegg_top(x)
    real,intent(in) :: x(3)
    dotegg_top = 0
  end function dotegg_top
  end program channel
  