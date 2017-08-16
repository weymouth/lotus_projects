program fish
  use bodyMod,    only: body
  use fluidMod,   only: fluid
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  use imageMod,   only: display
  use geom_shape
  implicit none
!
! -- Physical parameters
  real,parameter     :: Re = 5e3, f = 2
!
! -- Numerical parameters
  real,parameter     :: c = 320, m(3) = (/2.,1.6,0./)
  integer            :: b(3) = (/4,4,1/), box(4) = (/-0.5*c,-0.4*c,3*c,0.8*c/)
!
! -- Variables
  integer :: n(3)
  logical :: root, there = .false.
  type(fluid) :: flow
  type(body) :: geom
!
! -- Initialize
  call init_mympi(2,set_blocks=b)
  root = mympi_rank()==0
  if(root) print *,'Setting up the grid, body and fluid'
  if(root) print *,'-----------------------------------'
  n = composite(c*m,prnt=root)
  call xg(1)%stretch(n(1), -6*c, -0.2*c, 1.2*c, 7*c, h_min=2., h_max=16.,prnt=root)
  call xg(2)%stretch(n(2), -2.4*c, -0.2*c, 0.2*c, 2.4*c, prnt=root)
  geom = naca(chord=c,pivot=0.).map.init_warp(2,h,doth,dh)
  call flow%init(n/b,geom,V=(/1.,0.,0./),nu=c/Re,exit=.true.)
  call display(flow%velocity%vorticity_Z(), 'out_vort', lim = 0.25, box=box)
!
  if(root) print *,'Starting time update loop'
  if(root) print *,'-----------------------------------'
  if(root) print *,' -t- , -dt- '
  do while(flow%time*f/c<15 .and..not.there)
    call geom%update(flow%time)
    call flow%update(geom)
    write(9,'(f10.4,f8.4,8e16.8)') flow%time*f/c,flow%dt, &
       2./Re*geom%vforce(flow%velocity), &
      -2./c*geom%pforce(flow%pressure), &
       2./Re*geom%vpower(flow%velocity), &
      -2./c*geom%ppower(flow%pressure)
    flush(9)
    if(mod(flow%time,0.125*c/f)<flow%dt) then
      if(root) print '(f7.3,",",f6.3)',flow%time*f/c,flow%dt
      call display(flow%velocity%vorticity_Z(), 'out_vort', lim = 0.25, box=box)
      call write_line()
    end if
    inquire(file='.kill', exist=there)
  end do
  if(root) print *,'Loop complete: writing restart files and exiting'
  if(root) print *,'-----------------------------------'
  call flow%write()
  call mympi_end
contains
!
! -- fully prescribed carangiform motion (Maertens 2016 eq A1)
  real pure function h(x)
    real,intent(in) :: x(3)
    h = env(x(1))*sin(arg(x(1)))
  end function h
  real pure function doth(x)
    real,intent(in) :: x(3)
    doth = env(x(1))*cos(arg(x(1)))*2*pi*(-f)/c
  end function doth
  pure function dh(x)
    real,intent(in) :: x(3)
    real            :: dh(3)
    dh = 0
    dh(1) = denv(x(1))*sin(arg(x(1))) &
          + env(x(1))*cos(arg(x(1)))*2*pi/c
  end function dh
  real pure function arg(x)
    real,intent(in) :: x
    arg = 2*pi*(x/c-f*flow%time/c)
  end function arg
  real pure function env(x)
    real,intent(in) :: x
    real :: xp
    xp = min(max(x/c,0.),1.)
    env = c*(0.1-0.0825*(xp-1)+0.1625*(xp**2-1))
  end function env
  real pure function denv(x)
    real,intent(in) :: x
    real :: xp
    xp = min(max(x/c,0.),1.)
    denv = -0.0825+2.*0.1625*xp
  end function denv

  type(set) function naca(chord,thick,pivot,alpha)
    real,intent(in) :: chord
    real,intent(in),optional :: thick,pivot,alpha
    type(model_info) :: info
    ! the geometry is a NACA0012 defined from x=2.85-5.58, so
    real :: thick0=0.12, edge=2.8538, chord0=2.7303
    real :: t=0.12,piv=0.5,a=0

    if(present(thick)) t = thick
    if(present(pivot)) piv = pivot
    if(present(alpha)) a = alpha*pi/180. !convert to rad

    info%file = 'naca_square.IGS'
    info%x = (/-edge-piv*chord0,-10.271,-18.87/)
    info%r = (/a,0.,0./)
    info%s = chord/chord0*(/1.,t/thick0,-1./)
    info%xmax(1) = chord
    info%n = (/chord,chord,1./)
    ! surface_debug = .true.
    model_fill = .false.
    eps = 2.0
    naca = model_init(info)
  end function naca

  subroutine write_line
      real,allocatable :: fp(:,:)
      real :: x
      integer :: i
      fp = -2./c*geom%pforce_plane(flow%pressure)
      do i=1,n(1)
        x = xg(1)%x(i)+xg(1)%dx(i)/2
        if(x>=-4 .and. x<=c+4) write(10,'(3e16.8)') flow%time*f/c,x/c,fp(i,1)
      end do
      flush(10)
  end subroutine write_line
end program fish
