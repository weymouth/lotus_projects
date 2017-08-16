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
  real,parameter     :: Re = 5e3, f = 1, m_star = .5, lam = 1, zeta = 0.005, a_star = 0.1
  logical,parameter  :: flowing = .true., clamped = .true., free = .true.
!
! -- Numerical parameters
  real,parameter     :: c = 128, m(3) = (/2.,1.6,0./), h_min = 2, mu_a = 0*c
  integer            :: b(3) = (/4,4,1/), box(4) = (/-0.5*c,-0.4*c,3*c,0.8*c/)
!
! -- resultant parameters
  real,parameter     :: amp = a_star*c, mu = m_star*c, damp = (mu+mu_a)*zeta
  real,parameter     :: omega = 2*pi*f/c, k=2*pi/(lam*c), EI=(mu+mu_a)*omega**2/k**4
  integer,parameter  :: s = c/h_min
!
! -- Variables
  real :: y(s)=0,y0(s)=0,y00(s)=0,doty(s)=0,ddoty(s)=0,U=0
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
  call xg(1)%stretch(n(1), -6*c, -0.2*c, 1.2*c, 7*c, h_min=h_min, h_max=16.,prnt=root)
  call xg(2)%stretch(n(2), -2.4*c, -0.2*c, 0.2*c, 2.4*c, prnt=root)
  geom = rect().map.init_warp(2,h,doth,dh)
  call flow%init(n/b,geom,V=(/U,0.,0./),nu=c/Re,exit=.true.)
  flow%dt = 0.5
  call display(flow%velocity%vorticity_Z(), 'out_vort', lim = 0.25, box=box)
!
  if(root) print *,'Starting time update loop'
  if(root) print *,'-----------------------------------'
  if(root) print *,' -t- , -dt- '
  do while(flow%time*f/c<15 .and..not.there)
    call update_line()
    if(flowing) then
      call geom%update(flow%time+flow%dt) ! just sets a flag
      U = tanh((flow%time+flow%dt)*f/c)
      call flow%update(geom,V=(/U,0.,0./))
      flow%dt = min(flow%dt,0.5)
      write(9,'(f10.4,f8.4,8e16.8)') flow%time*f/c,flow%dt, &
         2./Re*geom%vforce(flow%velocity), &
        -2./c*geom%pforce(flow%pressure), &
         2./Re*geom%vpower(flow%velocity), &
        -2./c*geom%ppower(flow%pressure)
      flush(9)
      if(mod(flow%time,0.125*c/f)<flow%dt) then
        if(root) print '(f7.3,",",f6.3)',flow%time*f/c,flow%dt
        call display(flow%velocity%vorticity_Z(), 'out_vort', lim = 0.25, box=box)
      end if
    else
      flow%time = flow%time+flow%dt
      if(mod(flow%time,c/f)<flow%dt.and.root) print *,flow%time*f/c
    end if
    inquire(file='.kill', exist=there)
  end do
  if(root) print *,'Loop complete: writing restart files and exiting'
  if(root) print *,'-----------------------------------'
  call flow%write()
  call mympi_end
contains
  type(set) function rect() result(geom)
    geom = plane(norm=(/-1,0,0/),center=(/0,-3,0/)) &
         .and.plane(norm=(/0,-1,0/),center=(/0,-3,0/)) &
         .and.plane(norm=(/1,0,0/),center=(/c,3.,0./)) &
         .and.plane(norm=(/0,1,0/),center=(/c,3.,0./))
  end function
!
! -- free centerline motion
  subroutine update_line
    real,allocatable :: fp(:,:)
    real :: yn(s),q(s),y_w, new(s)
    integer :: i,j

    !! get stress and leading edge position
    if(flowing) then
    fp = -geom%pforce_plane(flow%pressure)/h_min
      do i=1,s
        j = xg(1)%hash(int(h_min*(i-0.5)))
        q(i) = fp(j,1)+mu_a*ddoty(i)
      end do
    else
      q = 0
    end if
    y_w = amp*sin(-2*pi*f*(flow%time+flow%dt)/c)*U

    !! update variables
    yn = EulerBernoulli(y_w,q)
    ddoty = (2*yn-5*y+4*y0-y00)/flow%dt**2
    doty = (3*yn-4*y+y0)/(2.*flow%dt)
    y00 = y0; y0 = y; y = yn

    !! print
    if(mod(flow%time+flow%dt,0.125*c/f)<flow%dt) then
      write(10,1) (flow%time*f/c,h_min/c*(i-0.5),q(i),y(i),i=1,s)
1     format(2e12.4,2e16.8)
      flush(10)
    end if
  end subroutine update_line
!
! -- implicit dynamic Euler-Bernoulli beam equation
  function EulerBernoulli(y_w,q) result(yn)
    real,intent(in) :: y_w,q(s)
    real :: b(s),A(5,s),yn(s)
    !! form bending matrix
    A = spread((/1,-4,6,-4,1/),2,s); b = 0

    !! apply boundary conditions (coefficients from stencil.py)
    if(clamped) then
      A(:,1) = (/0.,0.,48.,-10.67,1.92/); b(1) = 39.25*y_w
      A(:,2) = (/0.,-1.,5.67,-3.96,1./); b(2) = 1.71*y_w
    else ! pinned
      A(:,1) = (/0.,0.,16.70,-8.35,1.67/); b(1) = 10.02*y_w
      A(:,2) = (/0.,-4.913,5.956,-3.991,1./); b(2) = -1.948*y_w
    end if
    if(free) then
      A(:,s-1) = (/1.,-3.828,4.656,-1.828,0./)
      A(:,s) = (/0.828,-1.656,0.828,0.,0./)
    end if
    A = A*EI/h_min**4; b = b*EI/h_min**4

    !! add inertia, damping, and forcing
    A(3,:) = A(3,:)+2*(mu+mu_a)/flow%dt**2
    b = b-(mu+mu_a)*(-5*y+4*y0-y00)/flow%dt**2-damp*doty+q

    !! solve
    yn = penta(s,A,b)
  end function EulerBernoulli
!
! -- Solve pentadiagonal linear system
  function penta(n,A,b) result(x)
    implicit none
    integer,intent(in) :: n
    real,intent(in)    :: A(5,n),b(n)
    real :: alp(2,-1:n),bet(-1:n+2),mu,gam,x(n)
    integer :: i
    alp = 0; bet = 0
    do i=1,n
      gam = A(2,i)-alp(1,i-2)*A(1,i)
      mu = A(3,i)-alp(2,i-2)*A(1,i)-alp(1,i-1)*gam
      alp(1,i) = (A(4,i)-alp(2,i-1)*gam)/mu
      alp(2,i) = A(5,i)/mu
      bet(i) = (b(i)-bet(i-2)*A(1,i)-bet(i-1)*gam)/mu
    end do
    do i=n,1,-1
      bet(i) = bet(i)-alp(1,i)*bet(i+1)-alp(2,i)*bet(i+2)
    end do
    x = bet(1:n)
  end function penta
!
! -- interpolate displacement and velocity
  real pure function h(x)
    real,intent(in) :: x(3)
    h = interp(y,x(1))
  end function h
  real pure function doth(x)
    real,intent(in) :: x(3)
    doth = interp(doty,x(1))
  end function doth
  pure function dh(x)
    real,intent(in) :: x(3)
    real            :: dh(3)
    dh = 0
    dh(1) = interp(y,x(1)+0.5)-interp(y,x(1)-0.5)
  end function dh
  real pure function interp(y,x)
    real,intent(in) :: y(:),x
    integer :: i
    real :: r
    if(x>c+4 .or. x<-4) then
      interp = 0
    else
      i = floor(x/h_min+0.5)
      i = min(max(1,i),s-1)
      r = x/h_min-(i-0.5)
      interp = y(i)*(1-r)+y(i+1)*r
    end if
  end function interp
end program fish
