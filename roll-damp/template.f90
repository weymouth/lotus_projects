!-------------------------------------------------------!
!---------------- Rolling crossed plates ---------------!
!-------------------------------------------------------!
program crossed
  use fieldMod,   only: field
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  use imageMod,   only: display
  use geom_shape  ! to define geom (set,eps,plane, etc)
  implicit none
  ! physical parameters
  real,parameter     :: beta = 1e5            !a^2*omega/nu
  real,parameter     :: Theta = THETAVAL           !Max roll amplitude
  real,parameter     :: Per = 1               !# periods
  real,parameter     :: r = RVAL              ! central circle radius
  real,parameter     :: alpha = ALPHAVAL        ! plate crossing angle
  ! numerical parameters
  real,parameter     :: a = 140               !# cells across plate
  integer,parameter  :: ndims = 2             !# dimension
  real,parameter     :: m = 5*a               !approximate # points in x,y
  integer            :: bl(3) = (/4,4,1/)     !# blocks
  ! viz parameters
  real,parameter     :: dPrnt = 0.125         !how often to update image
  real,parameter     :: limit = 0.5           !vorticity contour limits
  integer,parameter  :: window(4) = int((/-2*a,-2*a,4*a,4*a/)) !image window
  ! resultant parameters
  real,parameter     :: omega = 1./a/Theta     !angular freq (radians)
  real,parameter     :: freq = omega/(2*pi)    !angular freq (cycles)
  real,parameter     :: nu = a**2*omega/beta   !viscosity
  ! variables
  integer            :: n(3),i                 !# points
  real               :: power,dt               !moment coefficient
  logical            :: root,there=.false.     !logical flags
  type(fluid)        :: flow                   !fluid
  type(body)         :: geom                   !solid
  type(field)        :: vort
!
! -- Initialize MPI (if MPI is ON)
#if MPION
  call init_mympi(ndims,set_blocks=bl(:ndims))
#else
  bl=1
#endif
  root = mympi_rank()==0
!
! -- Initialize array size
  n = composite((/m,m,32./),prnt=root)
  if(ndims==2) n(3) = 1
!
! -- Initialize and print grid
  call xg(1)%stretch(n(1),-10*a,-2*a,2*a,10*a,prnt=root)
  call xg(2)%stretch(n(2),-10*a,-2*a,2*a,10*a,prnt=root)
  if(ndims==3) xg(3)%h = 0.5*a/n(3)
!
! -- Initialize the geometry
  geom = make_geom().map.init_rigid(6,phi)
  if(root) print *,'omega=',omega,', nu=',nu
!
! -- Initialize fluid
  call flow%init(n/bl,geom,nu=nu)
!  flow%time = 0
  if(root) print *, '-- init complete --'
!
! -- Time update loop
  do while (flow%time*freq<Per.and..not.there)
     flow%dt = min(1.,flow%dt) ! limit time step (only needed initially)
     dt = flow%dt
     call geom%update(flow%time+flow%dt)
     call flow%update(geom)
!
! -- Print fields and force
     power = (nu*geom%vpower(flow%velocity)+geom%ppower(flow%pressure))/(omega**3*a**4)
     write(9,1) flow%time*freq,dt,power,phi(real(flow%time,8))
     flush(9)
     if(mod(flow%time,dPrnt/freq)<dt) then ! print
       vort = flow%velocity%vorticity_Z()
       call display(vort,'vort',lim=limit,box=window)
       if(ndims==3) call display(vort%average(),'vortAve',lim=limit,box=window)
       if(root) print 1,flow%time*freq,dt,power,phi(real(flow%time,8))
     end if
     inquire(file='.kill', exist=there)
1    format(f12.6,f8.4,2e14.6)
  end do
  call flow%write(geom)
  if(root) print *, '--- complete ---'
!
! -- Finalize MPI
#if MPION
  call mympi_end
#endif
contains
!
! -- motion definition
  real(8) pure function phi(t)
    real(8),intent(in) :: t
    phi = Theta*cos(omega*t)
  end function phi
!
! -- shape definition
  type(set) function make_geom() result(geom)
    geom = rect(-a,-3.5,2*a,7.,0.).or.rect(-a,-3.5,2*a,7.,alpha) &
            .or.cylinder(axis=3,radius=r*a,center=0)
  end function

  type(set) function rect(x,y,a,b,alpha)
    real,intent(in) :: x,y,a,b,alpha
    real :: sa,ca,xp,yp,xc,yc
    sa = sin(alpha); ca = cos(alpha)
    xc = ca*x-sa*y; yc = ca*y+sa*x
    xp = ca*(x+a)-sa*(y+b); yp = ca*(y+b)+sa*(x+a)
    rect = plane(norm=(/sa,-ca,0./),center=(/xc,yc,0./)) &
         .and.plane(norm=(/-ca,-sa,0./),center=(/xc,yc,0./)) &
         .and.plane(norm=(/-sa, ca,0./),center=(/xp,yp,0./)) &
         .and.plane(norm=(/ ca, sa,0./),center=(/xp,yp,0./))
  end function
end program crossed
