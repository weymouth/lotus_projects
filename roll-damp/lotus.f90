!-------------------------------------------------------!
!---------------- Hull with Sharp Bilges ---------------!
!-------------------------------------------------------!
program square_cyl
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  use imageMod,   only: display
  use geom_shape  ! to define geom (set,eps,plane, etc)
  implicit none
  ! physical parameters
  real,parameter     :: STK = sqrt(2.4e5)     !Stokes No.= sqrt(B^2*omega/(2*nu))
  real,parameter     :: BoD = 0.28/0.112      !Aspect ratio
  real,parameter     :: phizero = 0.143       !Roll amplitude
  real,parameter     :: Per = 10              !# periods
  ! numerical parameters
  real,parameter     :: D = 140               !depth as a function of cell size h
  real,parameter     :: rho = 1               !density
  real,parameter     :: Um = 1                !max veloctiy
  integer,parameter  :: ndims = 2             !dimensions
  integer            :: bl(3) = (/4,4,1/)     !blocks
  ! resultant parameters
  real,parameter     :: B = BoD*D             !beam in cells
  real,parameter     :: Omegam = Um*(2/B)     !angular velocity
  real,parameter     :: omega = Omegam/phizero!angular freq (rad)
  real,parameter     :: freq = omega/(2*pi)   !angular freq (Hz)
  real,parameter     :: nu = (B/STK)**2*omega/2.  ! viscosity
  real               :: m = 3*B               !approximate # points
  integer            :: window(4) = int((/-1.5*B,-3*D,3*B,6*D/))
  ! variables
  integer            :: n(3)                  !# points
  real               :: moment(3)             !moment coefficient
  logical            :: root,there=.false.     !logical flags
  type(fluid)        :: flow                  !fluid
  type(body)         :: square                !solid
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
  n = composite((/m,m,128./),prnt=root)
  if(ndims==2) n(3) = 1
!
! -- Initialize and print grid
  call xg(1)%stretch(n(1),-10*B,-B,B,10*B,prnt=root)
  call xg(2)%stretch(n(2),-10*B,-B,B,10*B,prnt=root)
!
! -- Initialize the square geometry
  square = make_geom().map.init_rigid(6,phi)
  if(root) print *,'omega=',omega,', nu=',nu
!
! -- Initialize fluid
  call flow%init(n/bl,square,nu=nu)
  if(root) print *, '-- init complete --'
!
! -- Time update loop
  do while (flow%time*freq<Per.and..not.there)
     flow%dt = min(0.5,flow%dt) ! limit time step (only needed initially)
     call square%update(flow%time+flow%dt)
     call flow%update(square,f=SGS(flow%velocity,C=0.2))
!
! -- Print fields and force
     moment = -square%pmoment(flow%pressure)/(rho*phizero*Omega**2*B**3*D)
     write(9,1) flow%time*freq,flow%dt,moment(3),phi(real(flow%time,8))
     flush(9)
     if(mod(flow%time,0.125/freq)<flow%dt) then ! print every 1/8 cycle
       call display(flow%velocity%vorticity_Z(),'vort',lim=0.5,box=window)
       if(root) print 1,flow%time*freq,flow%dt,moment(3),phi(real(flow%time,8))
     end if
     inquire(file='.kill', exist=there)
1    format(f10.4,f8.4,2e14.6)
  end do
  call display(flow%velocity%vorticity_Z(),'vort',lim=0.5,box=window)
  call flow%write(square)
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
    phi = phizero*cos(omega*t)
  end function phi
!
! -- shape definition
  type(set) function make_geom() result(geom)
    geom = rect(-0.5*B,0.,B,2*D,0.) !&
!          .or.bilges(0.5*B-0.1*D,D-0.1*D,0.1*B,2.,0.5*pi)
  end function

  type(set) function bilges(x,y,a,b,alpha)
    real,intent(in) :: x,y,a,b,alpha
    bilges = rect(x,y,a,b,alpha) &
        .or.rect(-x,y,a,b,pi-alpha) &
        .or.rect(x,-y,a,b,-alpha) &
        .or.rect(-x,-y,a,b,pi+alpha)
  end function bilges

  type(set) function rect(x,y,a,b,alpha)
    real,intent(in) :: x,y,a,b,alpha
    real :: sa,ca,xp,yp,xc,yc
    sa = sin(alpha); ca = cos(alpha)
    xc = x+sa*b/2.; yc = y-ca*b/2.
    xp = x+ca*a-sa*b/2.; yp = y+sa*a+ca*b/2.
    rect = plane(norm=(/sa,-ca,0./),center=(/xc,yc,0./)) &
         .and.plane(norm=(/-ca,-sa,0./),center=(/xc,yc,0./)) &
         .and.plane(norm=(/-sa, ca,0./),center=(/xp,yp,0./)) &
         .and.plane(norm=(/ ca, sa,0./),center=(/xp,yp,0./))
  end function
!
! -- LES SGS model
  type(vfield) function SGS(u,C) result(force)
    use gridMod, only: dxi,dxi2
    use fieldMod, only: field
    use vectorMod, only: vfield
    implicit none
    type(vfield),intent(in) :: u
    real,intent(in)         :: C !! literature suggests 0.1-0.25
    type(vfield) :: S(u%ndims)
    real         :: g(3,3),nu_t,tau(3,3)
    integer      :: ndims,is,ie,js,je,ks,ke,d1,d2,i,j,k,o(3),o2(3),ijk(3),id

    call u%e(1)%limits(is,ie,js,je,ks,ke)
    call force%init(n=u%e(1)%size())
    ndims = force%ndims
!
! -- get SGS stress at each point
    call u%gradTensor(S)
    do concurrent (i=is:ie,j=js:je,k=ks:ke)
      g = 0
      do concurrent(d1=1:ndims,d2=1:ndims)
        g(d1,d2) = S(d1)%e(d2)%p(i,j,k)
      end do
      nu_t = C**2*smagorinsky(g)
      do concurrent(d1=1:ndims,d2=1:ndims)
        S(d1)%e(d2)%p(i,j,k) = nu_t*(g(d1,d2)+g(d2,d1))
      end do
      ! tau = 0.5*(g-transpose(g))
      ! tau = C**2*g**2
      ! do concurrent(d1=1:ndims,d2=1:ndims)
      !   S(d1)%e(d2)%p(i,j,k) = tau(d1,d2)
      ! end do
    end do
    do d1=1,ndims
      call S(d1)%applyBC
    end do
!
! -- compute SGS force
    components: do d1=1,ndims
      faces: do d2=1,ndims
        o = 0; o(d1) = 1
        o2 = 0; o2(d2) = 1
        if(d1==d2) then ! in-line
           do concurrent (i=is:ie,j=js:je,k=ks:ke)
              ijk = (/i,j,k/); id=ijk(d1)
              force%e(d1)%p(i,j,k) = force%e(d1)%p(i,j,k)+dxi2(d1,id)* &
                  (S(d1)%e(d1)%p(i,j,k)-S(d1)%e(d1)%p(i-o(1),j-o(2),k-o(3)))
           end do
        else ! cross term
           do concurrent (i=is:ie,j=js:je,k=ks:ke)
              ijk = (/i,j,k/); id=ijk(d2)
              force%e(d1)%p(i,j,k) = force%e(d1)%p(i,j,k)+0.25*dxi(d2,id)* &
                  (S(d1)%e(d2)%p(i+o2(1),j+o2(2),k+o2(3)) &
                  +S(d1)%e(d2)%p(i+o2(1)-o(1),j+o2(2)-o(2),k+o2(3)-o(3)) &
                  -S(d1)%e(d2)%p(i-o2(1),j-o2(2),k-o2(3)) &
                  -S(d1)%e(d2)%p(i-o2(1)-o(1),j-o2(2)-o(2),k-o2(3)-o(3)))
           end do
        end if
       end do faces
    end do components
  end function SGS
  real pure function smagorinsky(g) result(nu)
    real,dimension(3,3),intent(in) :: g
    real,dimension(3,3) :: S
    S = 0.5*(g+transpose(g))
    nu = sqrt(sum(2.*S**2))
  end function smagorinsky
end program square_cyl
