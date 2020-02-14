program squeeze
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  use geom_shape
  use bodyMod,    only: body
  use fluidMod,   only: fluid
  use gridMod,    only: xg
  use imageMod,   only: display
  implicit none
!
! -- Define parameters, declare variables
  real,parameter    :: L = 170         ! major semi-axis size in cells
  real,parameter    :: beta0 = 0.25	   ! aspect ratio
  real,parameter    :: dm = DM_IN      ! fraction of mass expelled
  real,parameter    :: A0_Ae = AR_IN 	 ! area ratio
  real,parameter    :: per = PER_IN    ! periods of motion
  real,parameter    :: thk = 5         ! membrane half thickness in cells
  logical,parameter :: free = F_IN     ! Is the body free to surge?
  real              :: V = V_IN, a = a_IN   ! Background flow speed,accelation
  real,parameter    :: Re = 1e5        ! Approx reynolds number
						   ! Uj=10x(2L)/s with 2L = 0.05m and nu = 1*10^-6 m^2/s

  real,parameter    :: beta1 = dm*beta0/4       ! pulse amplitude
  real,parameter    :: xe = sqrt(1-1/A0_Ae)     ! exit location
  real,parameter    :: Vf = 2./3.+xe-xe**3/3.   ! volume factor V= Vf pi L**3
  real,parameter    :: T = Vf*dm/2.*A0_Ae*L     ! size-change timescale
  real,parameter    :: Af = pi*(beta0*L+thk)**2/4. ! mean frontal area
  real,parameter    :: A0 = pi*(beta0*L)**2/4. ! mean inner cross-section
  real,parameter    :: m =  2*Vf*pi*beta0**2*((L+thk/beta0)**3-L**3) ! mass
  real,parameter    :: nu = 2*L/Re              ! kinematic viscosity

  real,parameter    :: dprnt = 2*pi*T/DP_IN        ! how often to print?
  real,parameter    :: Tend  = per*2*pi*T          ! when should we stop?
  real,parameter    :: h_min = 3                    ! stretching in x
  real,parameter    :: f(3)  = (/2.2,0.75,0.75/)   ! approx grid factor
  integer           :: b(3)  = (/4,2,2/)           ! MPI domain cuts in ijk
  integer           :: n(3)                        ! number of cells in ijk
  integer           :: box(4) = (/-1.5*L,0.,5*L,L/) ! view for display
  logical           :: root, there = .FALSE., on = .FALSE.
  real              :: force(3)=0,pow,dt,ma
  real, allocatable :: u0(:,:)
  type(fluid)       :: flow
  type(body)        :: geom
!
! -- Initialize mpi, grid, body and fluid
  call init_mympi(ndims=3,set_blocks=b)
  root = mympi_rank()==0

  if(root) print *,'T',T,'xe',xe,'dm',dm,'Vf',Vf,'m',m,'nu',nu

  n = composite(L*f,prnt=root)                             ! n
  allocate(u0(1,n(1))); u0=0
  call xg(1)%stretch(n(1),-4.*L,-1.2*L,2*L,6.*L,h_min=h_min,h_max=6.,prnt=root) ! x
  call xg(2)%stretch(n(2),0.,0.,0.4*L,2*L,prnt=root)         ! y
  call xg(3)%stretch(n(3),0.,0.,0.4*L,2*L,prnt=root)         ! z

  geom = (sphere(radius=L+thk/beta0, center=0)-sphere(radius=L, center=0) &
         -(cylinder(axis=1, radius=L/sqrt(A0_Ae), center=0) &
         .and.plane(norm=(/-1,0,0/), center=0))) &
         .map.init_scale(2,beta,dbeta).map.init_scale(3,beta,dbeta)
  geom%dis_wall(2:3) = .true. ! ok to adjust velocity on +y,+z planes
  call flow%init(n/b, geom, V=(/V,0.,0./), nu=nu)
  flow%time = 0
  call flow%write(geom,lambda=.TRUE.)
  call line()

! -- Time update loop
  do while(flow%time<Tend .and. .not.there)

! -- update geom
    flow%dt = min(flow%dt,2.)
    dt = flow%dt
    if(on) then ! update velocity
      ma = beta(real(flow%time,8))
      ma = Vf*pi*L**3*ma**2+2./3.*pi*(L*ma)**3
      a = (a*ma+4*force(1))/(m+ma)
      V = V-dt*a
    else if(free) then
      on = .TRUE.
    end if
    call geom%update(flow%time+dt)

! -- update and measure flow
    call flow%update(geom,V=(/V,0.,0./))
    force = nu*geom%vforce(flow%velocity)-geom%pforce(flow%pressure) ! force
    pow = -nu*geom%vpower(flow%velocity)+geom%ppower(flow%pressure)  ! power

! -- write to file
    if(root) write(9,'(f10.4,f8.4,5e16.8)') flow%time/T,flow%dt,force/(0.5*A0),pow/(0.5*A0),V
    if(root) flush(9)
    if(mod(abs(flow%time),dprnt)<dt) then
      if(root) print '(f10.4,4f8.4)',flow%time/T,flow%dt, &
            beta(real(flow%time,8)),dbeta(real(flow%time,8))*T,V
      call display(flow%velocity%vorticity_Z(),name='01_out',lim=0.25)
      call flow%write(geom,lambda=.TRUE.)
    end if
    call line()
    inquire(file='.kill', exist=there)
  end do
  call mympi_end
contains
  real(8) pure function betao(t1)   ! pre-scale outer membrane
    real(8),intent(in) :: t1
    betao=(L+thk)/(L+thk/beta0)
  end function betao
  real(8) pure function betai(t1)  ! pre-scale inner membrane
     real(8),intent(in) :: t1
     betai=(L-thk)/(L-thk/beta0)
  end function betai
  real(8) pure function zero(t1)	   ! static
     real(8),intent(in) :: t1
     zero=0
  end function zero
  real(8) pure function beta(t1)  ! length scale of membrane
      real(8),intent(in) :: t1
      beta=beta0-beta1*cos(t1/T)
  end function beta
  real(8) pure function dbeta(t1)	   ! rate of change of length
    real(8),intent(in) :: t1
    dbeta=+beta1*sin(t1/T)/T
  end function dbeta

  subroutine line()
    real,allocatable :: pin(:,:),pout(:,:),a(:,:),au(:,:),au2(:,:),p(:,:)
    real :: u(1,n(1)),du(1,n(1))
    real :: cf(3),pdsi,pdse,dudv,udse
    integer :: i,s,e
    pin = geom%pforce_xslice(flow%pressure,0)/(0.5*A0)
    pout = geom%pforce_xslice(flow%pressure,1)/(0.5*A0)
    a = geom%flux_xslice(flow%pressure,0)/(0.5*A0)
    p = geom%flux_xslice(flow%pressure,1)/(0.5*A0)
    au = geom%flux_xslice(flow%velocity%e(1),1)/(0.5*A0)
    au2 = geom%flux_xslice(flow%velocity%e(1),2)/(0.5*A0)
    if(mod(abs(flow%time),dprnt)<dt) then
      s = xg(1)%hash(nint(-L-30)); e = xg(1)%hash(nint(L+15))
      write(10,'(f10.4,f8.4,10e16.8)') (flow%time/T,xg(1)%x(i)/L, &
        pin(:,i),pout(:,i),a(:,i),au(:,i),au2(:,i),p(:,i),i=s,e)
      flush(10)
    end if

    s = xg(1)%hash(nint(-L-3)); e = xg(1)%hash(nint(L))
    cf = (nu*geom%vforce(flow%velocity)-geom%pforce(flow%pressure))/(0.5*A0)
    u = merge(au/a,0.,a>0)
    du = (u-u0)/flow%dt
    u0 = u
    pdsi = sum(pin(1,s:e))
    pdse = p(1,e)
    dudv = sum(a(1,s:e)*du(1,s:e)*h_min)
    udse = au2(1,e)
    write(11,'(f10.4,4e16.8,2f12.2)') flow%time/T,pdsi,pdse,dudv,udse, &
          100*(-sum(pout(1,s:e))-pdsi-cf(1)), & ! pressure integral check
          100*(dudv+udse-pdsi-pdse)             ! momentum CV check
    flush(11)
  end subroutine line
end program squeeze
