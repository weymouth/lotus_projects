!-------------------------------------------------------!
!------------------- Tandem Cylinders ------------------!
!-------------------------------------------------------!
program tandem
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  use geom_shape  ! to define geom (set,eps,plane, etc)
  implicit none
  real,parameter     :: f = 2              ! scaling factor
  real,parameter     :: D = 100/f          ! length scale
  real,parameter     :: Re = 160           ! Reynolds number
  real,parameter     :: amp = 0.2*D        ! amplitude
  real,parameter     :: freq = 0.25        ! freqency
  integer,parameter  :: periods = 40       ! number of periods
  logical,parameter  :: pflow = .false.    ! use potential flow tangent velocity?
  logical,parameter  :: upstream = .false. ! place upstream body?
!
  integer,parameter  :: ndims = 2                       ! dimensions
  logical,parameter  :: p(2) = .false.                  ! periodic BCs
  real,parameter     :: nu = D/Re                       ! viscosity
  real,parameter     :: Ufric = sqrt(0.013/Re**(1./7.)) ! friction est.
  real,parameter     :: omega = 2*pi*freq/D             ! friction est.
  integer            :: b(2) = (/4,4/)                  ! blocks
  integer            :: n(3)
  real               :: t1,dt,dtPrint=1./freq,pforce(3),vforce(3)
!
  type(fluid)        :: flow
  type(set)          :: back
  type(body)         :: bodies
!
! -- Initialize MPI (if MPI is ON)
#if MPION
  call init_mympi(ndims,set_blocks=b(:ndims),set_periodic=p(:ndims))
#else
  b=1
#endif
!
! -- Print run info
  if(mympi_rank()==0) print *, '-- Tandem Test --'
  if(mympi_rank()==0) print '("   D=",f0.4,", nu=",f0.4,", y+=",f0.4)',D,nu,Ufric/nu
!
! -- Initialize array size
  n(:2) = composite(D*(/20,10/)/b); n(3) = 1
!
! -- Initialize and print grid
  call xg(1)%init(n(1)*b(1),2.7*D,12*D,1.0,f=f,r=1.02,d=4.)
  call xg(2)%init(n(2)*b(2),3.2*D,3.2*D,1.0,f=f,r=1.02)
  if(mympi_rank()==0) then
     call xg(1)%write(D)
     call xg(2)%write(D)
     print '("   total points=",i0)', product(n(1:2)*b)
  end if
!!$  call mympi_end()
!!$  stop
!
! -- Initialize the tandem geometry
  if(pflow) then
     back = (cylinder(1,1,3,D/2.,0.,0.,0.).map.init_velocity(circle).map.init_rigid(2,height,zip))
  else
     back = (cylinder(1,1,3,D/2.,0.,0.,0.).map.init_rigid(2,height,velocity))
  end if
  if(upstream) then
     bodies = back.or.cylinder(1,0,3,D/2.,(/-5*D,0.,0./),0.,0.) ! no force data
  else
     bodies = back
  end if
!
! -- Initialize fluid
  call flow%init(n,bodies,V=(/1.,0.,0./),nu=nu)
  if(flow%time==0) call flow%write(bodies)
  if(mympi_rank()==0) print *, '-- init complete --'
!
! -- Time update loop
  do while (flow%time/D<periods/freq)
!
! -- update body and fluid
     dt = flow%dt
     t1 = flow%time+dt
     call bodies%update(t1)
     call flow%update(bodies)
!
! -- print force
     pforce = -bodies%pforce(flow%pressure)
     vforce = nu*bodies%vforce(flow%velocity)
     write(9,'(f10.4,f8.4,6e16.8)') t1/D,dt,2.*pforce/D,2.*vforce/D
     flush(9)
!
! -- full output
     if(mod(t1,dtPrint*D)<dt) then
        call flow%write(bodies)
        if(mympi_rank()==0) print 1,t1/D,height(REAL(t1,8))/D,velocity(REAL(t1,8))
1       format("   t=",f0.4," height=",f7.4," velocity=",f7.4)
     end if
  end do
  if(mympi_rank()==0) write(6,*) '--- complete --- '
  !
! -- Finalize MPI
#if MPION
  call mympi_end
#endif

contains
!
! -- motion definitions
  real(8) pure function velocity(ts)
    implicit none
    real(8),intent(in) :: ts
    velocity = amp*omega*cos(omega*ts)
  end function velocity
  real(8) pure function height(ts)
    implicit none
    real(8),intent(in) :: ts
    height = amp*sin(omega*ts)
  end function height
  real(8) pure function zip(ts)
    implicit none
    real(8),intent(in) :: ts
    zip = 0
  end function zip
!
! -- pflow definitions
  pure function dipole(x) result(v)
    implicit none
    real(8),intent(in) :: x(3)
    real(8)            :: v(3)
    real(8)            :: r2,theta,ur,ut
    r2 = sum(x**2)/(0.5*D)**2
    theta = atan2(x(2),x(1))
    ur =  -cos(theta)/r2
    ut =  -sin(theta)/r2
    v(1) = cos(theta)*ur-sin(theta)*ut
    v(2) = sin(theta)*ur+cos(theta)*ut
    v(3) = 0.
  end function dipole
  pure function circle(x) result(v)
    implicit none
    real(8),intent(in) :: x(3)
    real(8)            :: v(3),dip(3),vel,m
    dip = dipole(x)
    vel = velocity(REAL(t1,8))
    v(1) = 1+dip(1)+vel*dip(2)
    v(2) = 0+dip(2)-vel*dip(1)
    v(3) = 0
    m = sqrt(sum(v**2)/(1+vel**2))*(1-eps/D)**2/2.
    if(m>1) v = v/m ! moderate the singularity
  end function circle
end program tandem
