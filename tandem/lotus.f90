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
  real,parameter     :: Re = 200           ! Reynolds number
  real,parameter     :: amp = D            ! amplitude
  real,parameter     :: freq = 0.1         ! freqency
!
  integer,parameter  :: ndims = 2                       ! dimensions
  logical,parameter  :: p(2) = .false.                  ! periodic BCs
  real,parameter     :: nu = D/Re                       ! viscosity
  real,parameter     :: Ufric = sqrt(0.013/Re**(1./7.)) ! friction est.
  real,parameter     :: omega = 2*pi*freq/D             ! friction est.
  integer            :: b(2) = (/4,4/)                  ! blocks
  integer            :: n(3)
  real               :: t1,dt,dtPrint=1./freq/20,force(3)
!
  type(fluid)        :: flow
  type(body)         :: back
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
  n(:2) = composite(D*(/7,6/)/b); n(3) = 1
!
! -- Initialize and print grid
  call xg(1)%init(n(1)*b(1),0.6*D,3.0*D,0.5,f=f,r=1.02,c=10.,d=4.)
  call xg(2)%init(n(2)*b(2),1.8*D,1.8*D,1.0,f=f,r=1.03,c=25.)
  if(mympi_rank()==0) then
     call xg(1)%write
     call xg(2)%write
     print '("   total points=",i0)', product(n(1:2)*b)
  end if
!!$  call mympi_end()
!!$  stop
!
! -- Initialize the foil geometry
  back = cylinder(1,1,3,D/2.,0.,0.,0.).map.init_velocity(circle)
!
! -- Initialize fluid
  call flow%init(n,back,V=(/1.,0.,0./),nu=nu)
  if(mympi_rank()==0) print *, '-- init complete --'
!
! -- Time update loop
  do while (flow%time/D<20/freq)
!
! -- update body and fluid
     dt = flow%dt
     t1 = flow%time+dt
     call back%update(t1)
     call flow%update(back)
!
! -- print force
     force = back%pforce(flow%pressure)
     write(9,'(f10.4,f8.4,3e16.8)') t1/D,dt,2.*force/D
     flush(9)
!
! -- full output
     if(mod(t1,dtPrint*D)<dt) then
        call flow%write(back)
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
    real(8),intent(in) :: ts
    velocity = amp*omega*cos(omega*ts)
  end function velocity
  real(8) pure function height(ts)
    real(8),intent(in) :: ts
    height = amp*sin(omega*ts)
  end function height
!
! -- pflow definitions
  pure function circle(u,x) result(v)
    real(8),intent(in) :: u(3),x(3)
    real(8)            :: v(3)
    real(8)            :: r,theta
    r = sqrt(sum(x**2))
    theta = atan2(x(2),x(1))
    ur = cos(theta)*(1-(0.5*D/r)**2)
    ut = -sin(theta)*(r+(0.5*D/r)**2)
    v(1) = cos(theta)*ur-sin(theta)*ut
    v(2) = sin(theta)*ur+cos(theta)*ut
  end function velocity
end program tandem
