!-------------------------------------------------------!
!------------------------ Cylinder ---------------------!
!-------------------------------------------------------!
program tandem
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  use imageMod,   only: display
  use geom_shape  ! to define geom (set,eps,plane, etc)
  implicit none
  real,parameter     :: D = 64             ! length scale
  real,parameter     :: Re = 160           ! Reynolds number
!
  integer,parameter  :: ndims = 2          ! dimensions
  logical,parameter  :: p(2) = .false.     ! periodic BCs
  real,parameter     :: nu = D/Re          ! viscosity
  integer            :: b(3) = (/4,4,1/)   ! blocks
  integer            :: n(3)               ! number of points
  logical            :: root               ! are you root?
!
  type(fluid)        :: flow
  type(body)         :: geom
!
! -- Initialize MPI (if MPI is ON)
#if MPION
  call init_mympi(ndims,set_blocks=b(:ndims),set_periodic=p(:ndims))
#else
  b=1
#endif
!
! -- Print run info
  root = mympi_rank()==0
  if(root) print '("   D=",f0.4,", nu=",f0.4)',D,nu
!
! -- Initialize array size
  n = composite((/6*D,4*D,0./),prnt=root)
!
! -- Initialize and print grid
  call xg(1)%stretch(n(1),-5*D,-0.6*D,2*D,10*D,h_max=5.,prnt=root)
  call xg(2)%stretch(n(2),-5*D,-D,D,5*D,prnt=root)
!
! -- Initialize the geometry and add the potential flow velocity BC
  geom = cylinder(axis=3,radius=D/2.,center=0.).map.init_velocity(circle)
!
! -- Initialize fluid
  call flow%init(n/b,geom,V=(/1.,0.,0./),nu=nu)
!
! -- Time update loop
  if(root) print *, '-- init complete --'
  do while (flow%time/D<40)
!
! -- update body and fluid
     call geom%update(flow%time)
     call flow%update(geom)
!
! -- print force
     write(9,'(f10.4,f8.4,6e16.8)') flow%time/D,flow%dt,& ! time and CFL
                -geom%pforce(flow%pressure)/(0.5*D),&   ! pressure force
              nu*geom%vforce(flow%velocity)/(0.5*D)     ! viscous force
     flush(9)
     if(mod(flow%time,D)<flow%dt) then
       if(mympi_rank()==0) print *,flow%time/D,flow%dt
       call display(flow%velocity%vorticity_Z(),name='01_out',lim=0.25)
       call flow%write()
     end if
  end do
  if(root) print *, '--- complete ---'
!
! -- Finalize MPI
#if MPION
  call mympi_end
#endif

contains
!
! -- pflow definitions
  pure function dipole(x) result(v)
    implicit none
    real,intent(in) :: x(3)
    real            :: v(3)
    real            :: r2,theta,ur,ut
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
    real,intent(in) :: x(3)
    real            :: v(3),m
    v = dipole(x)
    v(1) = 1+v(1)
    m = sqrt(sum(v**2))*(1-eps/D)**2/2.
    if(m>1) v = v/m ! moderate the singularity
  end function circle
end program tandem
