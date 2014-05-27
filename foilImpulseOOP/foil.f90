!-------------------------------------------------------!
!--------------- Foil Impulse Test case ----------------!
!-------------------------------------------------------!
program foil_impulse
  use mympiMod
  use fluidMod
  use vectorMod
  use bodyMod
  use geom_shape ! to define the geomtry
  implicit none
  integer,parameter  :: f=2**4             ! resolution  
  real,parameter     :: Re = 410           ! Reynolds number
  integer,parameter  :: ndims=3            ! dimensions
  integer,parameter  :: b(3) = (/4,1,4/)   ! blocks
  integer,parameter  :: d(3) = (/4,4,4/) ! domain size
!
  real(8),parameter  :: L = f       ! length
  integer,parameter  :: m(3) = f*d  ! points
  real,parameter     :: nu = L/Re   ! viscosity
  real(8),parameter  :: yc = m(2)/2 ! location
  real(8),parameter  :: zc = m(3)/2 ! location
  real(8),parameter  :: alpha = 10  ! AOA
  integer            :: n(3)
  real               :: force(3)
!
  type(fluid)        :: flow
  type(body)         :: foil
  type(model_info)   :: info
!
! -- Initialize MPI (if MPI is ON)
#if MPION
  call init_mympi(ndims,set_blocks=b(:ndims))
!
! -- Get grid size
  n = m/b
#else
  n = m
#endif
  if(ndims==2) n(3) = 1
!
! -- Initialize the foil geometry
  surface_debug = .true.
  model_fill = .false.
  info%file = 'naca_square.IGS'
  info%x = (/-4.219,-10.271,0./)
  info%s = 0.36626*L
  foil = (model_init(info) &
       .map.((init_affn()**(/alpha,0.D0,0.D0/))+(/zc,yc,0.D0/)))&
       .and.plane(4,1,(/0,0,-1/),(/zc,yc,zc/),0,0)
!
! -- Initialize fluid
  call flow%init(n,foil,V=(/1.,0.,0./),nu=nu)
  call flow%write
  stop
  if(mympi_rank()==0) print *, '-- Foil Impulse --'
  if(mympi_rank()==0) print '("   L=",i0," nu=",f0.4)', f,nu
!
! -- Run it
  do while (flow%time<100*L)
     call flow%update
     if(mod(flow%time,5*L)<flow%dt) call flow%write
!
! -- Print force on the square
     force = foil%pforce(flow%pressure)
     if(ndims==3) force = force/real(m(3))
     write(9,'(f10.4,f7.4,3e16.8)') flow%time/L,flow%dt,2.*force/L
     flush(9)
  end do
  if(mympi_rank()==0) write(6,*) '--- complete --- '
!
! -- Finalize MPI
#if MPION
  call mympi_end
#endif
end program foil_impulse
!
