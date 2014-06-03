!-------------------------------------------------------!
!--------------- Foil Impulse Test case ----------------!
!-------------------------------------------------------!
program foil_impulse
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use geom_shape  ! to create geometry
  implicit none
  integer,parameter  :: f=2**6             ! resolution  
  real,parameter     :: Re = 1000          ! Reynolds number
  real(8),parameter  :: alpha = 10         ! AOA
  integer,parameter  :: b(3) = (/2,4,2/)   ! blocks
  integer,parameter  :: d(3) = (/6,4,6/)   ! domain size
!
  integer,parameter  :: ndims = 3   ! dimensions
  real(8),parameter  :: L = f       ! length
  integer,parameter  :: m(3) = f*d  ! points
  real,parameter     :: nu = L/Re   ! viscosity
  real(8),parameter  :: yc = m(2)/2 ! location
  integer            :: n(3)
  real               :: force(3),area
!
  type(fluid)        :: flow
  type(body)         :: foil
  type(set)          :: geom
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
  if(mympi_rank()==0) print *, '-- Foil Impulse --'
  if(mympi_rank()==0) print '("   L=",i0," nu=",f0.4)', f,nu
!
! -- Initialize the foil geometry
!  surface_debug = .true.
  model_fill = .false.
  info%file = 'naca_square.IGS'
  info%x = (/-4.219,-10.271,0./)
  info%s = 0.36626*L
  geom = (model_init(info) &
       .map.((init_affn()**(/alpha,0.D0,0.D0/))+(/yc,yc,0.D0/)))
  if(ndims==3) geom = geom.and.plane(4,1,(/0,0,-1/),yc,0,0)
  foil = geom
  area = L
  if(ndims==3) area = L*(m(3)-yc)
!
! -- Initialize fluid
  call flow%init(n,foil,V=(/1.,0.,0./),nu=nu)
  call flow%write
  if(mympi_rank()==0) print *, '-- init complete --'
!
! -- Run it
  do while (flow%time<L)
     call flow%update
     if(mod(flow%time,0.05*L)<flow%dt) call flow%write
!
! -- Print force on the square
     force = foil%pforce(flow%pressure)
     write(9,'(f10.4,f7.4,3e16.8)') flow%time/L,flow%dt,2.*force/area
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
