program square_cyl
  use mympiMod
  use fluidMod
  use vectorMod
  use bodyMod
  use geom_shape ! to define the square
  implicit none
  integer,parameter  :: f=2**5             ! resolution  
  real,parameter     :: Re = 410           ! Reynolds number
  integer,parameter  :: ndims=3            ! dimensions
  integer,parameter  :: b(3) = (/4,4,1/)   ! blocks
  integer,parameter  :: d(3) = (/24,16,4/) ! domain size
  logical,parameter  :: p(3) = (/.false.,.false.,.true./)  ! periodic
!
  real,parameter     :: L = f       ! length
  logical,parameter  :: s(3) = p    ! random seed
  integer,parameter  :: m(3) = f*d  ! points
  real,parameter     :: nu = L/Re   ! viscosity
  real,parameter     :: yc = m(2)/2 ! location
  integer            :: n(3)
  real               :: force(3)
!
  type(fluid)        :: flow
  type(body)         :: square
!
! -- Initialize MPI (if MPI is ON)
#if MPION
    call init_mympi(ndims,set_blocks=b(:ndims),set_periodic=p(:ndims))
!
! -- Get grid size
    n = m/b
#else
    n = m
#endif
    if(ndims==2) n(3) = 1
!
! -- Initialize the square geometry
    square = plane(4,1,(/-1,0,0/),(/yc-L/2.,0.,0./),0,0) &
        .and.plane(4,1,(/ 1,0,0/),(/yc+L/2.,0.,0./),0,0) &
        .and.plane(4,1,(/0,-1,0/),(/0.,yc-L/2.,0./),0,0) &
        .and.plane(4,1,(/0, 1,0/),(/0.,yc+L/2.,0./),0,0)
!
! -- Initialize fluid
    call flow%init(n,square,V=(/1.,0.,0./),nu=nu,seed=s)
    call flow%write
    if(mympi_rank()==0) print *, '-- Square Cylinder --'
    if(mympi_rank()==0) print '("   N=",i0," L=",f0.0," nu=",f0.4)', ndims,L,nu
!
! -- Run it
    do while (flow%time<100*L)
       call flow%update
       if(mod(flow%time,5*L)<flow%dt) call flow%write
!
! -- Print force on the square
       force = square%pforce(flow%pressure)
       if(ndims==3) force = force/real(m(3))
       write(9,'(f10.4,f7.4,2e16.8)') flow%time/L,flow%dt,2.*force(1:2)/L
       flush(9)
    end do
    if(mympi_rank()==0) write(6,*) '--- complete --- '
!
! -- Finalize MPI
#if MPION
    call mympi_end
#endif
end program square_cyl
!
