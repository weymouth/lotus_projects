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
  type(fluid)        :: a           ! fluid
  type(body)         :: square      ! body geometry
  integer            :: n(3)
  real               :: fx,fy
!
! -- Init MPI if on
#if MPION
    n = m/b
    call init_mympi(ndims,set_blocks=b(:ndims),set_periodic=p(:ndims))
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
! -- initialize fluid
    call a%init(n,square,V=(/1.,0.,0./),nu=nu,seed=s)
!
! -- run it
    if(mympi_rank()==0) write(6,*) '-- Square Cylinder -- '
    if(mympi_rank()==0) write(6,'("   N=",i0," L=",f0.0," nu=",f0.4)') ndims,L,nu
    a%write = .true.
    do while (a%time<100*L .or. a%write)
       call a%update
       if(mod(a%time,5*L)<a%dt) a%write = .true.
!
! -- get force
       fx = a%pforce(1)
       fy = a%pforce(2)
       if(ndims==3) then
          fx = fx/real(m(3)); fy = fy/real(m(3))
       end if
       write(9,'(f10.4,f7.4,2e16.8)') a%time/L,a%dt,2.*fx/L,2.*fy/L
       flush(8); flush(9)
    end do
    if(mympi_rank()==0) write(6,*) '--- complete --- '
!
! -- Finalize MPI if on
#if MPION
    call mympi_end
#endif
end program square_cyl
!
