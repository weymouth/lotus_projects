program square_cyl
  use mympiMod
  use fluidMod
  use vectorMod
  implicit none
  integer,parameter  :: f=2**4             ! resolution  
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
  type(vfield)       :: beta,norm   ! body indicator and normal
  integer            :: n(3),i,j,k
  real               :: fx,fy,x(3)
  real,pointer,dimension(:,:,:) :: bx,by,bz,nx,ny
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
! -- Initialize the square
    call beta%init(n)
    bx => beta%e(1)%point()
    by => beta%e(2)%point()
    if(ndims==3) bz => beta%e(3)%point()
    call norm%init(n)
    nx => norm%e(1)%point()
    ny => norm%e(2)%point()
    do concurrent (i=1:n(1),j=1:n(2),k=1:n(3))
       x = beta%e(1)%pos(i,j,k,1)-yc               ! x component
       if(     x(1) < -0.5*L .or. x(1) > 0.5*L &
          .or. x(2) < -0.5*L .or. x(2) > 0.5*L) bx(i,j,k) = 1
       if(x(2) > -0.5*L .and. x(2) < 0.5*L) then
          if( x(1) ==  0.5*L   ) nx(i,j,k) =  1
          if( x(1) == -0.5*L-1 ) nx(i,j,k) = -1
       end if
       x = beta%e(2)%pos(i,j,k,2)-yc               ! y component
       if(     x(1) < -0.5*L .or. x(1) > 0.5*L &
          .or. x(2) < -0.5*L .or. x(2) > 0.5*L) by(i,j,k) = 1
       if(x(1) > -0.5*L .and. x(1) < 0.5*L) then
          if( x(2) ==  0.5*L   ) ny(i,j,k) =  1
          if( x(2) == -0.5*L-1 ) ny(i,j,k) = -1
       end if
       if(ndims==3) then
          x = beta%e(3)%pos(i,j,k,3)-yc            ! z component
          if(     x(1) < -0.5*L .or. x(1) > 0.5*L &
             .or. x(2) < -0.5*L .or. x(2) > 0.5*L) bz(i,j,k) = 1
       end if
    end do
    call beta%applyBC
!
! -- initialize fluid
    call a%init(n,beta,V=(/1.,0.,0./),nu=nu,seed=s)
    deallocate(beta%e)
!
! -- run it
    if(mympi_rank()==0) write(6,*) '-- Square Cylinder -- '
    if(mympi_rank()==0) write(6,'("   N=",i0," L=",f0.0," nu=",f0.4)') ndims,L,nu
    a%write = .true.
    do while (a%time<5*L .or. a%write)
       call a%update
       if(mod(a%time,L)<a%dt) a%write = .true.
!
! -- get force
       fx = a%pressure%inner(norm%e(1))
       fy = a%pressure%inner(norm%e(2))
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
