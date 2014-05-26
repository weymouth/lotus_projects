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
! -- Initialize the square
    call body(beta,norm)
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
contains
  subroutine body(beta,norm)
    use geom_shape
    implicit none
    type(set)    :: square
    type(prop)   :: p
    type(vfield),intent(inout) :: beta,norm
    real,pointer,dimension(:,:,:) :: bp,n1,n2
    integer :: d,i,j,k
    real(8) :: x(3)
!
! -- define body
    square = plane(4,1,(/-1,0,0/),(/yc-L/2.,0.,0./),0,0) &
        .and.plane(4,1,(/ 1,0,0/),(/yc+L/2.,0.,0./),0,0) &
        .and.plane(4,1,(/0,-1,0/),(/0.,yc-L/2.,0./),0,0) &
        .and.plane(4,1,(/0, 1,0/),(/0.,yc+L/2.,0./),0,0)
!
! -- fill arrays
    call beta%init(n)
    do d=1,ndims
       bp => beta%e(d)%point()
       do concurrent (i=1:n(1),j=1:n(2),k=1:n(3))
          x = beta%e(d)%pos(i,j,k)
          p = square%at(x)
          bp(i,j,k) = mu0(p%distance)
       end do
    end do
    call beta%applyBC

    call norm%init(n,centered=.true.)
    n1 => norm%e(1)%point() 
    n2 => norm%e(2)%point() 
    eps = .55
    do concurrent (i=1:n(1),j=1:n(2),k=1:n(3))
       x = norm%e(1)%pos(i,j,k)
       p = square%at(x)
       n1(i,j,k) = p%normal(1)*mu1(p%distance)
       n2(i,j,k) = p%normal(2)*mu1(p%distance)
    end do
  end subroutine body
!
  pure real function mu0(dis)
    implicit none 
    real(8),intent(in) :: dis
    mu0 = 0
    if(dis>0) mu0 = 1
  end function mu0
!
  pure real function mu1(dis)
    implicit none 
    real(8),intent(in) :: dis
    mu1 = 0
    if(abs(dis-0.5)<0.01) mu1 = 1
  end function mu1
end program square_cyl
!
