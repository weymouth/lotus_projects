module linearAlgebra
  real(8),parameter :: sqrt3 = sqrt(3.D0), pi = acos(-1.D0)
contains
!
! -- identify vortex cores with the \lambda_2 criterion from Jeong and Hussain JFM 1995 
!
  subroutine lambda2(u,lam2)
    use global
    use utility, only: median
    use grid,    only: ddx_face,ddx_cent
    use mympi,   only: mympi_scalar
    implicit none
    real(8),intent(in)  :: u(ndims,ni,nj,nk)
    real(8),intent(out) :: lam2(ni,nj,nk)
    real(8) :: A(3,3),lam(3),gradu(3,3),sym(3,3),anti(3,3),uij(3,3,ni,nj,nk)
    integer :: i,j,k,ii,jj
!
! -- get u_{i,j}
    do i=1,ndims
       call ddx_face(i,u(i,:,:,:),uij(i,i,:,:,:))    ! uii @cell center
       do j=1,ndims
          if(i.eq.j) cycle
          call ddx_cent(j,u(i,:,:,:),uij(i,j,:,:,:)) ! uij @cell corners
       end do
    end do
!
! -- loop through points
    do i=is,ie
    do j=js,je
    do k=ks,ke
!
! -- interpolate uij to the cell centers
       do ii=1,3
       do jj=1,3
          gradu(ii,jj) = uij(ii,jj,i,j,k)
          if(ii.eq.jj) cycle
          gradu(ii,jj) = 0.125*(gradu(ii,jj)+uij(ii,jj,i+1,j+1,k+1) &
               +uij(ii,jj,i+1,j,k)+uij(ii,jj,i,j+1,k)+uij(ii,jj,i,j,k+1) &
               +uij(ii,jj,i+1,j+1,k)+uij(ii,jj,i,j+1,k+1)+uij(ii,jj,i+1,j,k+1))
       end do
       end do
!
! -- get the symmetric and antisymmetric parts (scaling doesn't matter)
       sym = gradu+transpose(gradu)
       anti = gradu-transpose(gradu)
!
! -- contruct A and find the second eigenvalue
       A = matmul(sym,sym)+matmul(anti,anti)
       call sym_cubic_eigen(A,lam)
       call median(lam2(i,j,k),lam(1),lam(2),lam(3))
    end do
    end do
    end do
!
! -- update boundaries (and set ghost to 1, not 0)
    lam2 = lam2-1
    call mympi_scalar(lam2)
    lam2 = lam2+1
  end subroutine lambda2
!
! -- Eigenvalues of a symmetric 3x3 matrix from Smith 1961
!
  subroutine sym_cubic_eigen(A,lam)
    implicit none
    real(8),intent(in)  :: A(3,3)
    real(8),intent(out) :: lam(3)
    real(8) :: m,K(3,3),q,p,phi,cphi,sphi,sqrtp
    m = (A(1,1)+A(2,2)+A(3,3))/3.D0
    K = A-m*reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
    q = det(K)*0.5D0
    p = sum(K*K)/6.D0
    sqrtp = sqrt(p)
    phi = min(max(q/(sqrtp**3+tiny(1.D0)),-1.D0),1.D0)
    phi = acos(phi)/3.D0
    if(phi<0.D0) phi = phi+pi/3.D0
    cphi = cos(phi); sphi = sin(phi)
    lam(1) = m+2*sqrtp*cphi
    lam(2) = m-sqrtp*(cphi+sqrt3*sphi)
    lam(3) = m-sqrtp*(cphi-sqrt3*sphi)
  end subroutine sym_cubic_eigen
!
! -- Determinant of a 3x3 matrix
!
  real(8) function det(A)
    real(8),intent(in)  :: A(3,3)
    det =     sum(A(1,:)*cshift(A(2,:), 1)*cshift(A(3,:), 2))
    det = det-sum(A(1,:)*cshift(A(2,:),-1)*cshift(A(3,:),-2))
  end function det  
end module linearAlgebra
!!$program test
!!$  use linearAlgebra
!!$  real(8) :: B(3,3),v(3)
!!$
!!$  B = reshape((/0,1,-1,1,1,0,-1,0,1/),(/3,3/))
!!$  print 1,B
!!$  print *,'----'
!!$  print 1,det(B)
!!$  call sym_cubic_eigen(B,v)
!!$  print *,'----'
!!$  print 1,v
!!$
!!$1 format(3e12.4)
!!$end program test
