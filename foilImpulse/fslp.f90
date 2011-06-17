!-----------------------------------------------------------!
!------------------- Slip Model Routines -------------------!
!-----------------------------------------------------------!
!
! This module holds the slip model data and routines. The
!  primary function is to add the influence of the
!  immersed surface tangential conditions to the momentum
!  equation.
! The slip model in this code is the Nuemann condition
!
!  L{u} = \ddn{u_tau}\hat{tau}+\ddn{u_sigma}\hat{sigma} = 0
!
!  where n, tau, and sigma are the local normal-tangential
!  coordinates.
! L is imposed by adding it to the momentum equation, as
!
!             u+del*L{u}=R,  solve for u
!
!  where del is the immersed surface indicator function.
! L is a linear function, but it would be a pain to construct
!  and invert directly. Instead, we make L lower diagonal by 
!  sorting the points by increasing values of del, This
!  basically projects the fluid slip velocity onto the surface.
! The condition is normal-tangential, only active on
!  a small subset of points, and the derivatives are one-
!  sided. All of this requires using special routines and 
!  data structures in this module, making it less clean.
!
! The routine is parralel and is dependant on the 
!  background modules and the m_mrgrnk module. 
!
module slip
  use global
  implicit none
  private
  public slip_reset,slip_add,slip_sort,slip_operator,slip_inv_operator
  real(8),dimension(:,:,:,:),allocatable :: norm
  real(8),dimension(:),allocatable       :: L_del
  integer,dimension(:,:),allocatable     :: L_list
  logical,dimension(:),allocatable       :: L_pos
  integer :: n_list=0,L_level(0:10)=0
contains
!
!----------------------------------------------------------!
!
! At each time step the module data is reset and the points 
!  in the smoothing width are added one by one. Then the data 
!  is stored and sorted (actually, just ranked).
!
! -- Reset geometry influence arrays
!
  subroutine slip_reset
    if(.not.allocated(norm)) allocate(norm(ndims,ni,nj,nk))
    if(allocated(L_del)) deallocate(L_del,L_list,L_pos)
    allocate(L_del(ndims*ni*nj*nk),L_pos(ndims*ni*nj*nk) &
         ,L_list(ndims*ni*nj*nk,4))
    n_list = 0; norm = 0
    L_del = 0; L_list = 0; L_pos = .false.
  end subroutine slip_reset
!
! -- Add info from new point
!
  subroutine slip_add(point,positive,delta,normal,velocity,flag)
    integer,intent(in)    :: point(4),flag
    logical,intent(in)    :: positive
    real(8),intent(in)    :: delta,normal(ndims)
    real(8),intent(inout) :: velocity(ndims)  ! inout
    integer :: d
!
! -- save normal
    d = point(1)
    norm(d,point(2),point(3),point(4)) = normal(d)
!
! -- save info in a 1D array
    if(delta.ne.zero.and.(flag.eq.-2.or.&
         (flag.eq.-1.and.positive)).and.&
         point(2).ge.is.and.point(2).le.ie.and.&
         point(3).ge.js.and.point(3).le.je.and.&
         point(4).ge.ks.and.point(4).le.ke) then
       n_list = n_list+1
       L_list(n_list,:) = point
       L_del(n_list) = delta
       L_pos(n_list) = positive
!
! -- project the velocity onto the normal
       velocity(d) = sum(normal*velocity)*normal(d)
    end if
  end subroutine slip_add
!
! -- Sort points by increasing 'delta'
!
  subroutine slip_sort
    use m_mrgrnk, only: mrgrnk  !<-- content here
    real    :: del1d(n_list)
    integer :: i,rank(n_list),t_list(n_list,4)
    real(8) :: t_del(n_list)
!
! -- resize the arrays
    del1d = L_del(:n_list)  ! single precision for ranking
    t_del = L_del(:n_list)
    t_list = L_list(:n_list,:)
    deallocate(L_list,L_del)
    allocate(L_list(n_list,4),L_del(n_list))
!
! -- rank the 1d array and split into sets
    call mrgrnk(del1d,rank) ! ranking routine
    L_level(0) = 0
    do i=1,9
       L_level(i) = count(del1d.le.real(1-i*0.1))
    end do
    L_level(10) = n_list
!
! -- fill in the ranked list of 4D points
    do i=1,n_list
       L_list(i,:) = t_list(rank(i),:)
       L_del(i) = t_del(rank(i))
    end do
  end subroutine slip_sort
!
!----------------------------------------------------------!
!
! -- Apply the operator L{u} to points on the list
!
  subroutine slip_operator(u,Lu)
    real(8),dimension(ndims,ni,nj,nk),intent(in)  :: u
    real(8),dimension(ndims,ni,nj,nk),intent(out) :: Lu
    integer :: n,point(4)
    Lu = 0
    do n=1,n_list
       point = L_list(n,:)
       Lu(point(1),point(2),point(3),point(4))= L_del(n)* &
            L_operator(u,point,L_pos(n)) ! see below
    end do
  end subroutine slip_operator
!
!----------------------------------------------------------!
!
! -- Invert the operator by backsubstitution
!
! To 'invert' u+del*L{u}=R for u, we note that at point x:
!
!     L{u,x}=alpha*u{x}+beta{u,points near x away from surf}
!
!  where alpha and beta are unknown coefficients. beta only
!  references point further from the surface than x because
!  of the one-sided derivatives. Thus by solving points 
!  from the outside-in, beta depends only on known values 
!  of u. This performs forward substitution of lower diagonal L.
! To get beta, set u{x}=0 => L{u,x}=beta. To get alpha set
!  u{x} = 1 => L{u,x} = alpha+beta. Then solve for u{x}.
! Making this truly parralel would be a major pain. I use
!  an inexact method by updating the boundary values of u
!  at regular intervals. As long as we progress less than
!  dx/2 towards the surface between updates then we shouldn't
!  reference incorrect values. Since the width is usually less
!  than about 5 points, I do 10 updates.
!
  subroutine slip_inv_operator(rhs,u)
    use mympi, only: mympi_vector,mympi_all
    use utility, only: get
    use inout
    real(8),dimension(ndims,ni,nj,nk),intent(in)  :: rhs
    real(8),dimension(ndims,ni,nj,nk),intent(out),target :: u
    real(8) :: delta,r,alpha,beta,gamma
    real(8),pointer :: u_p
    integer :: l,n,point(4)
    u = rhs
    if(mympi_all(n_list.eq.0)) return
    do l=1,10
       call mympi_vector(u,all=.true.)
       do n=L_level(l-1)+1,L_level(l)
!
! -- get point props
          point = L_list(n,:)
          delta = L_del(n)
          r = get(rhs,point)
          u_p => u(point(1),point(2),point(3),point(4))
!
! -- get coeffs L(u)=alpha*u_p+beta in two steps
          u_p = zero
          beta  = L_operator(u,point,L_pos(n),.false.)
          u_p = one
          alpha = L_operator(u,point,L_pos(n),.true.)-beta
!
! -- solve u+del*L(u)=r
          gamma = one+delta*alpha
          u_p = (r-delta*beta)/gamma
          if(abs(gamma).lt.1e-8) then
             u_p = r
             write(9,*) 'L_inv problem'
          end if
       end do
    end do
  end subroutine slip_inv_operator
!
!----------------------------------------------------------!
!
! -- Apply the L_operator to 'u' at 'point'
!
! The grid and u are Cartesian (\hat{e}) but the equations 
!  are normal-tangential (\hat{n}), requiring manipulation.
!  The tangential projection of u in summation notation is 
!
!                 v1 = (e1 e2-n1 n2)u2
! 
!  where the numbers are the indexes, ie e1=e_{index1}.
!  The normal derivative of u is
!
!                    du1dn = n2 du1d2
!
!  where du1dn=\ddn{u1} and du1d2=\ddx{u1,2}. The operator is
!
!         L1 = (e1 e2-n1 n2)(h n3 du2d3-u2)
!
!  where the h is a scaling factor (dx) and the second u2 is
!  to maintain u_n = U_n while inverting.
! The operator is inverted indirectly by calling this routine
!  twice, once with u{point}=0 and once with u{point}=1. 
!  To save time, when second_time=.true. the d2.ne.d1
!  loops are skipped the values from the previous call (which
!  are unchanged since u1 is not referenced) are used.
!
  real(8) function L_operator(u,point,positive,second_time)
    use grid, only: grid_size
    use utility, only: get,ishift
    implicit none
    real(8),intent(in) :: u(ndims,ni,nj,nk)
    integer,intent(in) :: point(4)
    logical,intent(in) :: positive
    logical,intent(in),optional :: second_time
    integer :: i,j,k,d1,d2,d3
    integer,dimension(4) :: point2,point2p1,point2m1,point3
    real(8) :: direction,n1,n2,n3,du2d3,du2dn,v2
    real(8),save,dimension(3) :: L12
!
! -- get local values (index1)
    d1 = point(1); i = point(2)
    j  = point(3); k = point(4)
    n1 = get(norm,point)
    if(n1.eq.one) then  ! tangent projection is zero
       L_operator = 0.
       return
    end if
    direction = merge(1,-1,positive)
!
! -- loop through vector components (index2)
    loop2: do d2=1,ndims
       point2 = (/d2,i,j,k/)
!
! -- when d1=/=d2
       if(d1.ne.d2) then
! ... and second_time then use saved value of L12
          if(present(second_time)) then
             if(second_time) cycle loop2
          end if
! ... otherwise get local point2 further from surf than point1 (downwind)
          if(n1*direction.gt.zero) then  ! point2 is away from surf in d1
             point2 = point2
          else                           ! point2 is near, shift in d1
             point2 = ishift(point2,-1,d1)
          end if
          n2 = get(norm,point2)      ! possible n2 value
          if(n2*direction.le.zero) then ! point2 is away from surf in d2
             point2 = point2
          else                          ! point2 is near, shift in d2
             point2 = ishift(point2,1,d2)
          end if
       end if
       n2 = get(norm,point2)      ! point2 normal
!
! -- loop through derivative components (index3)
       du2dn = 0
       loop3: do d3=1,ndims
!
! -- get downwind points                  ! n3 always points away from surf
          point3 = ishift(point2,d3,0)    !  to give \ddn the correct sign
          n3 = get(norm,point3)*direction ! n1,n2 can point either way since
                                          !  they only project to the tangent
          if(n3.gt.0) then             
             point2m1 = point2
             point2p1 = ishift(point2, 1,d3)
             point3   = point3
          else
             point2m1 = ishift(point2,-1,d3)
             point2p1 = point2
             point3   = ishift(point3,-1,d3)
          end if
!
! -- get spacing
          du2d3 = grid_size(point3)
          if(d2.ne.d3) then
             point3 = ishift(point3, 1,d3)
             du2d3  = (du2d3+grid_size(point3))*0.5
          end if
!
! -- get derivative
          du2d3 = (get(u,point2p1)-get(u,point2m1))/du2d3
          du2dn = du2dn-du2d3*n3
       end do loop3
!
! -- compute L12 = [e1e2-n1n2][h du2dn-u2] (not summing over d2)
       v2 = grid_size(point2)*du2dn-get(u,point2)
       L12(d2) = -n1*n2*v2
       if(d1.eq.d2) L12(d2) = v2+L12(d2)  ! if(d1=d2) e1e2=1, else e1e2=0
    end do loop2
!
! -- sum over d2 to complete the operator
    L_operator = sum(L12(1:ndims))

  end function L_operator
end module slip
!
! -------------------------------------------------------
