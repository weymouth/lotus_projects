!-------------------------------------------------------!
!------------ NURBS Curve/Surface Modules --------------!
!-------------------------------------------------------!
! -- Gabriel Weymouth 
! -- MIT, 02/16/06 - Initial
! -- MIT, 05/01/07 - Object orriented rewrite
! -- MIT, 07/05/07 - Further OOP conversion
! -- MIT, 11/25/07 - Added 2D+t velocity calc
! -- MIT, 11/28/07 - Added NURBS fitting routine
! -- DSS, 10/13/08 - Added distance inversion routine
! -- DSS, 10/21/08 - Added NURBS curves
! -- DSS, 10/22/08 - Added trimmed NURBS
! -- DSS, 11/07/08 - Added array distance routines (yeah Obama)
!
!-------------------------------------------------------!
!
! -- Calculate 1D basis function using deBoor's algorithm
!
module nurbs_1D
contains
  subroutine nurbs_basis(n,order,knot,u,m,basis,dbasis)
    implicit none
    integer,intent(in) :: n,order
    real(8),intent(in) :: knot(n+order),u
    integer,intent(out) :: m
    real(8),intent(out) :: basis(order),dbasis(order)
    integer :: i,j
    real(8) :: f1,f2,bo(order+1),dbo(order+1)
    basis = 0 ; bo = 0 ; dbasis = 0 ; dbo = 0
    basis(order) = 1
!
! -- locate the active segment
!
    do m=0,n-order
       if(knot(m+order+1).ge.u) exit
    end do
!
! -- get interpolant and derivative
!
    recursive_loop: do j=2,order
        bo(1:order) =  basis
       dbo(1:order) = dbasis
       point_loop: do i=1,order
          f1 = knot(m+i+j-1)-knot(m+i)
          f2 = knot(m+i+j)-knot(m+i+1)
          if(f1.ne.0) f1 = 1./f1
          if(f2.ne.0) f2 = 1./f2
          basis(i) = (u-knot(m+i  ))*bo(i  )*f1 &
                    -(u-knot(m+i+j))*bo(i+1)*f2
          dbasis(i) = (bo(i  )+(u-knot(m+i  ))*dbo(i  ))*f1 &
                     -(bo(i+1)+(u-knot(m+i+j))*dbo(i+1))*f2
       end do point_loop
    end do recursive_loop
    m = m+1
  end subroutine nurbs_basis
!-------------------------------------------------------!
!
! -- Calculate NURBS-curve interpolant and derivative at u
!
  function nurbs_pvec_1D(n,order,knot,xcp,wght,u)
    implicit none    
    integer,intent(in) :: n,order
    real(8),intent(in) :: knot(n+order),xcp(2,n),wght(n),u
    real(8) :: nurbs_pvec_1D(2,2)
!
    integer :: d,m
    real(8),dimension(order) :: basis,dbasis,ww,xw
    real(8) :: den,num,p,den_du,num_du
    logical :: debug=.false.
!
! -- get segment locations and basis functions
!
    call nurbs_basis(n,order,knot,u,m,basis,dbasis)
    if(debug) print 1,n,order,knot,u,m,basis,dbasis
1   format(2i3,/,14f5.2,/,f6.2,i3,/,4e12.4,/,4e12.4)
!
! -- get denominator coeffs
!
    ww = wght(m:m+order-1)
    den =    sum( basis*ww)
    den_du = sum(dbasis*ww)
    if(debug) print 3,den,den_du
3   format(2e12.4)
!
! -- get numerator coeffs
!
    do d=1,2
       xw = ww*xcp(d,m:m+order-1)
       num =    sum( basis*xw)
       num_du = sum(dbasis*xw)
       p   = num/den
    if(debug) print 4,d,xw,num,num_du
4   format(i3,/,4e12.4,/,2e12.4)
!
! -- get p (value) and dp/du
!
       nurbs_pvec_1D(d,1) = p
       nurbs_pvec_1D(d,2) = (num_du-p*den_du)/den
    end do

  end function nurbs_pvec_1D
end module nurbs_1D
!
!-------------------------------------------------------!
! -- Calculate NURBS-surface interpolant and gradient at (u,v)
!
module nurbs_2D
contains
  function nurbs_pvec_2D(nu,nv,ordu,ordv, &
       knotu,knotv,xcp,wght,uv)
    use nurbs_1D
    implicit none    
    integer,intent(in) :: nu,nv,ordu,ordv
    real(8),intent(in) :: knotu(nu+ordu),knotv(nv+ordv)
    real(8),intent(in) :: xcp(3,nu,nv),wght(nu,nv),uv(2)
    real(8) :: nurbs_pvec_2D(3,3)
!
    integer :: d,m,n
    real(8),dimension(ordu) :: bu,du
    real(8),dimension(ordv) :: bv,dv
    real(8),dimension(ordu,ordv) :: ww,xw
    real(8) :: den,num,p,num_du,num_dv,den_du,den_dv
    logical :: debug=.false.
!
! -- get segment locations and basis functions
!
! in 'u'
    call nurbs_basis(nu,ordu,knotu,uv(1),m,bu,du)
    if(debug) print 1,nu,ordu,knotu,uv(1),m,bu,du
1   format(2i3,/,14f5.2,/,f6.2,i3,/,4e12.4,/,4e12.4)
!
! in 'v'
    call nurbs_basis(nv,ordv,knotv,uv(2),n,bv,dv)
    if(debug) print 2,nv,ordv,knotv,uv(2),n,bv,dv
2   format(2i3,/,4f5.2,/,f6.2,i3,/,4e12.4)
!
! -- get denominator coeffs
!
    ww = wght(m:m+ordu-1,n:n+ordv-1)
    den =    sum(bu*matmul(ww,bv))
    den_du = sum(du*matmul(ww,bv))
    den_dv = sum(bu*matmul(ww,dv))
    if(debug) print 3,den,den_du,den_dv
3   format(3e12.4)
!
! -- get numerator coeffs
!
    do d=1,3
       xw = ww*xcp(d,m:m+ordu-1,n:n+ordv-1)
       num =    sum(bu*matmul(xw,bv))
       num_du = sum(du*matmul(xw,bv))
       num_dv = sum(bu*matmul(xw,dv))
       p   = num/den
    if(debug) print 4,d,xw,num,num_du,num_dv
4   format(i3,/,4e12.4,/,4e12.4,/,3e12.4)
!
! -- get p (value), dp/du and dp/dv
!
       nurbs_pvec_2D(d,1) = p
       nurbs_pvec_2D(d,2) = (num_du-p*den_du)/den
       nurbs_pvec_2D(d,3) = (num_dv-p*den_dv)/den
    end do

  end function nurbs_pvec_2D
end module nurbs_2D
!
!-------------------------------------------------------!
!----------------- Main Module -------------------------!
!-------------------------------------------------------!
! -- Define nurbs data types for clean passing
!
module nurbs
  type nurbs_curve
     private
     integer :: ordu,nu,sign
     real(8),allocatable :: knotu(:)
     real(8),allocatable :: xcp(:,:),wght(:)
     real(8) :: u_save=-1,pvec_save(2,2)
  end type nurbs_curve
  type nurbs_surf
     private
     integer :: ordu,ordv,nu,nv,sign
     real(8),allocatable :: knotu(:),knotv(:)
     real(8),allocatable :: xcp(:,:,:),wght(:,:)     
     real(8) :: uv_save(2)=-1,pvec_save(3,3)
  end type nurbs_surf
  type nurbs_line
     private
     integer :: nb,nl
     real(8),allocatable :: xb(:,:)
  end type nurbs_line
  type nurbs_grid
     private
     integer :: nb,nl
     real(8),allocatable :: xb(:,:,:)
  end type nurbs_grid
  type nurbs_trim
     private
     integer :: nc
     type(nurbs_surf) :: surf
     type(nurbs_grid) :: grid
     type(nurbs_curve),allocatable :: curve(:)
     type(nurbs_line) ,allocatable :: line(:)
     real(8),allocatable :: dis(:,:),norm(:,:,:)
  end type nurbs_trim
  private rot_mat,gauss,grad_xcp,swap,switch
  real(8),private :: one=1,zero=0,pi=3.14159
contains
!
! ------------------------------------------------------!
!
! -- Initialize the curve by reading from a CRV file
!
  subroutine init_curve(fnum,curve)
    implicit none
    integer,intent(in) :: fnum
    type(nurbs_curve)  :: curve
    integer :: i
!
! -- read dims and allocate
    read(fnum,'(i6)',end=1) curve%ordu
    read(fnum,'(i6)',end=2) curve%nu  
    allocate(curve%knotu(  curve%nu+curve%ordu),&
             curve%wght (  curve%nu),  &
             curve%xcp  (2,curve%nu))
!
! -- read knots, points and weights
    read(fnum,*,end=3) curve%knotu
    curve%knotu = (curve%knotu-curve%knotu(curve%ordu))&
         /(curve%knotu(curve%nu+1)-curve%knotu(curve%ordu))
    do i=1,curve%nu
       read(fnum,*,end=4) curve%xcp(:,i),curve%wght(i)
    end do
!
! -- reset CRV file
    rewind(fnum)
!
! -- default positive orientation
    curve%sign = 1

    return
1   stop 'ordu'
2   stop 'nu'
3   stop 'knot'
4   stop 'xcg,wght'
  end subroutine init_curve
!
! ------------------------------------------------------!
!
! -- Deallocate a curve
!
  subroutine curve_kill(curve)
    implicit none
    type(nurbs_curve),intent(inout) :: curve
    deallocate(curve%knotu,curve%wght,curve%xcp)
  end subroutine curve_kill
!
! ------------------------------------------------------!
!
! -- Copy a curve
!
  subroutine curve_copy(curve1,curve2)
    implicit none
    type(nurbs_curve),intent(in)  :: curve1
    type(nurbs_curve),intent(out) :: curve2
!
! -- copy dims and allocate
    curve2%ordu = curve1%ordu; curve2%nu = curve1%nu
    allocate(curve2%knotu(  curve2%nu+curve2%ordu),&
             curve2%wght (  curve2%nu),  &
             curve2%xcp  (2,curve2%nu))
!
! -- copy knots, points and weights and orientation
    curve2%knotu = curve1%knotu
    curve2%xcp   = curve1%xcp;   curve2%wght  = curve1%wght
    curve2%sign  = curve1%sign
  end subroutine curve_copy
!
! ------------------------------------------------------!
!
! -- Initialize the surface by reading from a SRF file
!
  subroutine init_surf(fnum,surf)
    implicit none
    integer,intent(in) :: fnum
    type(nurbs_surf)   :: surf
    integer :: i,j
    real(8) :: hold(4)
!
! -- read dims and allocate
    read(fnum,*) surf%ordu,surf%ordv
    read(fnum,*) surf%nu  ,surf%nv
    allocate(surf%knotu(  surf%nu+surf%ordu),&
             surf%knotv(  surf%nv+surf%ordv),&
             surf%wght (  surf%nu,surf%nv),  &
             surf%xcp  (3,surf%nu,surf%nv))
!
! -- read knots, points and weights
    read(fnum,*) surf%knotu
    read(fnum,*) surf%knotv
    surf%knotu = (surf%knotu-surf%knotu(surf%ordu))&
         /(surf%knotu(surf%nu+1)-surf%knotu(surf%ordu))
    surf%knotv = (surf%knotv-surf%knotv(surf%ordv))&
         /(surf%knotv(surf%nv+1)-surf%knotv(surf%ordv))
    do i=1,surf%nu
       do j=1,surf%nv
          read(fnum,*) surf%xcp(:,i,j),surf%wght(i,j)
       end do
    end do
!
! -- reset SRF file
    rewind(fnum)
!
! -- default positive orientation
    surf%sign = 1
  end subroutine init_surf
!
! ------------------------------------------------------!
!
! -- Deallocate a surface
!
  subroutine surf_kill(surf)
    implicit none
    type(nurbs_surf),intent(inout) :: surf
    deallocate(surf%knotu,surf%knotv,surf%wght,surf%xcp)
  end subroutine surf_kill
!
! ------------------------------------------------------!
!
! -- Copy a surface
!
  subroutine surf_copy(surf1,surf2)
    implicit none
    type(nurbs_surf),intent(in)  :: surf1
    type(nurbs_surf),intent(out) :: surf2
!
! -- copy dims and allocate
    surf2%ordu = surf1%ordu; surf2%nu = surf1%nu
    surf2%ordv = surf1%ordv; surf2%nv = surf1%nv
    allocate(surf2%knotu(  surf2%nu+surf2%ordu),&
             surf2%knotv(  surf2%nv+surf2%ordv),&
             surf2%wght (  surf2%nu,surf2%nv),  &
             surf2%xcp  (3,surf2%nu,surf2%nv))
!
! -- copy knots, points and weights and orientation
    surf2%knotu = surf1%knotu; surf2%knotv = surf1%knotv 
    surf2%xcp   = surf1%xcp;   surf2%wght  = surf1%wght
    surf2%sign  = surf1%sign
  end subroutine surf_copy
!
! ------------------------------------------------------!
!
! -- Initialize the trimmed surface by reading from a TRM file
!
  subroutine init_trim(fnum,trim)
    implicit none
    integer,intent(in) :: fnum
    type(nurbs_trim)   :: trim
    integer :: n,fsurf,fcurv
!
! -- read and init surf
    read(fnum,*,end=2) fsurf
    call init_surf(fsurf,trim%surf)
!
! -- read number of curves and allocate
    read(fnum,*,end=1) trim%nc
    allocate(trim%curve(trim%nc),trim%line(trim%nc))
!
! -- read and init curves
    do n=1,trim%nc
       read(fnum,*,end=3) fcurv
       call init_curve(fcurv,trim%curve(n))
       call curve_init_line(3,trim%curve(n),trim%line(n))
    end do
!
! -- reset TRM file
    rewind(fnum)

    return
1   stop 'nc'
2   stop 'surf'
3   stop 'curve'
  end subroutine init_trim
!
!-------------------------------------------------------!
! -- Create uniform (u) line
! note: search depth 'isrch' sets the line resolution
!
  subroutine curve_init_line(isrch,curve,line)
    implicit none
    integer,intent(in) :: isrch
    type(nurbs_curve),intent(inout)  :: curve
    type(nurbs_line),intent(out) :: line
    integer :: i
    real(8) :: u
!
! -- set search depth
!
    line%nl = isrch
!
! -- compute the number of points and allocate
!
    line%nb = 2**(line%nl+1)+1
    allocate(line%xb(2,line%nb))
!
! -- fill the points
!
    do i=1,line%nb
       u = (i-1.)/(line%nb-1.)
       line%xb(:,i) = curve_loc(curve,u)
    end do
    
  end subroutine curve_init_line
!
! ------------------------------------------------------!
!
! -- Deallocate a line
!
  subroutine curve_kill_line(line)
    implicit none
    type(nurbs_line),intent(inout) :: line
    deallocate(line%xb)
  end subroutine curve_kill_line
!
!-------------------------------------------------------!
! -- Create uniform (u,v) grid
! note: search depth 'isrch' sets the grid resolution
!
  subroutine surf_init_grid(isrch,surf,grid)
    implicit none
    integer,intent(in) :: isrch
    type(nurbs_surf),intent(inout)  :: surf
    type(nurbs_grid),intent(out) :: grid
    integer :: i,j
    real(8) :: uv(2)
!
! -- set search depth
!
    grid%nl = isrch
!
! -- compute the number of points and allocate
!
    grid%nb = 2**(grid%nl+1)+1
    allocate(grid%xb(3,grid%nb,grid%nb))
!
! -- fill the points
!
    do i=1,grid%nb
       uv(1) = (i-1.)/(grid%nb-1.)
       do j=1,grid%nb
          uv(2) = (j-1.)/(grid%nb-1.)
          grid%xb(:,i,j) = surf_loc(surf,uv)
       end do
    end do
    
  end subroutine surf_init_grid
!
! ------------------------------------------------------!
!
! -- Deallocate a grid
!
  subroutine surf_kill_grid(grid)
    implicit none
    type(nurbs_grid),intent(inout) :: grid
    deallocate(grid%xb)
  end subroutine surf_kill_grid
!
!-------------------------------------------------------!
! -- Create uniform (u,v) grid and distance and normal
! note: search depth 'isrch' sets the grid resolution
!
  subroutine trim_init_grid(isrch,trim)
    implicit none
    integer,intent(in) :: isrch
    type(nurbs_trim),intent(inout)  :: trim
    integer :: nb,n,i,j
    real(8) :: r(2)
!
! -- init grid
    call surf_init_grid(isrch,trim%surf,trim%grid)
!
! -- allocate and fill the distance function
    nb = trim%grid%nb
    allocate(trim%dis(nb,nb))
    trim%dis = 1e6
    do i=1,nb
    do j=1,nb
       r = (/i-1.,j-1./)/(nb-1.)
       call curve_array_distance(r,trim%line,trim%curve,trim%dis(i,j))
    end do
    end do
!
! -- allocate and fill the normal function
    allocate(trim%norm(3,nb,nb))
    do i=1,nb
    do j=1,nb
       r = (/i-1.,j-1./)/(nb-1.)
       trim%norm(:,i,j) = -surf_norm(trim%surf,r)
    end do
    end do
   
  end subroutine trim_init_grid
!
!-------------------------------------------------------!
! -- Basic operations on the curve geometry
!
! -- flip it
!
  subroutine curve_flip(curve)
    implicit none
    type(nurbs_curve),intent(inout) :: curve
    curve%sign = -1*curve%sign
  end subroutine curve_flip
!
! -- scale it
!
  subroutine curve_scale(factor,curve)
    implicit none
    real(8),intent(in) :: factor(2)
    type(nurbs_curve),intent(inout) :: curve
    integer :: d
    do d=1,2
       curve%xcp(d,:) = curve%xcp(d,:)*factor(d)
    end do
  end subroutine curve_scale
!
! -- shift it
!
  subroutine curve_shift(factor,curve)
    implicit none
    real(8),intent(in) :: factor(2)
    type(nurbs_curve),intent(inout) :: curve
    integer :: d
    do d=1,2
       curve%xcp(d,:) = curve%xcp(d,:)+factor(d)
    end do
  end subroutine curve_shift
!
! -- rotate it
!
  subroutine curve_rotate(center,angle,curve)
    implicit none
    real(8),intent(in) :: center(2),angle
    type(nurbs_curve),intent(inout) :: curve
    real(8) :: full(3,3)
    integer :: i
    call curve_shift(-center,curve)
    full = rot_mat(angle,3)
    do i=1,curve%nu
       curve%xcp(:,i) = matmul(curve%xcp(:,i),full(:2,:2))
    end do
    call curve_shift(center,curve)
  end subroutine curve_rotate
!
!-------------------------------------------------------!
! -- Basic operations on the surface geometry
!
! -- flip it
!
  subroutine surf_flip(surf)
    implicit none
    type(nurbs_surf),intent(inout) :: surf
    surf%sign = -1*surf%sign
  end subroutine surf_flip
!
! -- scale it
!
  subroutine surf_scale(factor,surf)
    implicit none
    real(8),intent(in) :: factor(3)
    type(nurbs_surf),intent(inout) :: surf
    integer :: d
    do d=1,3
       surf%xcp(d,:,:) = surf%xcp(d,:,:)*factor(d)
    end do
  end subroutine surf_scale
!
! -- shift it
!
  subroutine surf_shift(factor,surf)
    implicit none
    real(8),intent(in) :: factor(3)
    type(nurbs_surf),intent(inout) :: surf
    integer :: d
    do d=1,3
       surf%xcp(d,:,:) = surf%xcp(d,:,:)+factor(d)
    end do
  end subroutine surf_shift
!
! -- rotate it
!
  subroutine surf_rotate(center,angles,surf)
    implicit none
    real(8),intent(in) :: center(3),angles(3)
    type(nurbs_surf),intent(inout) :: surf
    real(8) :: full(3,3)
    integer :: i,j
    call surf_shift(-center,surf)
    full = rot_mat(angles(1),3)
    full = matmul(rot_mat(angles(2),1),full)
    full = matmul(rot_mat(angles(3),3),full)
    do i=1,surf%nu
    do j=1,surf%nv
       surf%xcp(:,i,j) = matmul(surf%xcp(:,i,j),full)
    end do
    end do
    call surf_shift(center,surf)
  end subroutine surf_rotate
!
  function rot_mat(x,i) ! rotation matrix
    implicit none
    real(8),intent(in) :: x
    integer,intent(in) :: i
    real(8) :: rot_mat(3,3)
    integer :: j,k
    j = mod(i,3)+1; k=mod(j,3)+1
    rot_mat = 0; rot_mat(i,i) = 1.
    rot_mat(j,j) = cos(x); rot_mat(k,k) =  rot_mat(j,j)
    rot_mat(k,j) = sin(x); rot_mat(j,k) = -rot_mat(k,j)
  end function rot_mat
!
!-------------------------------------------------------!
! -- Basic operations on the trimmed surface
!
! -- flip it
!
  subroutine trim_flip(trim)
    implicit none
    type(nurbs_trim),intent(inout) :: trim
    call surf_flip(trim%surf)
    trim%norm = -trim%norm
  end subroutine trim_flip
!
! -- scale it
!
  subroutine trim_scale(factor,trim)
    implicit none
    real(8),intent(in) :: factor(3)
    type(nurbs_trim),intent(inout) :: trim
    call surf_scale(factor,trim%surf)
  end subroutine trim_scale
!
! -- shift it
!
  subroutine trim_shift(factor,trim)
    implicit none
    real(8),intent(in) :: factor(3)
    type(nurbs_trim),intent(inout) :: trim
    call surf_shift(factor,trim%surf)
  end subroutine trim_shift
!
! -- rotate it
!
  subroutine trim_rotate(center,angles,trim)
    implicit none
    real(8),intent(in) :: center(3),angles(3)
    type(nurbs_trim),intent(inout) :: trim
    call surf_rotate(center,angles,trim%surf)
  end subroutine trim_rotate
!
!-------------------------------------------------------!
! -- Inquiry functions on the curve
!
! -- curve extents
!
  function curve_max(curve)
    type(nurbs_curve),intent(in) :: curve
    real(8),dimension(2) :: curve_max
    curve_max = maxval(curve%xcp,2)
  end function curve_max
  function curve_min(curve)
    type(nurbs_curve),intent(in) :: curve
    real(8),dimension(2) :: curve_min
    curve_min = minval(curve%xcp,2)
  end function curve_min
!
! -- pvec(u) [public version of nurbs_pvec_1D with save]
!
  function curve_nurbs_pvec(curve,u)
    use nurbs_1D
    implicit none
    type(nurbs_curve), intent(inout) :: curve
    real(8),intent(in)               :: u
    real(8) :: curve_nurbs_pvec(2,2)

    if(curve%u_save.eq.u) then
       curve_nurbs_pvec = curve%pvec_save ! use saved vector
    else
       curve_nurbs_pvec = nurbs_pvec_1D(curve%nu,curve%ordu,&
            curve%knotu,curve%xcp,curve%wght,u)
       curve%u_save = u                   ! save this point
       curve%pvec_save = curve_nurbs_pvec !  and vector
    end if
  end function curve_nurbs_pvec
!
! -- p(u)
!
  function curve_loc(curve,u)
    implicit none
    type(nurbs_curve), intent(inout) :: curve
    real(8),intent(in)               :: u
    real(8) :: curve_loc(2),vec(2,2)
    vec = curve_nurbs_pvec(curve,u)
    curve_loc = vec(:,1)
  end function curve_loc
!
! -- n(u) [signed curve unit normal]
!
  function curve_norm(curve,u)
    implicit none
    type(nurbs_curve), intent(inout) :: curve
    real(8),intent(in)           :: u
    real(8) :: curve_norm(2),vec(2,2),n(2),a,one = 1.
    vec = curve_nurbs_pvec(curve,u)
    n = (/-vec(2,2),vec(1,2)/)
    a = sqrt(sum(n**2))
    curve_norm = n/merge(one,a,abs(a).lt.1e-14)*curve%sign
  end function curve_norm
!
! -- plot out a tecplot line to file number 'fnum'
!
  subroutine curve_plot(fnum,curve,lxcp)
    implicit none
    integer, intent(in) :: fnum
    type(nurbs_curve), intent(inout) :: curve
    logical, intent(in),optional :: lxcp
    type(nurbs_line) :: line
    integer i
    call curve_init_line(4,curve,line)
    write(fnum,*)'VARIABLES=x,y,u'
    write(fnum,1) 33
    write(fnum,2) (line%xb(:,i),(i-1.)/32.,i=1,33)
    call curve_kill_line(line)
1   format('ZONE, I =',i5,', F=POINT')
2   format(3e14.6)
    if(present(lxcp))then
    if(lxcp) then
       write(fnum,*)'VARIABLES=x,y,u'
       write(fnum,1) curve%nu
       write(fnum,2) (curve%xcp(:,i),(i-1.)/(curve%nu-1.) &
            ,i=1,curve%nu)
    end if
    end if
  end subroutine curve_plot
!
! -- plot out the CRV file for a curve
!
  subroutine curve_crv(fnum,curve)
    implicit none
    integer, intent(in) :: fnum
    type(nurbs_curve), intent(in) :: curve
    integer i
    write(fnum,'(i6)') curve%ordu
    write(fnum,'(i6)') curve%nu
    if(curve%sign.eq.1) then
       write(fnum,'(99f6.3)') curve%knotu*(curve%nu-curve%ordu+1)
       do i=1,curve%nu
          write(fnum,'(3e16.8)') curve%xcp(:,i),curve%wght(i)
       end do
    else
       i=curve%nu+curve%ordu
       write(fnum,'(99f6.3)') (1.-curve%knotu(i:1:-1))*(curve%nu-curve%ordu+1)
       do i=curve%nu,1,-1
          write(fnum,'(3e16.8)') curve%xcp(:,i),curve%wght(i)
       end do
    end if
  end subroutine curve_crv
!
!-------------------------------------------------------!
! -- Inquiry functions on the surface
!
! -- surface extents
!
  function surf_max(surf)
    type(nurbs_surf),intent(in) :: surf
    real(8),dimension(3) :: surf_max
    surf_max = maxval(maxval(surf%xcp,3),2)
  end function surf_max
  function surf_min(surf)
    type(nurbs_surf),intent(in) :: surf
    real(8),dimension(3) :: surf_min
    surf_min = minval(minval(surf%xcp,3),2)
  end function surf_min
!
! -- pvec(u,v) [public version of nurbs_pvec_1D with save]
!
  function surf_nurbs_pvec(surf,uv)
    use nurbs_2D
    implicit none
    type(nurbs_surf), intent(inout) :: surf
    real(8),intent(in)              :: uv(2)
    real(8) :: surf_nurbs_pvec(3,3)

    if(all(surf%uv_save.eq.uv)) then
       surf_nurbs_pvec = surf%pvec_save ! use saved vector
    else
       surf_nurbs_pvec = nurbs_pvec_2D(surf%nu,surf%nv,surf%ordu,&
            surf%ordv,surf%knotu,surf%knotv,surf%xcp,&
            surf%wght,uv)
       surf%uv_save = uv                ! save this point
       surf%pvec_save = surf_nurbs_pvec !  and vector
    end if
  end function surf_nurbs_pvec
!
! -- p(u,v)
!
  function surf_loc(surf,uv)
    implicit none
    type(nurbs_surf), intent(inout) :: surf
    real(8),intent(in)          :: uv(2)
    real(8) :: surf_loc(3),vec(3,3)
    vec = surf_nurbs_pvec(surf,uv)
    surf_loc = vec(:,1)
  end function surf_loc
!
! -- n(u,v) [signed surface unit normal]
!
  function surf_norm(surf,uv,tdpt)
    implicit none
    type(nurbs_surf), intent(inout) :: surf
    real(8),intent(in)           :: uv(2)
    logical,intent(in),optional  :: tdpt
    real(8) :: surf_norm(3),vec(3,3),t1(3),t2(3)
    real(8) :: n(3),a,one = 1.
    vec = surf_nurbs_pvec(surf,uv)
    t1 = vec(:,2)
    t2 = vec(:,3)
    n = cshift(t1,1)*cshift(t2,-1) &
       -cshift(t2,1)*cshift(t1,-1)  ! n = t1 x t2
    if(present(tdpt)) then 
	if(tdpt) n(1) = 0 ! clear the x-component
    end if
    a = sqrt(sum(n**2))
    surf_norm = n/merge(one,a,abs(a).lt.1e-14)*surf%sign
  end function surf_norm
!
! -- 2D+t velocity in the direction (m1,m2) using U=1
!
  function surf_tdpt_velo(surf,uv,m)
    implicit none
    type(nurbs_surf), intent(inout) :: surf
    real(8),intent(in)          :: uv(2),m(2)
    real(8) :: surf_tdpt_velo(2),a(3,3),b(3),y(3)
    real(8),parameter :: tol = 1e-9
    integer :: err
    a = surf_nurbs_pvec(surf,uv)
    a(:,1) = a(:,2); a(:,2) = a(:,3)
    a(1,3) =  0   ; b(1) = 1 ! set up
    a(2,3) = -m(1); b(2) = 0 ! constraint
    a(3,3) = -m(2); b(3) = 0 ! equations
    if(abs(a(1,1)).lt.tol.and.abs(a(1,2)).lt.tol) then
       y(3) = 2.        ! blunt body, take y(3)=2U
    else
       call gauss(3,a,b,y,err)
       if(err.ne.0) then
          write(9,*) "error in tdpt_velo: a,y'"
          write(9,'(3e12.4)') transpose(a),y
       end if
    end if
    surf_tdpt_velo = m*y(3)
  end function surf_tdpt_velo
!
! -- plot out a tecplot grid to file number 'fnum'
!
  subroutine surf_plot(fnum,surf,lxcp)
    implicit none
    integer, intent(in) :: fnum
    type(nurbs_surf), intent(inout) :: surf
    logical, intent(in),optional :: lxcp
    type(nurbs_grid) :: grid
    integer i,j
    call surf_init_grid(4,surf,grid)
    write(fnum,*)'VARIABLES=x,y,z,u,v'
    write(fnum,1) 33,33
    write(fnum,2) ((grid%xb(:,i,j),(i-1.)/32.,(j-1.)/32. &
         ,i=1,33),j=1,33)
    call surf_kill_grid(grid)
1   format('ZONE, I =',i5,', J =',i5,', F=POINT')
2   format(5e14.6)
    if(present(lxcp))then
    if(lxcp) then
       write(fnum,*)'VARIABLES=x,y,z,u,v'
       write(fnum,1) surf%nu,surf%nv
       write(fnum,2) ((surf%xcp(:,i,j),(i-1.)/(surf%nu-1.) &
            ,(j-1.)/(surf%nv-1.),i=1,surf%nu),j=1,surf%nv)
    end if
    end if
  end subroutine surf_plot
!
! -- plot out the srf file for a surface
!
  subroutine surf_srf(fnum,surf)
    implicit none
    integer, intent(in) :: fnum
    type(nurbs_surf), intent(in) :: surf
    integer i,j
    write(fnum,'(i6,",",i6)') surf%ordu,surf%ordv
    write(fnum,'(i6,",",i6)') surf%nu,surf%nv
    write(fnum,'(99f6.3)') surf%knotu*(surf%nu-surf%ordu+1)
    write(fnum,'(99f6.3)') surf%knotv*(surf%nv-surf%ordv+1)
    do i=1,surf%nu
       do j=1,surf%nv
          write(fnum,'(4e16.8)') surf%xcp(:,i,j),surf%wght(i,j)
       end do
    end do    
  end subroutine surf_srf
!
!-------------------------------------------------------!
! -- Inquiry functions on the trimmed surface
!
! -- surface extents
!
  function trim_max(trim)
    type(nurbs_trim),intent(in) :: trim
    real(8),dimension(3) :: trim_max
    trim_max = surf_max(trim%surf)
  end function trim_max
  function trim_min(trim)
    type(nurbs_trim),intent(in) :: trim
    real(8),dimension(3) :: trim_min
    trim_min = surf_min(trim%surf)
  end function trim_min
!
! -- plot out a tecplot grid to file number 'fnum'
!
  subroutine trim_plot(fnum,trim)
    implicit none
    integer, intent(in) :: fnum
    type(nurbs_trim), intent(inout) :: trim
    integer n,i,j
    write(fnum,*)'VARIABLES=x,y,z,dis,u,v'
    n = trim%grid%nb
    write(fnum,1) n,n
    write(fnum,2) ((trim%grid%xb(:,i,j),trim%dis(i,j),&
         (i-1.)/(n-1),(j-1.)/(n-1),i=1,n),j=1,n)
1   format('ZONE, I =',i5,', J =',i5,', F=POINT')
2   format(6e14.6)
  end subroutine trim_plot
!
!-----------------------------------------------------------!
!
! -- Optimization
!
!-----------------------------------------------------------!
! -- Global Nested Grid Search for closest coordinate to 'r'
!
! Note: This search is nested which means it is faster than 
! the full global search but not guaranteed to return the 
! best global value. The chance of suboptimum results is 
! minimized by a smooth error function
!
  function curve_global_search(line,r)
    implicit none
    type(nurbs_line),intent(in) :: line
    real(8), intent(in)        :: r(2)
    integer :: nb,nl,uc,l,skip,i,j,ind(5)
    real(8) :: curve_global_search,phi(5),dis,dx
    nb = line%nb ; nl = line%nl ; uc = 1
!
! loop from coarse to fine
    do l = 1,nl
!
! -- get the index skip
       skip = 2**(nl-l)
!
! -- fill the squared distance array
       do i=1,5
          j = uc+skip*(i-3)
          j = mod(j-1+nb,nb)+1
          phi(i) = sum((r-line%xb(:,j))**2)
          ind(i) = j
       end do
!
! -- set the center of the new search window
       uc = ind(sum(minloc(phi)))
    end do
    curve_global_search = (uc-1.)/(nb-1.)
  end function curve_global_search
!-----------------------------------------------------------!
  function surf_global_search(grid,r,check)
    implicit none
    type(nurbs_grid),intent(in) :: grid
    real(8), intent(in)        :: r(3)
    integer :: nb,nl,uvc(2),uv(2),l,skip,i,j,loc(2)
    real(8) :: surf_global_search(2),phi(5,5),dis,dx
    real(8),optional :: check

    nb = grid%nb ; nl = grid%nl ; uvc = 1
!
! loop from coarse to fine
    do l = 1,nl
!
! -- get the index skip
       skip = 2**(nl-l)
!
! -- set the corner of the search window
       uv = max(min(uvc-2*skip,nb-4*skip),1)
!
! -- fill the squared distance array
       forall(i=1:5 , j=1:5) phi(i,j) = sum((r-grid% &
            xb(:,uv(1)+skip*(i-1),uv(2)+skip*(j-1)))**2)
!
! -- set the center of the new search window
       loc = minloc(phi)
       uvc = uv+skip*(loc-1)
    end do
    surf_global_search = (uvc-1.)/(nb-1.)
!
! -- check distance and grid spacing
    if(present(check)) then
       dis = phi(loc(1),loc(2))
       i = merge(loc(1)+1,loc(1)-1,loc(1).eq.1)
       j = merge(loc(2)+1,loc(2)-1,loc(2).eq.1)
       dx = dis/sum( &
            (grid%xb(:,uv(1)+skip*(loc(1)-1) &
                      ,uv(2)+skip*(loc(2)-1)) &
            -grid%xb(:,uv(1)+skip*(  i   -1) &
                      ,uv(2)+skip*(  j   -1)))**2)
!
! -- signal if the point is too far
       if(dx.gt.1..and.dis.gt.check) check = -1.
    end if
  end function surf_global_search
!-------------------------------------------------------!
  function trim_search(r,trim)
    implicit none
    type(nurbs_trim),intent(in) :: trim
    real(8), intent(in)        :: r(3)
    integer :: nb,nl,l,skip,i,j,ind(2,5,5)
    integer,dimension(2) :: trim_search,ij,ijc
    real(8) :: d2,del,phi(5,5)
    nb = trim%grid%nb ; nl = trim%grid%nl ; ijc = nb/2
!
! loop from coarse to fine
    do l = 1,nl
!
! -- get the index skip
       skip = 2**(nl-l)
!
! -- fill the squared distance array
       do i=1,5          
       do j=1,5
          ij = ijc+skip*(/i-3,j-3/)
          ij = mod(ij-1+nb,nb)+1
          d2 = sum((r-trim%grid%xb(:,ij(1),ij(2)))**2)
          del = switch(trim%dis(ij(1),ij(2)),dble(skip/(nb-1.)))
          phi(i,j) = d2*(1-del)+del
          ind(:,i,j) = ij
       end do
       end do
!
! -- set the center of the new search window
       ij = minloc(phi)
       ijc = ind(:,ij(1),ij(2))
    end do
    trim_search = ijc
  end function trim_search
  real(8) function switch(d,s)
    implicit none
    real(8),intent(in) :: d,s
    switch = -min(max(d/s,dble(-1)),dble(0))
  end function switch
!
!-------------------------------------------------------!
! -- Local Sequential Quadratic Programing (SQP) Search
!
! Note: 'tol' sets the stopping criterion of the search
!       'tdpt' is a switch which enforces an x_1 = X_1 
!       constraint on the optimum. used in 2D+t sims.
!
  function curve_local_search(curve,us,r)
    implicit none
    type(nurbs_curve), intent(inout) :: curve
    real(8), intent(in)         :: us,r(2)
    integer :: i
    real(8),parameter :: tol = 1e-6
    real(8) :: curve_local_search,u,u0
    real(8) :: A(2,2),phi(2),du,length
    real(8) :: zero=0,one=1
    logical :: debug = .false.,violate
!
! -- get function value and gradient
!
    u0 = us
    search_loop: do i=1,20
       if(debug) write(11,'("iteration",i3)') i
       A = curve_nurbs_pvec(curve,u0)
       phi = A(:,1)-r
       length = sqrt(sum(phi**2))
       if(debug) write(11,*)'pvec, phi and length'
       if(debug) write(11,'(2e12.4)') A,phi,length
!
! -- iterate with Gauss-Newton method
! note: adding damping proportional to the residual length
!       enforces gradient descent for distant points
!
       phi(1) = 2.*sum(A(:,2)*phi)
       A(1,1) = 2.*sum(A(:,2)*A(:,2))+length*10.
       du = phi(1)/A(1,1)
       if(debug) write(11,'("du",/,e12.4)') du
!
! -- limit the step length to 0.2
!
       du = max(min(du,dble(0.2)),dble(-0.2))
       if(debug) write(11,'("limited du",/,e12.4)') du
!
! -- if violating u limits try periodic
!
       violate = (u0.eq.1.and.du.lt.0).or.&
                 (u0.eq.0.and.du.gt.0)
       if(violate) then
          u = mod(u0-du+one,one)
          A = curve_nurbs_pvec(curve,u)
          phi = A(:,1)-r
          if(debug) write(11,'("periodic u",/,2e12.4)') &
               u,sqrt(sum(phi**2))
!
! -- if that doesn't work, exit
          if(sqrt(sum(phi**2)).gt.length) then
             u = u0
             if(debug) write(11, &
                  '("constrained:exiting u",2e12.4)') u
             exit search_loop
          end if
       else
!
! -- else step and check convergence
!
          u = min(max(u0-du,zero),one)
          if(debug) write(11,'("u",/,e12.4)') u
          if((u-u0)**2.lt.tol) exit search_loop
       end if
       u0 = u
    end do search_loop
    curve_local_search = u
  end function curve_local_search
!-------------------------------------------------------!
  function surf_local_search(tdpt,surf,uvs,r)
    implicit none
    logical, intent(in)         :: tdpt
    type(nurbs_surf), intent(inout) :: surf
    real(8), intent(in)         :: uvs(2),r(3)
    integer :: i,d,err,k
    real(8),parameter :: tol = 1e-6
    real(8) :: surf_local_search(2),uv(2),uv0(2)
    real(8) :: A(3,3),phi(3),duv(3),length
    real(8) :: zero=0,one=1
    logical :: debug = .false.,violate(2)
!
! -- get function value and gradient
!
    d = merge(3,2,tdpt)
    uv0 = uvs
    search_loop: do i=1,20
       if(debug) write(11,'("iteration",i3)') i
       A = surf_nurbs_pvec(surf,uv0)
       phi = A(:,1)-r
       length = sqrt(sum(phi**2))
       if(debug) write(11,*)'pvec, phi and length'
       if(debug) write(11,'(3e12.4)') A,phi,length
!
! -- construct Lagrangian using Gauss-Newton method
!
       A(1:2,:) = transpose(A(:,2:3))
       A(3,:) = (/.5,0.,0./)
       phi = 2.*matmul(A,phi)
       A = 2.*matmul(A,transpose(A)) ; A(3,3) = 0
!
! -- Add damping proportional to the residual length
! note: enforces gradient descent for distant points
!
       forall(k=1:2) A(k,k) = A(k,k)+length*10.
!
! --  invert for the step
!
       call gauss(d,A(1:d,1:d),phi(1:d),duv(1:d),err)
       if(debug) write(11,'("duv, err",/,3e12.4,i3)') duv,err
!
! -- constrain the variables if violating uv box limits
!
       violate(1) = (uv0(1).eq.1.and.duv(1).lt.0).or.&
                    (uv0(1).eq.0.and.duv(1).gt.0)
       violate(2) = (uv0(2).eq.1.and.duv(2).lt.0).or.&
                    (uv0(2).eq.0.and.duv(2).gt.0)
       if(all(violate)) then
          uv = uv0
          if(debug) write(11, &
               '("constrained corner:exiting uv",2e12.4)') uv
          exit search_loop
       end if
       if(any(violate)) then
          k = merge(1,2,violate(1))
          A(k,:) = 0.; A(k,k) = 1.; phi(k) = 0.
          call gauss(d,A(1:d,1:d),phi(1:d),duv(1:d),err)
          if(debug) write(11,'("constrained duv, err",/,3e12.4,i3)') duv,err
       end if
!
! -- limit the step length to 0.2
!
       length = max(sqrt(sum(duv(1:2)**2)),tol)
       duv(1:2) = duv(1:2)/length*min(length,0.2*one)
       if(debug) write(11,'("limited duv",/,2e12.4,i3)') duv(1:2)
!
! -- error check
!
       if(err.ne.0) then
          write(11,*) 'error in local_search: r,pvec,uv0'
          write(11,'(3e12.4)') r,surf_nurbs_pvec(surf,uv0),uv0
       end if
!
! -- step and check convergence
!
       uv = min(max(uv0-duv(1:2),zero),one)
       if(debug) write(11,'("uv",/,2e12.4)') uv
       if(sum((uv-uv0)**2).lt.tol) exit search_loop
       uv0 = uv
    end do search_loop
    surf_local_search = uv
  end function surf_local_search
!
!-----------------------------------------------------------!
! -- Sequential Least-Squares Surface Fitting
!
  subroutine surf_fitted(surf_in,data,surf_out)
    implicit none
    type(nurbs_surf), intent(inout)  :: surf_in
    real(8)         , intent(in)  :: data(:,:)
    type(nurbs_surf), intent(out) :: surf_out
    type(nurbs_grid) :: grid
    real(8),allocatable,dimension(:,:) :: data_uv,A,ATA,phi,c,dc
    real(8),allocatable,dimension(:)   :: data_err
    real(8),parameter :: tol = 1e-9
    real(8) :: r(3),uv_old(2),uv(2),dc_mag
    integer :: data_size,c_size,d,j,err,iteration
    logical :: debug = .true.
!
! -- Allocate arrays and copy surf info
!
    data_size = size(data,2)
    if(data_size.ne.3) stop 'surf_fitted: data looks wrong'
    data_size = size(data,1)
    c_size = surf_in%nu*surf_in%nv
    allocate(data_uv(data_size,2),data_err(data_size),A(data_size,c_size))
    allocate(ATA(c_size,c_size),phi(c_size,3))
    allocate(c(c_size,3),dc(c_size,3))
    if(debug) write(9,'("sizes:",2i5)') data_size,c_size
    call surf_copy(surf_in,surf_out)
    do d=1,3
       c(:,d) = reshape(surf_out%xcp(d,:,:),(/c_size/))
    end do
!
! -- Initialize uv coord and jacobian for every data point
!
    call surf_init_grid(3,surf_out,grid)
    init_loop: do j=1,data_size
       r = data(j,:)
       uv_old = surf_global_search(grid,r)
       uv = surf_local_search(.false.,surf_out,uv_old,r)
       data_uv(j,:) = uv
       A(j,:) = grad_xcp(surf_out,uv)
       r = matmul(A(j,:),c)-r
       data_err(j) = sqrt(sum(r**2))
    end do init_loop
    call surf_kill_grid(grid)
    if(debug) then
       call surf_plot(7,surf_out,.true.)
       call data_plot(8,data_size,data,data_uv,data_err)
    end if
!
! -- Iterate the EM-stlye optimization
!
    uv_loop: do iteration=1,100
       if(debug) write(9,'("iteration",i4)') iteration
!
! -- Solve the L-S problem for each coordinate
!
       phi = matmul(transpose(A),data)
       ATA = matmul(transpose(A),A)
       phi = phi-matmul(ATA,c)
!!$       forall(j=1:c_size) ATA(j,j) = ATA(j,j)+10.
!!$       coord_loop: do d=1,3
!!$          call gauss(c_size,ATA,phi(:,d),dc(:,d),err)
!!$          if(err.ne.0) then
!!$             write(9,*) 'error inverting',err
!!$             stop 'shit'
!!$          end if
!!$       end do coord_loop
!
! use grad desc instead
       dc = 0.01*phi
!
! -- Check convergence of control points and update
!
       c = c+dc
       do d=1,3
          surf_out%xcp(d,:,:) = reshape(c(:,d),&
               (/surf_out%nu,surf_out%nv/))
       end do
       dc_mag = sqrt(sum(dc**2))/(3.*c_size)
       if(debug) write(9,'(i4,4e16.8)') iteration,dc_mag &
            ,sum(data_err)/real(data_size,8) &
            ,sqrt(sum(data_err**2)/real(data_size,8)) &
            ,maxval(data_err)
       if(dc_mag.lt.tol) exit uv_loop
!
! -- Get new uv coord and jacobian for every data point
!
       data_loop: do j=1,data_size
          r = data(j,:)
          uv_old = data_uv(j,:)
          uv = surf_local_search(.false.,surf_out,uv_old,r)
          data_uv(j,:) = uv
          A(j,:) = grad_xcp(surf_out,uv)
          r = matmul(A(j,:),c)-r
          data_err(j) = sqrt(sum(r**2))
       end do data_loop
    end do uv_loop
    if(debug) then
       call surf_plot(7,surf_out,.true.)
       call data_plot(8,data_size,data,data_uv,data_err)
    end if

  end subroutine surf_fitted
!
  function grad_xcp(surf,uv)
    use nurbs_1D
    implicit none
    type(nurbs_surf),intent(in) :: surf
    real(8),intent(in) :: uv(2)
    real(8) :: grad_xcp(surf%nu*surf%nv)
    integer :: nu,nv,ordu,ordv
    integer :: m,n,i,j
    real(8) :: bu(surf%ordu),du(surf%ordu)
    real(8) :: bv(surf%ordv),dv(surf%ordv)
    real(8) :: local(surf%ordu,surf%ordv)
    real(8) :: global(surf%nu,surf%nv),denom
    nu = surf%nu; nv = surf%nv; ordu = surf%ordu; ordv = surf%ordv
!
! -- get local and global basis function
!
    call nurbs_basis(nu,ordu,surf%knotu,uv(1),m,bu,du)
    call nurbs_basis(nv,ordv,surf%knotv,uv(2),n,bv,dv)
    forall(i=1:ordu,j=1:ordv) local(i,j) = bu(i)*bv(j)
    global = 0
    global(m:m+ordu-1,n:n+ordv-1) = local
!
! -- get wieghted demonimator and xcp gradient
!
    local = surf%wght(m:m+ordu-1,n:n+ordv-1)
    denom = sum(bu*matmul(local,bv))
    grad_xcp = reshape(global,(/nu*nv/))/denom
  end function grad_xcp
!
  subroutine data_plot(fnum,n,xyz,uv,err)
    implicit none
    integer,intent(in) :: fnum,n
    real(8),intent(in) :: xyz(n,3),uv(n,2),err(n)
    integer :: i
    write(fnum,*)'VARIABLES=x,y,z,u,v,err'
    write(fnum,1) n
    write(fnum,2) (xyz(i,:),uv(i,:),err(i),i=1,n)
1   format('ZONE , I =',i5,', F=POINT')
2   format(6e14.6)
  end subroutine data_plot
!
!-------------------------------------------------------!
! -- Utility functions
!
! -- Gaussian elemination with pivoting
!
  subroutine gauss(n,xi,b,y,err,debug_in)
    implicit none
    integer, intent(in)  :: n
    real(8), intent(in)  :: xi(n,n)
    real(8), intent(in)  :: b(n)
    real(8), intent(out) :: y(n)
    integer, intent(out) :: err
    logical,intent(in),optional :: debug_in
    integer,target  :: pivot(2)
    integer,pointer :: pi,pj
    integer :: i,j,k,record(n)
    real(8) :: c(n),x(n,n)
    logical :: debug = .false.

    x = xi; y = b; err = 0
    pi => pivot(1); pj => pivot(2)
    if(present(debug_in)) debug = debug_in

    if(debug) then
       write(9,*) 'initial'
       do j=1,n
          write(9,'(12e12.4)') x(j,:),y(j)
       end do
    end if

    do i=1,n-1
       pivot = maxloc(abs(x(i:n,i:n)))+i-1
       record(i) = pj
       call swap(x(i,:),x(pi,:)) ! exchange rows
       call swap(x(:,i),x(:,pj)) ! exchange collumns
       call swap(y(i),y(pi))     ! exchange source
       
       if(debug) then
          write(9,*) 'step',i
          write(9,*) 'pivot',pivot
          do j=i,n
             write(9,'(12e12.4)') x(j,i:n),y(j)
          end do
       end if

       if(abs(x(i,i)).le.1e-10) then
          err = i
          return
       end if

       c(i+1:n) = x(i+1:n,i)/x(i,i)
       forall(j=i+1:n , k=i+1:n) &
            x(j,k) = x(j,k)-c(j)*x(i,k)
       y(i+1:n) = y(i+1:n)-c(i+1:n)*y(i)
    enddo

    if(abs(x(n,n)).le.1e-10) then
       err = n
       return
    end if
    do i=n,1,-1
       y(i) = (y(i)-sum(x(i,i+1:n)*y(i+1:n)))/x(i,i)       
    enddo

    do i=n-1,1,-1
       pj = record(i)
       call swap(y(pj),y(i))    ! unswap solution

       if(debug) then
          write(9,*) 'back step',i
          write(9,*) 'pivot',pj
          write(9,'(12e12.4)') y
       end if

    end do

  end subroutine gauss
  elemental subroutine swap(a,b)
    implicit none
    real(8),intent(inout) :: a,b
    real(8) :: c
    c = a; a = b; b = c
  end subroutine swap
!
!-----------------------------------------------------------!
!
! -- Batch routine for the minimum distance from a point
!
  subroutine curve_distance(point,line,curve,d,n_out,loc_out)
    implicit none 
    real(8),intent(in) :: point(2)
    type(nurbs_line),intent(in) :: line
    type(nurbs_curve),intent(inout) :: curve
    real(8),intent(out) :: d
    real(8),intent(out),optional :: n_out(2),loc_out(2)
    real(8) :: coord,loc(2),n(2)
    logical :: debug = .false.
!
! --  Find nearest curve coordinate
!
    coord = curve_global_search(line,point)
    if(debug) write(11,'("global",f8.4)') coord
    coord = curve_local_search(curve,coord,point)
    if(debug) write(11,'("local",f8.4)') coord
!
! -- Get the location and unit normal vector
!
    loc = curve_loc(curve,coord)
    if(debug) write(11,'("location",2e12.4)') loc
    n = curve_norm(curve,coord)
    if(debug) write(11,'("normal",2e12.4)') n
!
! -- Get the signed distance
!
    loc = point-loc
    if(debug) write(11,'("displacement",2e12.4)') loc
    d = sign(sqrt(sum(loc**2)),sum(n*loc))
    if(debug) write(11,'("distance",e12.4)') d

    if(present(n_out)) then
       n_out = n
    end if
    if(present(loc_out)) then
       loc_out = loc
    end if

  end subroutine curve_distance
!-----------------------------------------------------------!
  subroutine curve_array_distance(point,line,curve,dist)
    implicit none 
    real(8),intent(in) :: point(2)
    type(nurbs_line) ,intent(inout) :: line(:)
    type(nurbs_curve),intent(inout) :: curve(:)
    real(8),intent(out) :: dist
!!$    real(8),intent(out),optional :: n_out(2)
    integer :: i,j
    real(8) :: d,n(2),v(2),s,norm(2),vect(2),sgn
!
! -- loop through array of curves
!
    j = size(curve)
    dist = huge(one)
    do i=1,j
!
! -- get props of the closest point on the curve
!
       call curve_distance(point,line(i),curve(i),d,n,v)
       s = sign(one,d)
       d = abs(d)
!
! -- if the curve has equal distance but different sign average the props
! -- to get the sign
!
       if(d.eq.dist.and.s.ne.sgn) then
          norm = norm+n
          vect = vect+v
          sgn  = sign(one,sum(norm*vect))
!
! -- if the curve is closer, keep its info
!
       else if(d.lt.dist) then
          dist = d
          vect = v
          norm = n
          sgn  = s
       end if
    end do
!
! -- get the signed distance and unit vector
!
    dist = dist*sgn
!!$    if(present(n_out)) then
!!$       d = sqrt(sum(vect**2))
!!$       d = merge(d,one,d.lt.1e-6)
!!$       n_out = vect/d*sgn
!!$    end if
    
  end subroutine curve_array_distance
!-----------------------------------------------------------!
  subroutine surf_distance(point,grid,surf,d,n_out,loc_out)
    implicit none 
    real(8),intent(in) :: point(3)
    type(nurbs_grid),intent(in) :: grid
    type(nurbs_surf),intent(inout) :: surf
    real(8),intent(out) :: d
    real(8),intent(out),optional :: n_out(3),loc_out(3)
    real(8) :: coord(2),loc(3),n(3)
    logical :: debug = .false.
!
! --  Find nearest surface coordinate
!
    coord = surf_global_search(grid,point)
    if(debug) write(11,'("global",2f8.4)') coord
    coord = surf_local_search(.false.,surf,coord,point)
    if(debug) write(11,'("local",2f8.4)') coord
!
! -- Get the location and unit normal vector
!
    loc = surf_loc(surf,coord)
    if(debug) write(11,'("location",3e12.4)') loc
    n = -surf_norm(surf,coord)
    if(debug) write(11,'("normal",3e12.4)') n
!
! -- Get the signed distance
!
    loc = point-loc
    if(debug) write(11,'("displacement",3e12.4)') loc
    d = sign(sqrt(sum(loc**2)),sum(n*loc))
    if(debug) write(11,'("distance",e12.4)') d

    if(present(n_out)) then
       n_out = n
    end if
    if(present(loc_out)) then
       loc_out = loc
    end if

  end subroutine surf_distance
!!$!-----------------------------------------------------------!
!!$  subroutine surf_array_distance(point,grid,surf,dist,n_out)
!!$    implicit none 
!!$    real(8),intent(in) :: point(3)
!!$    type(nurbs_grid),intent(in) :: grid(:)
!!$    type(nurbs_surf),intent(inout) :: surf(:)
!!$    real(8),intent(out) :: dist
!!$    real(8),intent(out),optional :: n_out(3)
!!$    integer :: i
!!$    real(8) :: d,n(3),v(3),s,norm(3),vect(3),sgn
!!$!
!!$! -- loop through array of surfaces
!!$!
!!$    dist = huge(one)
!!$    do i=1,size(surf)
!!$!
!!$! -- get props of the closest point on the surface
!!$!
!!$       call surf_distance(point,grid(i),surf(i),d,n,v)
!!$       s = sign(one,d)
!!$       d = abs(d)
!!$!
!!$! -- if the surface has equal distance but different sign average the props
!!$! -- to get the sign
!!$!
!!$       if(d.eq.dist.and.s.ne.sgn) then
!!$          norm = norm+n
!!$          vect = vect+v
!!$          sgn  = sign(one,sum(norm*vect))
!!$!
!!$! -- if the surface is closer, keep its info
!!$!
!!$       else if(d.lt.dist) then
!!$          dist = d
!!$          vect = v
!!$          norm = n
!!$          sgn  = s
!!$       end if
!!$    end do
!!$!
!!$! -- get the signed distance and unit vector
!!$!
!!$    dist = dist*sgn
!!$    if(present(n_out)) then
!!$       d = sqrt(sum(vect**2))
!!$       d = merge(d,one,d.lt.1e-6)
!!$       n_out = vect/d*sgn
!!$    end if
!!$
!!$  end subroutine surf_array_distance
!-----------------------------------------------------------!
  subroutine trim_distance(point,trim,d,n_out,loc_out)
    implicit none 
    real(8),intent(in) :: point(3)
    type(nurbs_trim),intent(inout) :: trim
    real(8),intent(out) :: d
    real(8),intent(out),optional :: n_out(3),loc_out(3)
    integer :: ij(2)
    real(8) :: loc(3),n(3)
    logical :: debug = .false.
!
! --  Find nearest trim coordinate
!
    ij = trim_search(point,trim)
    if(debug) write(11,'("global",2i4)') ij
!
! -- Get the location and unit normal vector
!
    loc = trim%grid%xb(:,ij(1),ij(2))
    if(debug) write(11,'("location",3e12.4)') loc
    n = trim%norm(:,ij(1),ij(2))
    if(debug) write(11,'("normal",3e12.4)') n
!
! -- Get the signed distance
!
    loc = point-loc
    if(debug) write(11,'("displacement",3e12.4)') loc
    d = sign(sqrt(sum(loc**2)),sum(n*loc))
    if(debug) write(11,'("distance",e12.4)') d

    if(present(n_out)) then
       n_out = n
    end if
    if(present(loc_out)) then
       loc_out = loc
    end if
  end subroutine trim_distance
!-----------------------------------------------------------!
  subroutine trim_array_distance(point,tri,dist,n_out)
    implicit none 
    real(8),intent(in) :: point(3)
    type(nurbs_trim),intent(inout) :: tri(:)
    real(8),intent(out) :: dist
    real(8),intent(out),optional :: n_out(3)
    integer :: i
    real(8) :: d,n(3),v(3),s,norm(3),vect(3),sgn
!
! -- loop through array of trim surfaces
!
    dist = huge(one)
    do i=1,size(tri)
!
! -- get props of the closest point on the trim surface
!
       call trim_distance(point,tri(i),d,n,v)
       s = sign(one,d)
       d = abs(d)
!
! -- if the trim surface has equal distance but different sign average the
! -- props to get the sign
!
       if(d.eq.dist.and.s.ne.sgn) then
          norm = norm+n
          vect = vect+v
          sgn  = sign(one,sum(norm*vect))
!
! -- if the trim surface is closer, keep its info
!
       else if(d.lt.dist) then
          dist = d
          vect = v
          norm = n
          sgn  = s
       end if
    end do
!
! -- get the signed distance and unit vector
!
    dist = dist*sgn
    if(present(n_out)) then
       d = sqrt(sum(vect**2))
       d = merge(d,one,d.lt.1e-6)
       n_out = vect/d*sgn
    end if

  end subroutine trim_array_distance
end module nurbs
!
!-------------------------------------------------------!
! -- Test routines for NURBS module
!
!!$program test_gauss
!!$  use mympi
!!$  use nurbs
!!$  implicit none
!!$  real(8) :: a(3,3),b(3),x(3)
!!$  integer :: err
!!$
!!$  call init_mympi
!!$  if(myid.eq.0) then
!!$     a(1,:) = (/0.5,0., 0./)
!!$     a(2,:) = (/0. ,0.,-1./)
!!$     a(3,:) = (/0. ,5., 0./)
!!$     b = (/0.25,0.,0./)
!!$     call gauss(3,a,b,x,err)
!!$     write(9,'(3e12.4,i3)') x,err
!!$  end if
!!$
!!$end program test_gauss
!!$program test
!!$  use nurbs
!!$  implicit none
!!$  type(nurbs_surf) :: surf
!!$  type(nurbs_grid) :: grid
!!$  real(8) :: r(3),uv(2),norm(3),distance
!!$!
!!$  type(nurbs_curve) :: curve,curve2
!!$  type(nurbs_line) :: line
!!$  real(8) :: r2(2),u,norm2(2)
!!$!
!!$  type(nurbs_trim) :: trim
!!$
!!$  open(7,file='wedge.srf')
!!$  open(8,file='wedge.dat')
!!$  call init_surf(7,surf)
!!$  call surf_plot(8,surf)
!!$
!!$  call surf_init_grid(3,surf,grid)
!!$  print 1,surf_max(surf)
!!$  print 1,surf_min(surf)
!!$
!!$  r = (/.3328,0.1053,0.0297/)
!!$  call surf_distance(r,grid,surf,distance,norm,uv)
!!$  print 1,r,norm,distance,uv
!!$1 format(3e12.4)
!!$
!!$  open(9,file='circle.crv')
!!$  open(10,file='circle.dat')
!!$  call init_curve(9,curve)
!!$  call curve_plot(10,curve)
!!$
!!$  call curve_init_line(3,curve,line)
!!$  print 2,curve_max(curve)
!!$  print 2,curve_min(curve)
!!$
!!$  r2 = 0.5
!!$  call curve_distance(r2,line,curve,distance,norm2,u)
!!$  print 2,r2,norm2,distance,u
!!$2 format(2e12.4)
!!$
!!$  call curve_copy(curve,curve2)
!!$  r2 = 0.25
!!$  call curve_scale(r2,curve2)
!!$  r2 = 0.3
!!$  call curve_shift(r2,curve2)
!!$  call curve_flip(curve2)
!!$  open(11,file='circle2.cvr')
!!$  call curve_crv(11,curve2)
!!$  rewind(11)
!!$
!!$  r2 = 0.5
!!$  call curve_distance(r2,line,curve2,distance,norm2,u)
!!$  print 2,r2,norm2,distance,u
!!$
!!$  open(12,file='trimmed_wedge.trm')
!!$  open(13,file='trimmed_wedge.dat')
!!$  call init_trim(12,trim)
!!$  call trim_init_grid(6,trim)
!!$  call trim_plot(13,trim)
!!$
!!$  r = (/.3328,0.1053,0.0297/)
!!$  call trim_distance(r,trim,distance,norm)
!!$  print 1,r,norm,distance
!!$
!!$  r = (/.2,.07,.033/)
!!$  call trim_distance(r,trim,distance,norm)
!!$  print 1,r,norm,distance
!!$
!!$end program test
