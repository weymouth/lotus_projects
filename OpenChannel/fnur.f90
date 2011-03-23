!-------------------------------------------------------!
!------------ NURBS Curve/Surface Modules --------------!
!-------------------------------------------------------!
module nurbs
  use analytic
  private
  integer,parameter :: nl = 6, np = 2**(nl+1)+1
  logical,parameter :: trueDistance = .false.

  type pvec2d
     real(8) :: u,pos(2),dpdu(2)
  end type pvec2d
  type pvec
     real(8) :: uv(2),pos(3),dpdu(3),dpdv(3)
  end type pvec
  type bspline
     integer             :: n,order
     real(8),allocatable :: knot(:)
  end type bspline
  type curve
     type(bspline)       :: bu
     type(pvec2d)        :: pv(np)
     real(8),allocatable :: xcp(:,:),wght(:)
  end type curve
  type trim
     type(bspline)           :: bu,bv
     type(curve),allocatable :: crv(:)
     type(prop)              :: pr(np,np)
     type(pvec)              :: pv(np,np)
     real(8),allocatable     :: xcp(:,:,:),wght(:,:)
     integer                 :: sign = 1
  end type trim
  type nurbsSet
     private
     real(8)                :: nmax(3),nmin(3)
     type(trim),allocatable :: trm(:)
     logical,allocatable    :: blank(:)
  end type nurbsSet

  interface operator(.at.)
     module procedure curve_prop,curveSet_prop,trim_prop,nurbs_prop
  end interface
  interface operator(.of.)
     module procedure curve_pgu,surf_pgu
  end interface
  interface operator(.find.)
     module procedure curve_pvec,trim_pvec
  end interface

  public :: nurbsSet,init_nurbs,nurbs_plot,operator(.at.)
contains
!-------------------------------------------------------!
!
! -- Initialize types
!
  subroutine init_bspline(b)
    implicit none
    type(bspline) :: b
    logical :: debug = .false.
    b%knot = (b%knot-b%knot(b%order))/(b%knot(b%n+1)-b%knot(b%order))
    ! if(debug) print'(2i3,/,99f5.2)',b%n,b%order,b%knot
  end subroutine init_bspline

  subroutine init_curve(crv,fnum)
    implicit none
    type(curve),intent(out)  :: crv
    integer,intent(in)       :: fnum
    integer :: i,order,n
!
! -- read dims and allocate
    read(fnum,*,end=1) order; crv%bu%order = order
    read(fnum,*,end=1) n;     crv%bu%n = n
    allocate(crv%bu%knot(n+order),crv%wght(n),crv%xcp(2,n))
!
! -- read knots, points and weights
    read(fnum,*,end=1) crv%bu%knot
    do i=1,n
       read(fnum,*,end=1) crv%xcp(:,i),crv%wght(i)
    end do
!
! -- normalize knot and save data points
    call init_bspline(crv%bu)
    do i=1,np
       crv%pv(i) = crv.of.(i-1.D0)/(np-1.D0)
    end do
!
! -- finalize
    rewind(fnum)
    return
1   stop 'init_curve: eof'
  end subroutine init_curve

  subroutine init_surf(srf,fnum)
    implicit none
    type(trim),intent(out) :: srf
    integer,intent(in)     :: fnum
    integer :: i,j,ordu,nu,ordv,nv
!
! -- read dims and allocate
    read(fnum,*,end=1) ordu,ordv
    read(fnum,*,end=1) nu,nv
    srf%bu%order = ordu; srf%bv%order = ordv 
    srf%bu%n = nu      ; srf%bv%n = nv
    allocate(srf%bu%knot(nu+ordu),srf%bv%knot(nv+ordv))
    allocate(srf%xcp(3,nu,nv),srf%wght(nu,nv))
!
! -- read knots, points and weights
    read(fnum,*,end=1) srf%bu%knot
    read(fnum,*,end=1) srf%bv%knot
    do i=1,nu
    do j=1,nv
       read(fnum,*,end=1) srf%xcp(:,i,j),srf%wght(i,j)
    end do
    end do
!
! -- normalize knots and save data points
    call init_bspline(srf%bu)
    call init_bspline(srf%bv)
    do i=1,np
    do j=1,np
       srf%pv(i,j) = srf.of.(/ (i-1.D0)/(np-1.D0), (j-1.D0)/(np-1.D0) /)
    end do
    end do
!
! -- finalize
    rewind(fnum)
    return
1   stop 'init_surf: eof'
  end subroutine init_surf

  subroutine init_trim(trm,fnum)
    implicit none
    type(trim),intent(out) :: trm
    integer,intent(in) :: fnum
    integer    :: i,j,n,fnum2
    type(pvec2d) :: pv
!
! -- read and init surf
    read(fnum,*,end=1) fnum2
    call init_surf(trm,fnum2)
!
! -- read number of curves and allocate
    read(fnum,*,end=1) n
    if(n==0) stop 'init_trim: no curves'
    allocate(trm%crv(n))
!
! -- read and init curves
    do i=1,n
       read(fnum,*,end=1) fnum2
       call init_curve(trm%crv(i),fnum2)
    end do
!
! -- fill the property array
    do i=1,np
    do j=1,np
       trm%pr(i,j) = trm%crv.at.(/i-1.D0,j-1.D0/)/(np-1.D0)
    end do
    end do
!
! -- finalize
    rewind(fnum)
    return
1   stop 'init_trim: eof'
  end subroutine init_trim

  subroutine init_nurbs(nrb,figs,fnor)
    use iges
    implicit none
    type(nurbsSet),intent(out) :: nrb
    integer,intent(in)         :: figs,fnor
    integer :: i,n,fnum=3000,fnum2
!
! -- convert IGES to SRF
    call iges_to_srf(figs,fnum,n)
    if(n.ne.144) stop 'init_nurbs: no trimmed surfaces'
!
! -- read pointer file and allocate
    read(fnum,*,end=1) n
    allocate(nrb%trm(n))
!
! -- init trimmed surfaces
    do i=1,n
       read(fnum,*,end=1) fnum2
       call init_trim(nrb%trm(i),fnum2)
    end do
!
! -- normalize
    call nurbs_norm(nrb,fnor)
!
! -- finalize
    rewind(fnum)
    return
1   stop 'init_trim: eof'
  end subroutine init_nurbs

!-----------------------------------------------------------!
! -- normalize nurbs  (scale, shift, rotate)
!
  subroutine nurbs_norm(nrb,fnum)
    implicit none
    type(nurbsSet),intent(inout) :: nrb
    integer,intent(in)           :: fnum
    integer :: n,i,j,d
    real(8) :: val(3),w
    n = size(nrb%trm)
    allocate(nrb%blank(n))
    nrb%blank = .false.
!
! -- scale
    call nurbs_extents
    read(fnum,*,end=100) val
    if(all(val.ne.0.D0)) then       ! scale each dimension
       val = val/(nrb%nmax-nrb%nmin)
    else if(all(val.eq.0.D0)) then  ! don't scale
       val = 1.D0
    else                            ! scale uniformly
       i = sum(maxloc(val))
       val = val(i)/(nrb%nmax(i)-nrb%nmin(i))
    end if
    call nurbs_scale
!
! -- center, rotate, shift
    call nurbs_extents
    val = -(nrb%nmin+nrb%nmax)*0.5
    call nurbs_shift
    read(fnum,*,end=100) val
    val = val/180*acos(-1.D0)
    call nurbs_rotate
    read(fnum,*,end=100) val
    call nurbs_shift
!
! -- read in optional trimming and surface switches
    read(fnum,*,IOSTAT=i) w
    if(i.ne.0) w=1               ! keep trimming on
    read(fnum,*,IOSTAT=i) nrb%blank
!
! -- recompute pv and set extents with 25% margin
    call nurbs_pv
    call nurbs_extents
    w = median(nrb%nmax-nrb%nmin)*0.25 ! note: might break
    nrb%nmax = nrb%nmax+w; nrb%nmin = nrb%nmin-w
!
! -- finalize
    call nurbs_flip
    rewind(fnum)
    return
100 stop 'nurbs_norm: eof'
  contains

    real(8) function median(r) result(m)
      real(8),intent(in) :: r(3)
      m = min(max(r(1),r(2)),max(r(2),r(3)),max(r(1),r(3)))
    end function median
    subroutine nurbs_extents
      nrb%nmin = huge(1.D0); nrb%nmax = -huge(1.D0)
      do i=1,n
         if(nrb%blank(i)) cycle
         nrb%nmax = max(nrb%nmax,maxval(maxval(nrb%trm(i)%xcp,3),2))
         nrb%nmin = min(nrb%nmin,minval(minval(nrb%trm(i)%xcp,3),2))
      end do
    end subroutine nurbs_extents
    subroutine nurbs_scale
      do i=1,n
         if(nrb%blank(i)) cycle
         do d=1,3
            nrb%trm(i)%xcp(d,:,:) = nrb%trm(i)%xcp(d,:,:)*val(d)
         end do
      end do
    end subroutine nurbs_scale
    subroutine nurbs_shift
      do i=1,n
         if(nrb%blank(i)) cycle
         do d=1,3
            nrb%trm(i)%xcp(d,:,:) = nrb%trm(i)%xcp(d,:,:)+val(d)
         end do
      end do
    end subroutine nurbs_shift
    subroutine nurbs_rotate
      real(8) :: full(3,3)
      full = rot_mat(val(1),3)
      full = matmul(rot_mat(val(2),1),full)
      full = matmul(rot_mat(val(3),3),full)
      do d=1,n
         if(nrb%blank(d)) cycle
         do i=1,nrb%trm(d)%bu%n
         do j=1,nrb%trm(d)%bv%n
            nrb%trm(d)%xcp(:,i,j) = matmul(nrb%trm(d)%xcp(:,i,j),full)
         end do
         end do
      end do
    end subroutine nurbs_rotate
    subroutine nurbs_pv
      do d=1,n
         if(nrb%blank(d)) cycle
         do i=1,np
         do j=1,np
            nrb%trm(d)%pv(i,j) = nrb%trm(d).of.(/i-1.D0,j-1.D0/)/(np-1.D0)
         end do
         end do
         if(w.eq.0) then ! turn off trimming
            nrb%trm(d)%pr = prop(0,huge(1.D0),0,0,0)
         end if
      end do         
    end subroutine nurbs_pv
    subroutine nurbs_flip
      type(prop) :: pr
      do d=1,n
         if(nrb%blank(d)) cycle
         pr = nrb%trm(d).at.val
         if(pr%distance.gt.0) nrb%trm(d)%sign = -1 ! note: could break
      end do
    end subroutine nurbs_flip
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
  end subroutine nurbs_norm

!-------------------------------------------------------!
!
! -- Calculate basis function using deBoor's algorithm
!
  subroutine get_basis(b,u,m,basis,dbasis)
    implicit none
    type(bspline),intent(in) :: b
    real(8),intent(in)       :: u
    integer,intent(out)      :: m
    real(8),intent(out)      :: basis(b%order),dbasis(b%order)
    integer :: i,j
    real(8) :: f1,f2,bo(b%order+1),dbo(b%order+1)
    basis = 0 ; bo = 0 ; dbasis = 0 ; dbo = 0
    basis(b%order) = 1
!
! -- locate the active segment
!
    do m=0,b%n-b%order
       if(b%knot(m+b%order+1).ge.u) exit
    end do
!
! -- get interpolant and derivative
!
    recursive_loop: do j=2,b%order
        bo(1:b%order) =  basis
       dbo(1:b%order) = dbasis
       point_loop: do i=1,b%order
          f1 = b%knot(m+i+j-1)-b%knot(m+i)
          f2 = b%knot(m+i+j)-b%knot(m+i+1)
          if(f1.ne.0) f1 = 1./f1
          if(f2.ne.0) f2 = 1./f2
          basis(i) = (u-b%knot(m+i  ))*bo(i  )*f1 &
                    -(u-b%knot(m+i+j))*bo(i+1)*f2
          dbasis(i) = (bo(i  )+(u-b%knot(m+i  ))*dbo(i  ))*f1 &
                     -(bo(i+1)+(u-b%knot(m+i+j))*dbo(i+1))*f2
       end do point_loop
    end do recursive_loop
    m = m+1
  end subroutine get_basis

!-------------------------------------------------------!
!
! -- Calculate pvec given u,v on a curve and surface
!
  type(pvec2d) function curve_pgu(c,u) result(p)
    implicit none    
    type(curve),intent(in) :: c
    real(8),intent(in)     :: u
    integer :: d,m,mm
    real(8),dimension(c%bu%order) :: basis,dbasis,ww,xw
    real(8) :: den,num,den_du,num_du
    logical :: debug = .false.
!
! -- get segment locations and basis functions
!
    call get_basis(c%bu,u,m,basis,dbasis)
    mm = m+c%bu%order-1
    ! if(debug) print 1,u,m,mm,basis,dbasis
1   format(f6.2,2i3,/,4e12.4,/,4e12.4)
!
! -- get denominator coeffs
!
    ww = c%wght(m:mm)
    den =    sum( basis*ww)
    den_du = sum(dbasis*ww)
    ! if(debug) print 3,den,den_du
3   format(2e12.4)
!
! -- get numerator coeffs
!
    do d=1,2
       xw = ww*c%xcp(d,m:mm)
       num = sum( basis*xw)
       num_du = sum(dbasis*xw)
       ! if(debug) print 4,d,xw,num,num_du
4      format(i3,/,4e12.4,/,2e12.4)
!
! -- get pvec
!
       p%u       = u
       p%pos(d)  = num/den
       p%dpdu(d) = (num_du-p%pos(d)*den_du)/den
    end do
  end function curve_pgu

  type(pvec) function surf_pgu(srf,uv) result(p)
    implicit none    
    type(trim),intent(in) :: srf
    real(8),intent(in)    :: uv(2)
    integer :: d,m,mm,n,nn
    real(8),dimension(srf%bu%order) :: bu,du
    real(8),dimension(srf%bv%order) :: bv,dv
    real(8),dimension(srf%bu%order,srf%bv%order) :: ww,xw
    real(8) :: den,num,den_du,num_du,den_dv,num_dv
    logical :: debug = .false.
!
! -- get segment locations and basis functions
!
    call get_basis(srf%bu,uv(1),m,bu,du)
    mm = m+srf%bu%order-1
    ! if(debug) print 1,uv(1),m,mm,bu,dv
!
    call get_basis(srf%bv,uv(2),n,bv,dv)
    nn = n+srf%bv%order-1
    ! if(debug) print 1,uv(2),n,nn,bv,dv
1   format(f6.2,2i3,/,4e12.4,/,4e12.4)
!
! -- get denominator coeffs
!
    ww = srf%wght(m:mm,n:nn)
    den =    sum(bu*matmul(ww,bv))
    den_du = sum(du*matmul(ww,bv))
    den_dv = sum(bu*matmul(ww,dv))
    ! if(debug) print 3,den,den_du,den_dv
3   format(3e12.4)
!
! -- get numerator coeffs
!
    do d=1,3
       xw = ww*srf%xcp(d,m:mm,n:nn)
       num =    sum(bu*matmul(xw,bv))
       num_du = sum(du*matmul(xw,bv))
       num_dv = sum(bu*matmul(xw,dv))
       ! if(debug) print 4,d,xw,num,num_du,num_dv
4      format(i3,/,4e12.4,/,4e12.4,/,4e12.4,/,4e12.4,/,3e12.4)
!
! -- get p (value), dp/du and dp/dv
!
       p%uv      = uv
       p%pos(d)  = num/den
       p%dpdu(d) = (num_du-p%pos(d)*den_du)/den
       p%dpdv(d) = (num_dv-p%pos(d)*den_dv)/den
    end do
  end function surf_pgu

!-------------------------------------------------------!
!
! -- Print out in tecplot format
!
  subroutine curve_plot(crv,fnum,fxcp)
    implicit none
    type(curve),intent(in) :: crv
    integer,intent(in) :: fnum
    integer,intent(in),optional :: fxcp
    integer :: i,n
!
! -- grid points
    write(fnum,*)'VARIABLES=x,y,u'
    write(fnum,1) np
    write(fnum,2) (crv%pv(i)%pos,crv%pv(i)%u,i=1,np)
!
! -- control points
    if(present(fxcp))then
       n = crv%bu%n
       write(fxcp,*)'VARIABLES=x,y,u'
       write(fxcp,1) n
       write(fxcp,2) (crv%xcp(:,i),(i-1.D0)/(n-1.D0),i=1,n)
    end if
1   format('ZONE, I =',i5,', F=POINT')
2   format(3e14.6)
  end subroutine curve_plot

  subroutine trim_plot(trm,fnum,fcrv)
    implicit none
    type(trim),intent(in) :: trm
    integer,intent(in)    :: fnum
    integer,intent(in),optional :: fcrv
    integer i,j
!
! -- srf with distance
    write(fnum,*)'VARIABLES=x,y,z,dis,u,v'
    write(fnum,1) np,np
    write(fnum,2)  ((trm%pv(i,j)%pos,trm%pr(i,j)%distance, &
         trm%pv(i,j)%uv,i=1,np),j=1,np)
!
! -- curves
    if(present(fcrv)) then
       do i=1,size(trm%crv)
          call curve_plot(trm%crv(i),fcrv)
       end do
    end if
1   format('ZONE, I =',i5,', J =',i5,', F=POINT')
2   format(6e14.6)
  end subroutine trim_plot

  subroutine nurbs_plot(nrb,fnum,fcrv)
    implicit none
    type(nurbsSet),intent(in) :: nrb
    integer,intent(in)    :: fnum
    integer,intent(in),optional :: fcrv
    integer :: i
    do i=1,size(nrb%trm)
       if(nrb%blank(i)) cycle
       call trim_plot(nrb%trm(i),fnum,fcrv)
    end do
  end subroutine nurbs_plot

!-----------------------------------------------------------!
! -- Global Nested Grid Search for pvec closest to 'r'
!
! Note: This search is nested which means it is faster than 
! the full global search but not guaranteed to return the 
! best global value. The chance of suboptimum results is 
! minimized by a smooth error function
!
  type(pvec2d) function curve_pvec(crv,r) result(pv)
    implicit none
    type(curve),intent(in) :: crv
    real(8), intent(in)    :: r(2)
    integer :: uc,l,skip,i,j,ind(5)
    real(8) :: phi(5),dis,dx
    uc = 1
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
          if(j< 1) j = j+np+skip-1
          if(j>np) j = j-np-skip+1
          phi(i) = sum((r-crv%pv(j)%pos)**2)
          ind(i) = j
       end do
!
! -- set the center of the new search window
       uc = ind(sum(minloc(phi)))
    end do
    pv = crv%pv(uc)
  end function curve_pvec

  type(pvec) function trim_pvec(trm,r) result(pv)
    implicit none
    type(trim),intent(in) :: trm
    real(8),intent(in)    :: r(3)
    integer :: l,skip,i,j,ind(2,5,5)
    integer,dimension(2) :: ij,ijc,fx
    real(8) :: d2,del,phi(5,5)
    ijc = 1
!
! loop from coarse to fine
    do l = 1,nl
!
! -- get the index skip
       skip = 2**(nl-l)
       fx = np+skip-1
!
! -- fill the squared distance array
       do i=1,5          
       do j=1,5
          ij = ijc+skip*(/i-3,j-3/)
          where(ij< 1) ij = ij+fx
          where(ij>np) ij = ij-fx
          d2 = sum((r-trm%pv(ij(1),ij(2))%pos)**2)
          del = switch(trm%pr(ij(1),ij(2))%distance,dble(skip)/(np-1.D0))
          phi(i,j) = d2*(1.D0+del)
          ind(:,i,j) = ij
       end do
       end do
!
! -- set the center of the new search window
       ij = minloc(phi)
       ijc = ind(:,ij(1),ij(2))
    end do
    pv = trm%pv(ijc(1),ijc(2))
  contains
    real(8) function switch(d,s)
      implicit none
      real(8),intent(in) :: d,s
      switch = -min(max(d/s,-1.D0),0.D0)
    end function switch
  end function trim_pvec

!-----------------------------------------------------------!
! -- Get props nearest r
!
  type(prop) function curve_prop(crv,r) result(p)
    implicit none
    type(curve),intent(in) :: crv
    real(8), intent(in)    :: r(2)
    type(pvec2d) :: pv
    real(8)      :: v(2),n(2),a
    pv = crv.find.r
    v = r-pv%pos
    n = (/-pv%dpdu(2),pv%dpdu(1)/)
    a = sum(n*n); n = n/merge(sqrt(a),1.D0,a>1e-9) ! unit norm
    a = sum(n*v) ! signed projected distance
    p%distance = merge(sign(sqrt(sum(v*v)),a),a,trueDistance)
    p%normal   = (/ n , 0.D0 /)
    p%velocity = 0 ! unused
    p%kappa    = pv%u ! unused
    p%flag     = 1 ! unused
  end function curve_prop

  type(prop) function curveSet_prop(crv,r) result(p)
    implicit none
    type(curve),intent(in) :: crv(:)
    real(8),intent(in)     :: r(2)
    type(prop) :: q
    integer    :: i
    p = crv(1).at.r
    do i=2,size(crv)
       q = crv(i).at.r
       p = p.or.q
    end do    
  end function curveSet_prop

  type(prop) function trim_prop(trm,r) result(p)
    implicit none
    type(trim),intent(in) :: trm
    real(8),intent(in)    :: r(3)
    type(pvec) :: pv
    real(8)    :: v(3),n(3),a
    pv = trm.find.r
    v = r-pv%pos
    n(1) = pv%dpdu(2)*pv%dpdv(3)-pv%dpdu(3)*pv%dpdv(2)
    n(2) = pv%dpdu(3)*pv%dpdv(1)-pv%dpdu(1)*pv%dpdv(3)
    n(3) = pv%dpdu(1)*pv%dpdv(2)-pv%dpdu(2)*pv%dpdv(1)
    a = sum(n*n) 
    n = -trm%sign*n/merge(sqrt(a),1.D0,a>1e-9) ! unit norm
    a = sum(v*n) ! signed projected distance
    p%distance = merge(sign(sqrt(sum(v*v)),a),a,trueDistance)
    p%normal   = n
    p%velocity = 0 ! note: unknown
    p%kappa    = 0 ! note: this requires d2/du2(p)
    p%flag     = 1 ! note: solid
  end function trim_prop

  type(prop) function nurbs_prop(nrb,r) result(p)
    implicit none
    type(nurbsSet),intent(in) :: nrb
    real(8),intent(in)        :: r(3)
    type(prop) :: q
    integer    :: i
    p = prop(0,huge(1.D0),0,0,0) ! default to far far away
    if(.not.any(r>nrb%nmax.or.r<nrb%nmin)) then ! in the box
       do i=1,size(nrb%trm)
          if(nrb%blank(i)) cycle
          q = nrb%trm(i).at.r
          if(p%flag.eq.0) p = q
          p = p.and.q
       end do
    end if
  end function nurbs_prop
end module nurbs
!!$program test
!!$  use nurbs
!!$  call nurbs_test  
!!$end program test
