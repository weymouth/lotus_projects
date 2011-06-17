!------------------------------------------------------------!
!------------ Clustered Orthogonal Grid Routines ------------!
!------------------------------------------------------------!
module grid
  real(8),allocatable,dimension(:),private :: &
       xpos,ypos,zpos,dx,dy,dz
contains
!
! -------------------------------------------------------
!
! -- Read grid input and set up grid metrics
!
  subroutine init_grid(gis,gie,gjs,gje,gks,gke)
    use global
    use mympi
    use utility, only: nrerror,log,log_print
    implicit none
    integer,intent(out) :: gis,gie,gjs,gje,gks,gke
    real(8),allocatable,dimension(:) :: ds,spos,dsfi,dsci
    integer :: d,nseg,s,n,ss,se,i
    real(8) :: a,r,in(3),dl
    gis = ie; gie = is; gjs = je; gje = js; gks = ke; gke = ks
!
! -- Allocate arrays
!
    allocate(xpos(ni),ypos(nj),zpos(nk),dx(ni),dy(nj),dz(nk))
    call mympi_read(input_num)    
    do d=1,ndims
       allocate(  ds(1-thk:mympi_points(d)+thk))
       allocate(spos(1-thk:mympi_points(d)+thk))
!
! -- Read starting position
!
       call mympi_read(input_num)    
       call mympi_read(input_num,r=spos(1))    
!
! -- Read #segments and loops
!
       call mympi_read(input_num,i=nseg)    
       se = 0
       do s=1,nseg
!
! -- Reads #points, inital spacing, expansion ratio
!
          call mympi_read(input_num,ra=in)    
          n = in(1); a = in(2); r = in(3)
          ss = se+1
          se = se+n
!
! -- Compute spacing
!          
          forall(i=0:n-1) ds(ss+i) = a*r**i
          write(log,'("  grid params",3e12.4)') &
               ds(ss),ds(se),sum(ds(ss:se))
          call log_print
!
! -- Check extents and fill ghosts
!
          if(s.eq.nseg) then
             ds(1-thk:0) = ds(1)
             if(se.ne.mympi_points(d)) call nrerror('mympi_grid: #points')
             ds(se+1:se+thk) = ds(se)
          end if
       end do
       dl = minval(ds)
!
! -- Compute position
!              
       do i=0,1-thk,-1
          spos(i) = spos(i+1)-ds(i)
       end do
       do i=2,mympi_points(d)+thk
          spos(i) = spos(i-1)+ds(i-1)
       end do
!
! -- Save a window on each processesor in each direction 
!
       if(d.eq.1) then
          ss = mympi_coords(1)*(ni-2*thk)+1-thk
          xpos = spos(ss:ss+ni-1)
          dx   =   ds(ss:ss+ni-1)
          do i=ie-1,is,-1
             if(dx(i).eq.dl) gis = i
          end do
          do i=is,ie
             if(dx(i).eq.dl) gie = i
          end do
       else if(d.eq.2) then
          ss = mympi_coords(2)*(nj-2*thk)+1-thk
          ypos = spos(ss:ss+nj-1)
          dy   =   ds(ss:ss+nj-1)
          do i=je-1,js,-1
             if(dy(i).eq.dl) gjs = i
          end do
          do i=js+1,je
             if(dy(i).eq.dl) gje = i
          end do
       else if(d.eq.3) then
          ss = mympi_coords(3)*(nk-2*thk)+1-thk
          zpos = spos(ss:ss+nk-1)
          dz   =   ds(ss:ss+nk-1)
          do i=ke-1,ks,-1
             if(dz(i).eq.dl) gks = i
          end do
          do i=ks+1,ke
             if(dz(i).eq.dl) gke = i
          end do
       end if
       deallocate(ds,spos)
    end do
!
! -- Set values for 2D simulations
!
    if(ndims.eq.2) then
       zpos = 0.
       dz = 1.
       gks = ks; gke = ke
    end if

  end subroutine init_grid
!
! -------------------------------------------------------
!
! -- Get grid cell volume at (i,j,k)
!
  real(8) pure function grid_vol(i,j,k)
    implicit none
    integer,intent(in) :: i,j,k
    grid_vol = dx(i)*dy(j)*dz(k)
  end function grid_vol
!
! -- Get grid size length (d,i,j,k)
!
  real(8) pure function grid_size(point)
    implicit none
    integer,intent(in) :: point(4)
    if(point(1).eq.1) grid_size = dx(point(2))
    if(point(1).eq.2) grid_size = dy(point(3))
    if(point(1).eq.3) grid_size = dz(point(4))
  end function grid_size
!
! -------------------------------------------------------
!
! -- Get position of variable (d,i,j,k) (d=0 for center)
!
  pure function grid_pos(d,i,j,k)
    use global
    implicit none
    integer,intent(in):: d,i,j,k
    real(8) :: grid_pos(3)
    grid_pos(1) = xpos(i)+merge(zero,dx(i)*0.5,d.eq.1)
    grid_pos(2) = ypos(j)+merge(zero,dy(j)*0.5,d.eq.2)
    grid_pos(3) = zpos(k)+merge(zero,dz(k)*0.5,d.eq.3.or.ndims.eq.2)
  end function grid_pos
!
! -- Get position of cell corner (i,j,k)
!
  pure function grid_corner(i,j,k)
    use global
    implicit none
    integer,intent(in):: i,j,k
    real(8) :: grid_corner(3)
    grid_corner = (/xpos(i),ypos(j),zpos(k)/)
  end function grid_corner
!
! -------------------------------------------------------
!
! -- Get 1D grid coeffs spread to full array size
!
  pure function grid_spacing_array(d) result(ds)
    use global
    implicit none
    integer,intent(in)          :: d
    real(8),dimension(ni,nj,nk) :: ds
    if(d.eq.1) then
       ds = spread(spread(dx,2,nj),3,nk)
    else if(d.eq.2) then
       ds = spread(spread(dy,1,ni),3,nk)
    else
       ds = spread(spread(dz,1,ni),2,nj)
    end if
  end function grid_spacing_array
!
  pure function grid_area_array(d) result(da)
    use global
    implicit none
    integer,intent(in)          :: d
    real(8),dimension(ni,nj,nk) :: da
    if(d.eq.1) then       
       da = spread(spread(dy,2,nk)*spread(dz,1,nj),1,ni)
    else if(d.eq.2) then
       da = spread(spread(dx,2,nk)*spread(dz,1,ni),2,nj)
    else
       da = spread(spread(dx,2,nj)*spread(dy,1,ni),3,nk)
    end if
  end function grid_area_array
!
  pure function grid_position_array(d,d2) result(s)
    use global, only: ni,nj,nk,ndims,zero
    implicit none
    integer,intent(in)          :: d,d2
    real(8),dimension(ni,nj,nk) :: s
    real(8) :: s1d(max(ni,nj,nk))
    if(d.eq.1) then
       s1d(:ni) = xpos+merge(zero,dx*0.5,abs(d2).eq.1)
       s = spread(spread(s1d(:ni),2,nj),3,nk)
    else if(d.eq.2) then
       s1d(:nj) = ypos+merge(zero,dy*0.5,abs(d2).eq.2)
       s = spread(spread(s1d(:nj),1,ni),3,nk)
    else
       s1d(:nk) = zpos+merge(zero,dz*0.5,abs(d2).eq.3.or.ndims.eq.2)
       s = spread(spread(s1d(:nk),1,ni),2,nj)
    end if

  end function grid_position_array
!
! -------------------------------------------------------
!
! -- Scale the normal for each cell
!
  pure subroutine grid_vof_scale(i,j,k,n)
    implicit none
    integer,intent(in) :: i,j,k
    real(8),intent(inout),dimension(3) :: n
    n(1) = n(1)*dx(i)
    n(2) = n(2)*dy(j)
    n(3) = n(3)*dz(k)
  end subroutine grid_vof_scale
!
! -------------------------------------------------------
!
! -- Determine the convection limit number
!
  pure subroutine grid_CFL(u,C)
    use global
    implicit none
    real(8),intent(in),dimension(ndims,ni,nj,nk) :: u
    real(8),intent(out),dimension(ni,nj,nk)      :: C
    integer :: i,j,k

    C = 0
    forall(i=is:ie,j=js:je,k=ks:ke)
       C(i,j,k) = max(u(1,i,j,k),-u(1,i+1,j,k),&
            abs(u(1,i,j,k)-u(1,i+1,j,k)))/dx(i)&
                 +max(u(2,i,j,k),-u(2,i,j+1,k),&
            abs(u(2,i,j,k)-u(2,i,j+1,k)))/dy(j)
    end forall
    if(ndims.eq.3) then
    forall(i=is:ie,j=js:je,k=ks:ke)
       C(i,j,k) = C(i,j,k)&
                 +max(u(3,i,j,k),-u(3,i,j,k+1),&
            abs(u(3,i,j,k)-u(3,i,j,k+1)))/dz(k)
    end forall
    end if

  end subroutine grid_CFL
!
! -------------------------------------------------------
!
! -- Linear interpolation
! d = variable position, d2 = interpolation direction
! Note: The inline case simplifies to an average
!
  subroutine grid_linear(d,d2,a,b)
    use global
    use utility, only: shift
    implicit none
    integer,intent(in) :: d,d2
    real(8),intent(in),dimension(ni,nj,nk) :: a
    real(8),intent(out),dimension(ni,nj,nk) :: b
    real(8),dimension(max(ni,nj,nk)) :: r
    integer :: i,j,k
    b = 0
    if(d.eq.d2) then
       call shift(-1,d2,a,b)
       b = (a+b)*.5
    else
       forall(i=2:ni , d2.eq.1)
          r(i) = dx(i)/(dx(i)+dx(i-1))
          b(i,:,:) = r(i)*a(i,:,:)+(1-r(i))*a(i-1,:,:)
       end forall
       forall(j=2:nj , d2.eq.2)
          r(j) = dy(j)/(dy(j)+dy(j-1))
          b(:,j,:) = r(j)*a(:,j,:)+(1-r(j))*a(:,j-1,:)
       end forall
       forall(k=2:nk , d2.eq.3)
          r(k) = dz(k)/(dz(k)+dz(k-1))
          b(:,:,k) = r(k)*a(:,:,k)+(1-r(k))*a(:,:,k-1)
       end forall
    end if

  end subroutine grid_linear
!
! -------------------------------------------------------
!
! -- Get the derivative for values held on faces
!    
  pure subroutine ddx_face(d,a,b)
    use global
    implicit none
    integer,intent(in) :: d
    real(8),dimension(ni,nj,nk),intent(in) :: a
    real(8),dimension(ni,nj,nk),intent(out) :: b
    integer :: i,j,k
    b = 0
    forall(i=1:ni-1 , d.eq.1) &
         b(i,:,:) = (a(i+1,:,:)-a(i,:,:))/dx(i)
    forall(j=1:nj-1 , d.eq.2) &
         b(:,j,:) = (a(:,j+1,:)-a(:,j,:))/dy(j)
    forall(k=1:nk-1 , d.eq.3) &
         b(:,:,k) = (a(:,:,k+1)-a(:,:,k))/dz(k)
  end subroutine ddx_face
!
! -------------------------------------------------------
!
! -- Get the derivative for values held at the centers
!    
  pure subroutine ddx_cent(d,a,b)
    use global
    implicit none
    integer,intent(in) :: d
    real(8),dimension(ni,nj,nk),intent(in) :: a
    real(8),dimension(ni,nj,nk),intent(out) :: b
    integer :: i,j,k
    b = 0
    forall(i=2:ni , d.eq.1) &
         b(i,:,:) = (a(i,:,:)-a(i-1,:,:))*2./(dx(i)+dx(i-1))
    forall(j=2:nj , d.eq.2) &
         b(:,j,:) = (a(:,j,:)-a(:,j-1,:))*2./(dy(j)+dy(j-1))
    forall(k=2:nk , d.eq.3) &
         b(:,:,k) = (a(:,:,k)-a(:,:,k-1))*2./(dz(k)+dz(k-1))
  end subroutine ddx_cent
!
! -------------------------------------------------------
!
! -- Compute the gradient vector of a scalar
!
  pure subroutine grid_gradient(a,b)
    use global
    implicit none
    integer :: d
    real(8),intent(in),dimension(ni,nj,nk) :: a
    real(8),intent(out),dimension(ndims,ni,nj,nk) :: b
    do d=1,ndims
       call ddx_cent(d,a,b(d,:,:,:))
    end do
  end subroutine grid_gradient
!
! -------------------------------------------------------
!
! -- Take the integrated divergence of a vector field
!
  pure subroutine grid_divergence(a,b)
    use global
    implicit none
    real(8),intent(in),dimension(ndims,ni,nj,nk) :: a
    real(8),intent(out),dimension(ni,nj,nk) :: b
    integer :: i,j,k
    b = 0
    forall(i=is:ie , j=js:je , k=ks:ke) &
         b(i,j,k) = (a(1,i+1,j,k)-a(1,i,j,k))*dy(j)*dz(k) &
                   +(a(2,i,j+1,k)-a(2,i,j,k))*dx(i)*dz(k)
    forall(i=is:ie , j=js:je , k=ks:ke, ndims.eq.3) &
         b(i,j,k) = b(i,j,k)+(a(3,i,j,k+1)-a(3,i,j,k))*dx(i)*dy(j)
  end subroutine grid_divergence
!
! -------------------------------------------------------
!
! -- Get grid coeffs for Poisson solver
!
  pure subroutine grid_poisson(d,a)
    use global
    implicit none
    integer,intent(in) :: d
    real(8),intent(out),dimension(ni,nj,nk) :: a
    integer :: i,j,k
    a = 0
    forall(i=2:ni , j=1:nj , k=1:nk , d.eq.1) &
         a(i,j,k) = 2./(dx(i)+dx(i-1))*dy(j)*dz(k)
    forall(i=1:ni , j=2:nj , k=1:nk , d.eq.2) &
         a(i,j,k) = 2./(dy(j)+dy(j-1))*dx(i)*dz(k)
    forall(i=1:ni , j=1:nj , k=2:nk , d.eq.3) &
         a(i,j,k) = 2./(dz(k)+dz(k-1))*dx(i)*dy(j)
  end subroutine grid_poisson
!
! -------------------------------------------------------
!
! -- Plot out cell-centered grid
!
  subroutine grid_print(file,s,iss,iee,jss,jee,kss,kee)
    use global
    use utility, only: shift
    implicit none
    integer,intent(in) :: file,s,iss,iee,jss,jee,kss,kee
    integer :: i,j,k
    real(8) :: xcen(ni),ycen(nj),zcen(nk)

    xcen = xpos+dx*0.5
    ycen = ypos+dy*0.5
    zcen = zpos+dz*0.5

    write(file,3) (((xcen(i),i=iss,iee,s),j=jss,jee,s),k=kss,kee,s)
    write(file,3) (((ycen(j),i=iss,iee,s),j=jss,jee,s),k=kss,kee,s)
    if(ndims.eq.3) &
    write(file,3) (((zcen(k),i=iss,iee,s),j=jss,jee,s),k=kss,kee,s)
3   format(15e14.6)
  
  end subroutine grid_print
!
!----------------------------------------------------------------
end module grid
