!------------------------------------------------------------!
!---------------- Volume of Fluid Routines ------------------!
!------------------------------------------------------------!
!
! This module holds the free interface data and routines. The
!  primary function is to supply the unsteady density field
!  to the momentum solver. This is done using VOF interface
!  tracking of the volume fraction f, a descrete field having
!  it's value defined as the fraction of each cell's volume 
!  filled with `dark' fluid.
! The module is parralel and dependant on the background and 
!  geom and blob modules.
!
! VOF boundary types (btypes)
! 0 => interface symmetry, df/dn=0
! 1 => fixed interface   , f = f0
!
module freeint
  use global
  use inout
  private
  public init_freeint,freeint_update,freeint_get
  real(8),parameter :: tol=1e-9,lam_rho=1.207e-3,lam_mu=1.863e-2
  real(8),dimension(:,:,:) ,allocatable         :: f0,alpha,fs
  real(8),dimension(:,:,:,:),allocatable,target :: norm
  logical,public      :: freeint_on=.false.
  integer,allocatable :: btypes(:)
  type(print_flags)   :: pflags
  integer             :: res_num = 2000,cutoff=0
contains
!-------------------------------------------------------------
!
! -- Set the starting interface
!
  subroutine init_freeint
    use mympi, only: mympi_read,mympi_id
    implicit none
    integer :: d
    real(8) :: t0
    allocate(btypes(-ndims:ndims))
    btypes = 0
!
! -- Check for free interface input info
    call io_read
    call io_read(l=freeint_on)  ! free interface activation flag
    call io_read(f=pflags)      ! print flags
    call io_read(ia=btypes)     ! boundary types (see above)
    call io_read(i=cutoff)      ! 'filter width' in grid cells
    if(.not.freeint_on) return
    write(io_log,'("  vof btypes=",7i2)') btypes
    call log_print
!
! -- Initialize f0,norm,alpha
    res_num = res_num+mympi_id()
    allocate(f0(ni,nj,nk),norm(3,ni,nj,nk),alpha(ni,nj,nk),fs(ni,nj,nk))
    call set_geom(f0)
    if(lres) then
       read(res_num) t0,fs
       if(t0.ne.time0) call io_error('init_freeint: timing mismatch')       
       call vof_domain_bounds(fs,all=.true.)
       f0 = fs
    end if    
    fs = f0
    call io_write(0,pflags,f0)

  end subroutine init_freeint
!
! -------------------------------------------------------
!
! -- Set the volume fraction using geom_fint
!
! geom_fint returns the interface distance function and 
!  interface normal given a point in space. at the cell 
!  corner, the distance is the negative intercept. From this
!  f{i,j,k} is given by Scardovelli below.
!
  subroutine set_geom(f)
    use grid, only: grid_corner,grid_vof_scale
    use geom, only: geom_fint  !<- content here
    implicit none
    real(8),intent(out) :: f(ni,nj,nk)
    real(8) :: fnc,n(3),a
    integer :: i,j,k

    do i=1,ni
    do j=1,nj
    do k=1,nk
       call geom_fint(grid_corner(i,j,k),a,n) ! get distance,normal
       call grid_vof_scale(i,j,k,n)           ! scale norm by dx
       f(i,j,k) = vof_vol(n(1),n(2),n(3),-a)  ! get f using Scardovelli
    end do
    end do
    end do

  end subroutine set_geom
!
! -------------------------------------------------------
!
! -- Get density/viscosity field at the faces
!
! \rho{x} = f{x}(1-\lamda_\rho)+\lamda_\rho
! \mu{x}  = f{x}(1-\lamda_\mu )+\lamda_\mu
!
  subroutine freeint_get(density,viscosity,vof)
    use mympi
    use grid
    implicit none
    real(8),dimension(ndims,ni,nj,nk),intent(out),optional :: density,viscosity
    real(8),dimension(ni,nj,nk),intent(out),optional       :: vof
    real(8) :: fsmooth(ni,nj,nk),fface(ndims,ni,nj,nk)
    integer :: d

    if(.not.freeint_on) then
       if(present(density)) density = 1
       if(present(viscosity)) viscosity = 1
       if(present(vof)) vof = 1
       return
    end if
    if(present(vof)) then
       vof = f0
    end if
    if(present(density).or.present(viscosity)) then
       call vof_smooth(2,fs,fsmooth)
       do d=1,ndims
          call grid_linear(0,d,fsmooth,fface(d,:,:,:))
       end do
       call mympi_vector(fface)
       if(present(density)) then
          density = fface*(one-lam_rho)+lam_rho
       end if
       if(present(viscosity)) then
          viscosity = fface*(one-lam_mu)+lam_mu
       end if
    end if

  end subroutine freeint_get
!
! -------------------------------------------------------
!
! -- Smooth the interface transition over a few cells
!
  subroutine vof_smooth(it,f,sf)
    use mympi
    implicit none
    integer :: it
    real(8),dimension(ni,nj,nk),intent(in) :: f
    real(8),dimension(ni,nj,nk),intent(out) :: sf
    integer :: d,i,j,k
    sf = f
    do d = 1,it
!!$       forall(i=is:ie , j=js:je , k=ks-ndims+2:ke+ndims-2) 
!!$          sf(i,j,k) = 4.*sf(i,j,k)+sf(i-1,j,k)+sf(i+1,j,k) &
!!$                                  +sf(i,j-1,k)+sf(i,j+1,k)
!!$       end forall
!!$       sf = sf*0.125
!!$       if(ndims.eq.3) then
!
! note: this is a uniform filter, whereas a triangular is used above
!
!!$          forall(i=is:ie , j=js:je , k=ks:ke) 
!!$             sf(i,j,k) = sf(i,j,k)+sf(i,j,k-1)+sf(i,j,k+1)
!!$          end forall
!!$          sf = sf/3.
!!$       end if
       forall(i=is:ie , j=js:je , k=ks:ke) 
          sf(i,j,k) = 6.*sf(i,j,k)+sf(i-1,j,k)+sf(i+1,j,k) &
               +sf(i,j-1,k)+sf(i,j+1,k)+sf(i,j,k-1)+sf(i,j,k+1)
       end forall
       sf = sf/12.
       call mympi_scalar(sf)
    end do
  end subroutine vof_smooth
!
! -------------------------------------------------------
!
! -- Update f with my operator-split advection method
!
!  f = f+\dt(\ddx{F{d},d}+g\ddx{u{d},d}) for d=1,ndims
!
! Where g = merge(0,1,f0.le..5), the cell center value 
!  of the color function (since the interface is linear
!  if the cell is more than half full the center is dark),
!  and F is the flux of dark fluid through the cell face.
!
  subroutine freeint_update(t,m,u,u0)
    use mympi, only: mympi_scalar,mympi_id
    use grid,  only: ddx_face
    use blob,  only: blob_find
    implicit none
    integer :: t,m,d,d2,d3
    real(8),dimension(ndims,ni,nj,nk),intent(in) :: u,u0
    real(8),dimension(ni,nj,nk) :: v,ddx_v,flux
    integer,dimension(ni,nj,nk) :: g0
    if(.not.freeint_on) return
!
! -- reset f and construct g
!
    fs = f0
    g0 = merge(0,1,f0.le..5)
!
! -- loop through the split-operator method
!
    dir_loop: do d2=1,ndims
       d = mod(t+d2,ndims)+1
!
! -- reconstruct the interface
!
       call vof_reconstruct(fs)
!
! -- get scaled velocity and it's derivative
!
       v = (u(d,:,:,:)+u0(d,:,:,:))*.5*dt
       where(abs(v).lt.tol) v = 0.
       call ddx_face(d,v,ddx_v)
!
! -- get flux
!
       call vof_flux(d,fs,v,flux)
!
! -- update volume fraction with net flux and dilitation
!
       call ddx_face(d,flux,v)
       fs = fs-v+g0*ddx_v
!
! -- clean, update the boundary values
!
       where(   fs.lt.tol) fs = 0.
       where(1.-fs.lt.tol) fs = 1.
       call mympi_scalar(fs,all=m.eq.2.and.d2.eq.ndims)
       call vof_domain_bounds(fs,all=m.eq.2.and.d2.eq.ndims)
    end do dir_loop
!
! -- clean up underresolved blobs 
!
! note: this should be part of an "SGS" model. the mass removed
! note:  here should be added to a secondary SGS VOF field with
! note:  its own governing dynamics. momentum should be adjusted too.
!
    if(m.eq.2) then
       if(cutoff.gt.0) then
          call blob_find(fs,cutoff)
          v = 1.-fs
          call blob_find(v,cutoff)
          fs = 1.-v
          call mympi_scalar(fs,all=.true.)
          call vof_domain_bounds(fs,all=.true.)
       end if
!
! -- Save the new solution
!
       f0 = fs
       call io_write(t,pflags,f0)
       if(mod(t,tmd).eq.0) then
          rewind(res_num)
          write(res_num) time0+t*dt,f0
          call flush(res_num)
       end if
    end if
  end subroutine freeint_update
!
! -------------------------------------------------------
!
! -- Reconstruct the interface using my method
!
! The interface is defined by alpha=sum{norm*x}, where 
!  alpha is the intercept. From this we can write
!
!     norm{d}  = -\ddx{height,d} for d=1,ndims; d.ne.dc
!     norm{dc} = 1
!
!  where, the interface height is the distance above 
!  the cell bottom in the dc direction and dc is chosen
!  as the dominant normal direction.
! A second order estimate of is height = f\dx{dc} if
!  the interface doesn't pass through the top or bottom
!  cell walls. We estimate dc using norm\sim-\grad{f}.
! With norm known, alpha is given by Scardovelli below.
!
  subroutine vof_reconstruct(f)
    use mympi
    use grid
    use utility, only: ishift,get
    implicit none
    real(8),dimension(ni,nj,nk),intent(in)      :: f
    real(8),pointer :: nc(:)
    integer :: i,j,k,d,dc,p(4),pu(4),pd(4)
    real(8) :: fc,hu,hd,hc,lu,ld,lc,xu,xd,xc
    norm = 0
!
! -- loop through points
    do i=is,ie
    do j=js,je
    grid_loop: do k=ks,ke
!
! -- check for interface
       fc = f(i,j,k)
       if(fc.eq.zero.or.fc.eq.1) then
          alpha(i,j,k) = fc
          cycle grid_loop
       end if
!
! -- first guess, n = -\grad{f}
       nc => norm(:,i,j,k)
       do d = 1,ndims
          p = (/d,i,j,k/)
          pu = ishift(p, 1,d)
          pd = ishift(p,-1,d)
          nc(d) = -(get(f,pu)-get(f,pd)) &
               &/(grid_size(pu)+grid_size(pd))
       end do
       call grid_vof_scale(i,j,k,nc)
!
! -- get the dominate direction
       dc = sum(maxloc(abs(nc)))
       if(dc.gt.ndims) then
          write(io_log,"(4i4,4e10.2)") mympi_id(),i,j,k,nc,fc
          call log_print
          dc = 1
       end if
!
! -- refine with n(d) = -\ddx{height,d}
       p = (/dc,i,j,k/)
       direction_loop: do d=1,ndims
          if(d.eq.dc) then
             nc(d) = sign(one,nc(d))
          else
!
! -- get info
             pu = ishift(p, 1,d)
             pd = ishift(p,-1,d)
             call vof_height(pu,f,hu,lu)
             call vof_height(pd,f,hd,ld)
             call vof_height(p ,f,hc,lc)
             xu = grid_size(ishift(pu,d,0))
             xd = grid_size(ishift(pd,d,0))
             xc = grid_size((/d,i,j,k/))
!
! -- central difference
             nc(d) = -(hu-hd)/(xc+0.5*(xu+xd))
!
! -- if boundary, use one-sided difference
             if(mympi_out(pu)) then
                nc(d) = -(hc-hd)*2./(xc+xd)
             else if(mympi_out(pd)) then
                nc(d) = -(hu-hc)*2./(xu+xc)
!
! -- if completely full or empty nearby; unresolved.
             else if(hu+hd.eq.0.or.hu+hd.eq.lu+ld) then
                nc = 0
                cycle grid_loop
!
! -- if too steep, the interface will pass through the 
! --  top or bottom of the u or d cells.
! -- use a one-sided difference to avoid referencing it
             else if(abs(nc(d)).gt.0.5) then
                if(nc(d)*(fc-0.5).ge.zero) then
                   nc(d) = -(hu-hc)*2./(xu+xc)
                else
                   nc(d) = -(hc-hd)*2./(xc+xd)
                end if
             end if
          end if
       end do direction_loop
!
! -- scale n by \dx to effectively give the cell unit length
       call grid_vof_scale(i,j,k,nc)
!
! -- get intercept using Scardovelli
       alpha(i,j,k) = vof_int(nc(1),nc(2),nc(3),fc)
    end do grid_loop
    end do
    end do
!
! -- set BVs
    do d=1,ndims
       call mympi_scalar(norm(d,:,:,:))
    end do
    call mympi_scalar(alpha)
    call vof_domain_bounds(fin=f)
  end subroutine vof_reconstruct
!
! -------------------------------------------------------
!
! -- Find the mean surface height and length of the region
! by summing over three cells in the dc direction
!
  subroutine vof_height(p,f,h,l)
    use global
    use mympi,   only: mympi_out
    use grid,    only: grid_size
    use utility, only: ishift,get
    implicit none
    integer,intent(in) :: p(4)
    real(8),dimension(ni,nj,nk),intent(in) :: f
    real(8),intent(out) :: h,l
    integer :: p2(4),j
    h = 0; l = 0
    p2 = ishift(p,-1,p(1))

    do j=1,3
       if(.not.mympi_out(p2)) then
          h = h+grid_size(p2)*get(f,p2) ! h=sum(f\dx)
          l = l+grid_size(p2)           ! l=sum(\dx)
       end if
       p2 = ishift(p2,1,p(1))
    end do
  end subroutine vof_height
!
! -------------------------------------------------------
!
! -- Compute dark volume fluxed through face of a unit cell
!
! The flux is computed geometrically using a donor-region method.
!  The donating region is in the cell upwind of the face. The 
!  width of the nondimensional region is \dl=v\dt/\dx. v is scaled
!  by \dt on input and flux is scaled by \dt on output
!
  subroutine vof_flux(d,f_in,v_in,flux)
    use global
    use grid

    use mympi

    use utility, only: ishift,get
    implicit none
    integer,intent(in) :: d
    real(8),dimension(ni,nj,nk),intent(in) :: f_in
    real(8),dimension(ni,nj,nk),intent(in) :: v_in
    real(8),dimension(ni,nj,nk),intent(out) :: flux
    real(8),target :: n_tar(3)
    real(8),pointer :: n1,n2,n3,nd
    real(8) :: v,f,dl,a
    integer :: i,j,k,p(4),plus(3)
    plus = 0; plus(d) = 1
!
! -- loop through cell faces and get velocity
!    
    flux = 0
    do i=is,ie+plus(1)
    do j=js,je+plus(2)
    do k=ks,ke+plus(3)
       v = v_in(i,j,k)
       if(v.eq.zero) cycle
!
! -- get upwind cell and volume fraction there
!     positive v means the cell is one step back
!
       p = (/d,i,j,k/)
       if(v.gt.zero) p = ishift(p,-1,d)
       f = get(f_in,p)
!
! -- for full cells, flux is v 
!
       if(f.eq.one) then
          flux(i,j,k) = v
!
! -- on interface, get normal vector and donation width
!
       else if(f.ne.zero) then
          n_tar = norm(:,p(2),p(3),p(4))
          if(sum(abs(n_tar)).eq.zero) cycle
          n1 => n_tar(1)
          n2 => n_tar(2)
          n3 => n_tar(3)
          nd => n_tar(d)
          dl = v/grid_size(p)
!
! -- get the intercept for the upwind cell
!     again, positive v requires a shift
!
          a = get(alpha,p)
          if(dl.gt.zero) a = a-nd*(1.-dl)
!
! -- scale the normal by the width, effectively unitizing 
!     the donating region, and compute volume using Scardovelli
!
          nd = abs(dl)*nd
          flux(i,j,k) = v*vof_vol(n1,n2,n3,a)
!!$
!!$          if(mympi_coords(1).eq.0.and.&
!!$               mympi_coords(2).eq.0.and.&
!!$               mympi_coords(3).eq.2) then
!!$             if(all(p.eq.(/3,10,5,nk-1/))) then
!!$                print *,'bound'
!!$                print '(3e12.4)',n_tar
!!$                print '(3e12.4)',nd,a,dl
!!$                print '(3e12.4)',v,flux(i,j,k)
!!$             else if(all(p.eq.(/3,10,5,nk-2/))) then
!!$                print *,'not bound'
!!$                print '(3e12.4)',n_tar
!!$                print '(3e12.4)',nd,a,dl
!!$                print '(3e12.4)',v,flux(i,j,k)
!!$             end if
!!$          end if
!!$          
       end if
    end do
    end do
    end do
  end subroutine vof_flux
!
! -------------------------------------------------------
!
! -- Set the domain boundary values of f or alpha,norm
!
  subroutine vof_domain_bounds(f,all,fin)
    use utility, only: pbound
    use mympi,   only: mympi_domain_bound
    implicit none
    real(8),intent(inout),optional   :: f(ni,nj,nk)
    logical,intent(in),optional      :: all
    real(8),intent(in),optional      :: fin(ni,nj,nk)
    real(8),pointer,dimension(:,:,:) :: get,put,n1,n2,n3
    integer                          :: d,d2
!
! -- loop through domain boundaries
!
    do d=-ndims,ndims
       if(.not.mympi_domain_bound(d)) cycle
!
! -- set f
!
       if(present(f)) then
          if(btypes(d).eq.1) then
             put => pbound(f0,d,all=all)
             get => pbound(f,d,all=all)
             get = put
          end if
       else
!
! -- set norm, alpha
!
          do d2=1,ndims
             get => pbound(norm(d2,:,:,:),d)             
! 
! -- set side wall fixed interface normal (vertical)
             if(btypes(d).eq.1.and.abs(d).lt.ndims) then
                get = 0
                if(d2.eq.ndims) get = 1
!
! -- set reflection interface normal 
             else if(btypes(d).eq.0) then
                if(d2.ne.ndims) get = -get
             end if
          end do
!
! -- set alpha based on new norm and f0
          n1 => pbound(norm(1,:,:,:),d)
          n2 => pbound(norm(2,:,:,:),d)
          n3 => pbound(norm(3,:,:,:),d)
          put => pbound(fin,d)
          get => pbound(alpha,d)
          get = vof_int(n1,n2,n3,put)
       end if
    end do

  end subroutine vof_domain_bounds
!
! -------------------------------------------------------
! -- See Scardovelli and Zaleski, JCP 164 (2000) for 
! details on the 3d linear interface geometric equations
! -------------------------------------------------------
!
  pure elemental function vof_int(n1,n2,n3,g)
    use global, only: one,zero
    implicit none
    real(8),intent(in) :: n1,n2,n3,g
    real(8) :: a,m1,m2,m3,t,vof_int
!
! -- scale and sort the slopes
!
    t = abs(n1)+abs(n2)+abs(n3)
    if(g.ne.0.5) then
       call sort3(abs(n1)/t,abs(n2)/t,abs(n3)/t,m1,m2,m3)
!
! -- get the intercept, scale and shift it 
!
       a = a3(m1,m2,m3,merge(g,1.-g,g.lt.0.5))
    else
       a = 0.5
    end if
100 vof_int = merge(a,1.-a,g.lt.0.5)*t &
         +min(n1,zero)+min(n2,zero)+min(n3,zero)

  end function vof_int
!
! -------------------------------------------------------
!
  pure elemental function vof_vol(n1,n2,n3,b)
    use global, only: one,zero
    implicit none
    real(8),intent(in) :: n1,n2,n3,b
    real(8) :: a,m1,m2,m3,t,vof_vol
!
! -- shift and scale the intercept
!    
    t = abs(n1)+abs(n2)+abs(n3)
    a = (b-min(n1,zero)-min(n2,zero)-min(n3,zero))/t
!
! -- check bounds
!
    if(a.le.zero.or.a.eq.0.5.or.a.ge.one) then
       vof_vol = min(max(a,zero),one)
    else
!
! -- scale and sort the slopes and get fraction
!
       call sort3(abs(n1)/t,abs(n2)/t,abs(n3)/t,m1,m2,m3)
       t = f3(m1,m2,m3,merge(a,1.-a,a.lt.0.5))
       vof_vol = merge(t,1.-t,a.lt.0.5)
    end if

  end function vof_vol
!
! -------------------------------------------------------
!
  pure elemental subroutine sort3(n1,n2,n3,m1,m2,m3)
    implicit none
    real(8),intent(in)  :: n1,n2,n3
    real(8),intent(out) :: m1,m2,m3
    integer :: l
    real(8) :: n(3),m(2)

    n = (/n1,n2,n3/)
    l = sum(minloc(n)) ; m1 = n(l)
    m(1:l-1) = n(1:l-1) ; m(l:2) = n(l+1:3)
    m2 = minval(m); m3 = maxval(m)

  end subroutine sort3
!
! -------------------------------------------------------
!
  pure elemental function a3(m1,m2,m3,v)
    implicit none
    real(8),intent(in) :: m1,m2,m3,v
    real(8) :: m12,p,v1,v2,v3,c0,c1,c2,c3,a3
    
    m12 = m1+m2
    p = 6.*m1*m2*m3
    v1 = m1**2/(6.*m2*m3+tol)
    v2 = v1+(m2-m1)/(2.*m3)
    v3 = merge((m3**2*(3.*m12-m3)+m1**2*(m1-3.*m3) &
         +m2**2*(m2-3.*m3))/p,m12*.5/m3,m3.lt.m12)
    if(v.lt.v1) then
       a3 = (p*v)**(1./3.)
    else if(v.lt.v2) then
       a3 = 0.5*(m1+sqrt(m1**2+8.*m2*m3*(v-v1)))
    else if(v.lt.v3) then
       c0 = m1**3+m2**3-p*v
       c1 = -3.*(m1**2+m2**2)
       c2 = 3.*m12 ; c3 = -1
       a3 = proot(c0,c1,c2,c3) ! defined below
    else if(m3.lt.m12) then
       c0 = m1**3+m2**3+m3**3-p*v
       c1 = -3.*(m1**2+m2**2+m3**2)
       c2 = 3 ; c3 = -2
       a3 = proot(c0,c1,c2,c3) ! defined below
    else
       a3 = m3*v+m12*.5
    end if

  end function a3
!
! -------------------------------------------------------
!
  pure elemental function f3(m1,m2,m3,a)
    implicit none
    real(8),intent(in) :: m1,m2,m3,a
    real(8) :: m12,f3

    m12 = m1+m2

    if(a.lt.m1) then
       f3 = a**3/(6.*m1*m2*m3)
    else if(a.lt.m2) then
       f3 = a*(a-m1)/(2.*m2*m3)+m1**2/(6.*m2*m3+tol)
    else if(a.lt.min(m3,m12)) then
       f3 = (a**2*(3.*m12-a)+m1**2*(m1-3.*a)+m2**2*(m2-3.*a)) &
            /(6*m1*m2*m3)
    else if(m3.lt.m12) then
       f3 = (a**2*(3.-2.*a)+m1**2*(m1-3.*a)+m2**2*(m2-3.*a) &
            +m3**2*(m3-3.*a))/(6*m1*m2*m3)
    else
       f3 = (2.*a-m12)/(2.*m3)
    end if

  end function f3
!
! -------------------------------------------------------
!
  pure elemental function proot(c0,c1,c2,c3)
    implicit none
    real(8),intent(in) :: c0,c1,c2,c3
    real(8) :: a0,a1,a2,p0,q0,t,proot

    a0 = c0/c3 ; a1 = c1/c3 ; a2 = c2/c3

    p0 = a1/3.-a2**2/9.
    q0 = (a1*a2-3.*a0)/6.-a2**3/27.
    t = acos(q0/sqrt(-p0**3))/3.

    proot = sqrt(-p0)*(sqrt(3.)*sin(t)-cos(t))-a2/3.

  end function proot
!
! -------------------------------------------------------
end module freeint
