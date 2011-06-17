!------------------------------------------------------------!
!--------------- Velocity Field Routines --------------------!
!------------------------------------------------------------!
module velo
  use global
  use inout
  private
  public init_velo,velo_RHS,velo_update
  real(8)           :: Rei=0
  logical           :: v_limit=.true.
  type(print_flags) :: pflags
contains
!-------------------------------------------------------------
  subroutine init_velo(u)
    use geom
    use grid
    implicit none
    real(8),intent(inout) :: u(ndims,ni,nj,nk)
    integer :: d,i,j,k

    if(thk.lt.2) call io_error('init_quick: wrong ghost thickness')
    call io_read
    call io_read(r=Rei)
    call io_read(l=v_limit)
    call io_read(f=pflags)
    write(io_log,'("  Rei=",e12.4," limit?=",l)') Rei,v_limit
    call log_print   
    
    if(.not.lres) then
       forall(d=1:ndims,i=1:ni,j=1:nj,k=1:nk) &
            u(d,i,j,k) = geom_velo(d,grid_pos(d,i,j,k))
    end if
    call io_write(0,pflags,u=u)

  end subroutine init_velo
!
! -------------------------------------------------------
!
  subroutine velo_RHS(u,R)
    use freeint
    implicit none
    real(8),intent(in)    :: u(ndims,ni,nj,nk)
    real(8),intent(inout) :: R(ndims,ni,nj,nk)
    real(8) :: v(ndims,ni,nj,nk),rho(ndims,ni,nj,nk)

    call convect(u,v)
    R = R-v
    if(Rei.ne.0) then
       if(freeint_on) then
          call freeint_get(density=rho,viscosity=v)
          call diffuse_mu(u,v)    ! v is in/out
          R = R+Rei*v/rho
       else
          call diffuse(u,v)
          R = R+Rei*v
       end if
    end if

  end subroutine velo_RHS
!
! -------------------------------------------------------
!
  subroutine velo_update(t,m,u)
    use freeint
    use mympi
    use grid
    use inout, only: io_write
    implicit none
    integer,intent(in)    :: t,m
    real(8),intent(inout) :: u(ndims,ni,nj,nk)
    real(8) :: C(ni,nj,nk),dt_CFL

    call mympi_vector(u,all=.true.)
    if(m.eq.2) call io_write(t,pflags,u=u)

    call grid_cfl(u,C)
    if(freeint_on) C = 2.*C ! vof requires halfing the CFl conditions
    dt_CFL = 1./mympi_max(maxval(C))
    if(m.eq.2) then
       write(io_log,'(i6,e12.4)') t,dt_CFL
       call log_print(CFL_num)
    end if

    if(dt_CFL.lt.dt) then
       if(v_limit) then
          where(spread(C,1,ndims).gt.dti) u = u/spread(C,1,ndims)
       else
          C = -log(merge(dti,dti/C,C.eq.zero))
          call io_write(t,p=C,error=.true.)
          call io_error('velo_update: CFL')
       end if
    end if

!!$    C = -log(merge(dti,dti/C,C.eq.zero))
!!$    C = min(max(C,-9.99),9.99)
!!$    if(t.eq.1.and.m.eq.1) call io_write(0,pflags,C)
!!$    if(m.eq.2) call io_write(t,pflags,C)

  end subroutine velo_update
!
! -------------------------------------------------------
!
! -- Compute diffusion term
!
!
! -- use f_d = lap(u_d)
!
  subroutine diffuse(u,fdiff)
    implicit none
    real(8),intent(in), dimension(ndims,ni,nj,nk)   :: u
    real(8),intent(inout),dimension(ndims,ni,nj,nk) :: fdiff
    real(8),dimension(ni,nj,nk) :: v,w
    integer :: d,d2
    fdiff = 0
    component: do d=1,ndims
       derivative_direction: do d2=1,ndims
          w = u(d,:,:,:)
          call ddx_mom(d,d2,w,v,.true.)
          call ddx_mom(d,d2,v,w,.false.)
          fdiff(d,:,:,:) = fdiff(d,:,:,:)+w
       end do derivative_direction
    end do component
  end subroutine diffuse
!
! -- use f_d = ppl{mu*[u_{d,d2}+u_{d2,d}]}{x_{d2}}
!
  subroutine diffuse_mu(u,fdiff)
    use grid, only: grid_linear
    implicit none
    real(8),intent(in), dimension(ndims,ni,nj,nk)   :: u
    real(8),intent(inout),dimension(ndims,ni,nj,nk) :: fdiff
    real(8),dimension(ni,nj,nk) :: v,w,hold
    real(8)                     :: mu(ndims,ni,nj,nk)
    integer                     :: d,d2
    mu = fdiff
    fdiff = 0
    component: do d=1,ndims
       derivative_direction: do d2=1,ndims
          call grid_linear(d,d2,mu(d,:,:,:),hold)
          call ddx_mom(d,d2,u(d ,:,:,:),v,.true.)
          call ddx_mom(d2,d,u(d2,:,:,:),w,.true.)
          v = hold*(v+w)
          call ddx_mom(d,d2,v,w,.false.)
          fdiff(d,:,:,:) = fdiff(d,:,:,:)+w
       end do derivative_direction
    end do component
  end subroutine diffuse_mu
!
! -------------------------------------------------------
!
! -- Compute convective force using QUICK scheme
!
  subroutine convect(u,fconv)
    use grid, only: grid_linear
    implicit none
    real(8),intent(in), dimension(ndims,ni,nj,nk) :: u
    real(8),intent(out),dimension(ndims,ni,nj,nk) :: fconv
    real(8),dimension(ni,nj,nk) :: mf,v
    integer :: d,d2
    fconv = 0
!
! -- loop through vector components (d) and faces (d2)
!
    do d=1,ndims
       do d2 = 1,ndims
!
! -- interpolate mass flux and velocity to face
!
          call grid_linear(d2,d,u(d2,:,:,:),mf)
          call quick(d2,mf,u(d ,:,:,:),v )
!
! -- get gradient of velocity flux, save
!
          v = mf*v
          call ddx_mom(d,d2,v,mf,.false.)
          fconv(d,:,:,:) = fconv(d,:,:,:)+mf
       end do
    end do

  end subroutine convect
!
! -------------------------------------------------------
!
! -- Derivatives for the momentum control volume.
! Note: 'p2m' goes from pressure CV to mom CV, else reverse
! Note: the inline derivative is shifted
!    
  subroutine ddx_mom(d,d2,a,b,p2m)
    use utility, only: shift
    use grid,    only: ddx_face,ddx_cent
    implicit none
    integer,intent(in) :: d,d2
    real(8),dimension(ni,nj,nk),intent(in)  :: a
    real(8),dimension(ni,nj,nk),intent(out) :: b
    logical,intent(in) :: p2m
    real(8),dimension(ni,nj,nk)             :: c
    if(p2m) then
       if(d.eq.d2) then
          call ddx_face(d2,a,c)
          call shift(-1,d2,c,b)
       else
          call ddx_cent(d2,a,b)
       end if
    else
       if(d.eq.d2) then
          call shift( 1,d2,a,c)
          call ddx_cent(d2,c,b)
       else
          call ddx_face(d2,a,b)
       end if
    end if
  end subroutine ddx_mom
!
! -------------------------------------------------------
!
! -- Flux-limited upwind quadratic interpolation
!
! This routine sets up pointers (including grid metrics)
! that are passed to the elemental algorithm below
!
  subroutine quick(d,mf,v,q)
    use grid,    only: grid_spacing_array
    use utility, only: pshift
    implicit none
    integer,intent(in) :: d
    real(8),intent(in),dimension(ni,nj,nk) :: mf,v
    real(8),intent(out),dimension(ni,nj,nk) :: q
    real(8),dimension(ni,nj,nk) :: di
    real(8),pointer,dimension(:,:,:) :: vp1,v0,vm1,vm2
    real(8),pointer,dimension(:,:,:) :: mf0,q0,dm2,dm1,d0
!
! -- get the inverse grid spacing coeffs and initalize q to 0
    di = 1./grid_spacing_array(d)
    q = 0
!
! -- point to neighboring values of v and di
    vp1 => pshift(v , 1,d)  ! plus one
    vm1 => pshift(v ,-1,d)  ! minus one
    vm2 => pshift(v ,-2,d)  ! minus two
    dm2 => pshift(di,-2,d)  ! minus two
    dm1 => pshift(di,-1,d)  ! minus one
!
! -- create unshifted pointers of the same size
    v0  => pshift(v ,0,d)
    q0  => pshift(q ,0,d)
    mf0 => pshift(mf,0,d)
    d0  => pshift(di,0,d)
!
! -- call elemental procedure
    call quick_elemental(mf0,vp1,v0,vm1,vm2,dm2,dm1,d0,q0)
  end subroutine quick
!
! -------------------------------------------------------
!
  elemental subroutine quick_elemental &
       (mf,vp1,v,vm1,vm2,dm2,dm1,d0,q)
    implicit none
    real(8),intent(in) :: mf,vp1,v,vm1,vm2,dm2,dm1,d0
    real(8),intent(out) :: q
    real(8) :: u,c,d,a,b

    if(mf.gt.0) then
       u = vm2
       c = vm1
       d = v
       a = dm2
    else
       u = vp1
       c = v
       d = vm1
       a = d0
    end if
    a = a/dm1
    b = 2.*a/(1.+a)
    q = (-a*b*u+(4.+2.*a)*c+(4.-b)*d)*0.125
    call limit(q,d,c,u)

  end subroutine quick_elemental
!
! -- Limit the slope
!
  elemental subroutine limit(funq,fund,func,funu)
    implicit none
    real(8),intent(in) :: fund,func,funu 
    real(8),intent(inout) :: funq
    real(8) :: ref1,ref2

    ref1=funu+10*(func-funu)
    call median(ref2,func,fund,ref1)
    call median(ref1,funq,func,ref2)
    funq=ref1

  end subroutine limit
!
! -- Find the median of three values
!
  elemental subroutine median(fun,a,b,c)
    implicit none
    real(8),intent(in) :: a,b,c
    real(8),intent(out) :: fun
    fun=max(min(a,b),min(max(a,b),c))
  end subroutine median
!
! -------------------------------------------------------
end module velo
