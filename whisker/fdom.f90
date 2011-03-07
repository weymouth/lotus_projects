!------------------------------------------------------------!
!---------------- Physical Domain Routines ------------------!
!------------------------------------------------------------!
module domain
  use global
  private
  real(8) :: U_max=0,grav=0,ramp_time=0,aref,vref
!
! Note: u,p boundary types (btypes)
!   0 -> imposed reference-frame velocity U=(vref,0,0)
!   1 -> no adjustment (BDIM or homogeneous pressure, P = 0)
!   2 -> set in geom_domain (as above, defaul to P=0)
!   3 -> u(x) = u0(x-dx) which is a convection condition with c=dx/dt
!   4 -> 3 with global flux correction
!
  integer,allocatable :: btypes(:),btypes4mgsolver(:)
  public init_domain,domain_update,domain_acceleration,btypes4mgsolver
  interface pdomain
     module procedure pdomain_s,pdomain_v
  end interface
contains
!----------------------------------------------------------------
!
! -- Initialize the domain info
!
  subroutine init_domain
    use inout
    use mympi, only: mympi_periodic
    implicit none
    integer :: d
    allocate(btypes(-ndims:ndims),btypes4mgsolver(-ndims:ndims))

    call io_read
    call io_read(ia=btypes)
    if(any(btypes.lt.0.or.btypes.gt.4)) &
         call io_error('init_domain: not a valid domain bytpe')
    write(io_log,'("  domain btypes=",7i2)') btypes
    call log_print

    call io_read(r=U_max)
    call io_read(r=ramp_time)
    call io_read(r=grav)
    write(io_log,'("  U,T,g=",3e12.4)') U_max,ramp_time,grav
    call log_print
!
! -- mgsovler uses a difference btype convention, translate to:
!  0 = imposed vectors, 1 = zeroed scalars, 2 = periodic
    btypes4mgsolver = btypes ! for 0,1
    where(btypes.gt.1) btypes4mgsolver = 1
    where(btypes.gt.2) btypes4mgsolver = 0
    forall(d=1:ndims,mympi_periodic(d))
       btypes4mgsolver( d) = 2
       btypes4mgsolver(-d) = 2
    end forall
  end subroutine init_domain
!
  function domain_btypes()
    implicit none
    integer :: domain_btypes(-ndims:ndims)
    domain_btypes = btypes
  end function domain_btypes
!
!----------------------------------------------------------------
!
! -- Domain Update
!
  subroutine domain_update(t,u0,del,ub)
    use mympi, only: mympi_domain_bound
    use grid,  only: grid_position_array
    use geom,  only: geom_domain,geom_domain_set ! <- content here
    implicit none
    integer,intent(in)                              :: t
    real(8),intent(in),dimension(ndims,ni,nj,nk)    :: u0
    real(8),intent(inout),dimension(ndims,ni,nj,nk) :: del,ub
    real(8),dimension(ni,nj,nk),target :: xt,yt,zt
    real(8),pointer,dimension(:,:)     :: pdel,pub,x,y,z,pu0
    real(8)                            :: time,a
    integer                            :: d
!
! -- Update reference-frame velocity and acceleration
    time = t*dt+time0
    if(ramp_time.eq.zero) then
       vref = U_max
       aref = zero
    else
       a = time/ramp_time
       vref = U_max*merge(3.*a**2-2.*a**3,one,a.lt.one)
       a = max(a-dt/ramp_time,zero)
       aref = U_max*merge(3.*a**2-2.*a**3,one,a.lt.one)
       aref = (vref-aref)*dti
    end if
!
! -- Update domain values of del and ub
    do d=-ndims,ndims
       if(.not.mympi_domain_bound(d)) cycle
       pdel => pdomain(d,del,0)
       pub  => pdomain(d,ub ,0)
!
! -- velocity bounds
       if(btypes(d).eq.0) then
          pdel = 1
          x => pdomain(d,del,1)
          pub  = (1.-x)*merge(vref,zero,abs(d).eq.1)
          if(d.eq.-1.and.time.lt.ramp_time) &
               call wavy_inlet_hack ! note temp
!
! -- geom bounds
       else if(btypes(d).eq.2) then
          xt = grid_position_array(1,d)
          yt = grid_position_array(2,d)
          zt = grid_position_array(3,d)
          x => pdomain(d,xt,0)
          y => pdomain(d,yt,0)
          z => pdomain(d,zt,0)
          call geom_domain_set(d)
          call geom_domain(x,y,z,pdel,pub)
!
! -- zero-grad
       else if(btypes(d).eq.3.or.btypes(d).eq.4) then
          pdel = 1
          pu0 => pdomain(d,u0,1)
          pub = pu0
          if(d.eq.1) call upwind_hack ! note temp ?
          if(btypes(d).eq.4) then ! correct for global flux
             x => pdomain(d,del,1)
             call domain_correct(d,merge(vref,zero,abs(d).eq.1),x,pub)
          end if
       end if
    end do
  contains
    subroutine wavy_inlet_hack
      use grid
      yt = grid_position_array(2,d)/grid_len(2)*2.*pi
      zt = grid_position_array(3,d)/grid_len(3)*2.*pi
      y => pdomain(d,yt,0)
      z => pdomain(d,zt,0)
      pub = pub*(1.D0+0.001*cos(z*5))
    end subroutine wavy_inlet_hack
    subroutine upwind_hack
      use grid
      real(8),pointer :: uU(:,:),uC(:,:),dx(:,:)
      xt = grid_spacing_array(d)
      dx => pdomain(d,xt,0)
      uU => pdomain(d,u0,1)         ! upwind
      uC => pdomain(d,u0,0)         ! central
      pub = uC-dt*U_max*(uC-uU)/dx  ! 1D convection equation
      if(d.eq.1) pub = pub+dt*aref  ! add body forces
      if(d.eq.ndims) pub = pub-dt*grav
    end subroutine upwind_hack
  end subroutine domain_update
!
!----------------------------------------------------------------
!
  subroutine domain_acceleration(R)
    implicit none
    real(8),intent(out) :: R(ndims,ni,nj,nk)
    R = 0
    R(1,:,:,:) = aref
    R(ndims,:,:,:) = -grav
  end subroutine domain_acceleration
!
!----------------------------------------------------------------
!
  subroutine domain_correct(d,uave,del,ub)
    use mympi, only: mympi_domain_sum
    use grid,  only: grid_area_array
    implicit none
    integer,intent(in)    :: d
    real(8),intent(in)    :: uave,del(:,:)
    real(8),intent(inout) :: ub(:,:)
    real(8)               :: da(ni,nj,nk),mine,area,goal,flux
    real(8),pointer       :: pda(:,:)
    da  = grid_area_array(abs(d))
    pda => pdomain(d,da,0)
    mine = sum(pda*(1-del))
    call mympi_domain_sum(abs(d),mine,area)
    mine = sum(pda*ub)
    call mympi_domain_sum(abs(d),mine,flux)
    mine = sum(pda*uave)
    call mympi_domain_sum(abs(d),mine,goal)
    ub = ub+(goal-flux)*(1-del)/area
  end subroutine domain_correct
!
!----------------------------------------------------------------
!
! -- Point to domain bound 'd', with shift into the domain 's'
!
  function pdomain_v(d,v,s)
    integer,intent(in)        :: d,s
    real(8),intent(in),target :: v(ndims,ni,nj,nk)
    real(8),pointer           :: pdomain_v(:,:)

    if     (d.eq.-1) then
       pdomain_v => v(1,is+s  ,js:je,ks:ke)
    else if(d.eq. 1) then
       pdomain_v => v(1,ie+1-s,js:je,ks:ke)
    else if(d.eq.-2) then
       pdomain_v => v(2,is:ie,js+s  ,ks:ke)
    else if(d.eq. 2) then
       pdomain_v => v(2,is:ie,je+1-s,ks:ke)
    else if(d.eq.-3) then
       pdomain_v => v(3,is:ie,js:je,ks+s  )
    else if(d.eq. 3) then
       pdomain_v => v(3,is:ie,js:je,ke+1-s)
    end if
    
  end function pdomain_v
  function pdomain_s(d,p,s)
    integer,intent(in)        :: d,s
    real(8),intent(in),target :: p(ni,nj,nk)
    real(8),pointer           :: pdomain_s(:,:)

    if     (d.eq.-1) then
       pdomain_s => p(is+s  ,js:je,ks:ke)
    else if(d.eq. 1) then
       pdomain_s => p(ie+1-s,js:je,ks:ke)
    else if(d.eq.-2) then
       pdomain_s => p(is:ie,js+s  ,ks:ke)
    else if(d.eq. 2) then
       pdomain_s => p(is:ie,je+1-s,ks:ke)
    else if(d.eq.-3) then
       pdomain_s => p(is:ie,js:je,ks+s  )
    else if(d.eq. 3) then
       pdomain_s => p(is:ie,js:je,ke+1-s)
    end if
    
  end function pdomain_s
!
!----------------------------------------------------------------
end module domain
