!-------------------------------------------------------!
!---------------- Immersed Body Routines ---------------!
!-------------------------------------------------------!
module body
  use global
  use inout
  private
  public init_body,body_update,body_unsteady,body_int_force
  real(8),dimension(:,:,:,:),allocatable :: kernel
  real(8) :: eps = 0
  logical :: body_on=.false.,body_move=.false.,body_slip=.true.
  type(print_flags) :: pflags
contains
  subroutine init_body(del,ub)    
    implicit none
    real(8),dimension(ndims,ni,nj,nk),intent(out) :: del,ub

    del = 0; ub = 0
    call io_read
    call io_read(l=body_on)
    call io_read(r=eps)
    call io_read(f=pflags)
    call io_read(l=body_move)
    call io_read(l=body_slip)
    if(.not.body_on) return
    if(eps.eq.0)  call io_error('body_init: must set eps')
    write(io_log,'("  eps,unsteady,slip=",e12.4,2l)') eps,body_move,body_slip
    call log_print

    allocate(kernel(ndims,ni,nj,nk))
    call body_update(0,del,ub)

  end subroutine init_body
  logical function body_unsteady()
    body_unsteady = body_on.and.body_move
  end function body_unsteady
!
! -------------------------------------------------------
!
! -- Update geometry influence arrays
!
! Note: body boundary type (flag)
! 0 -> no immersed surface influence
! 1 -> immersed body,     R^{ndims}
! 2 -> immersed manifold, R^{ndims-1}
!
  subroutine body_update(t,del,ub)
    use grid, only: grid_pos,grid_vol
    use geom, only: geom_body
    use slip, only: slip_reset,slip_add,slip_sort
    implicit none
    integer,intent(in)                              :: t
    real(8),dimension(ndims,ni,nj,nk),intent(inout) :: del,ub
    real(8),dimension(ni,nj,nk) :: dis
    real(8),dimension(ndims)    :: velocity,normal
    real(8)                     :: distance,delta,kappa,time,dv,dx
    integer                     :: d,i,j,k,flag
!
! -- check if we need to update
    if((.not.body_on).or.((.not.body_move).and.t.ne.0)) return
    kernel = 0
    if(body_slip) call slip_reset
!
! -- loop through grid
    do i=1,ni
    do j=1,nj
    do k=1,nk
!
! -- get properties on each face
       do d=1,ndims
          call geom_body(grid_pos(d,i,j,k),distance,normal,velocity,flag)
          call smooth_switch(distance/eps,flag,delta)
          if(body_slip) call slip_add &                 ! velocity is inout !
               ((/d,i,j,k/),distance.gt.zero,delta,normal,velocity,flag)
          del(d,i,j,k) = delta
          ub(d,i,j,k) = velocity(d)
       end do
!
! -- get properties on cell-center
       call geom_body(grid_pos(0,i,j,k),distance,normal,velocity,flag,kappa)
       dis(i,j,k) = min(max(distance/eps,-9.99),9.99)
       if(flag.ne.0) then
!!$          dv = grid_vol(i,j,k)
!!$          dx = dv**(1./3.)
!!$          call surf_del(distance/dx-one,delta)
!!$          kernel(:,i,j,k) = normal*delta*dv/(dx+distance*kappa)
          call surf_del(distance/eps-one,delta)
          kernel(:,i,j,k) = normal*delta*grid_vol(i,j,k) &
               /eps/(one+distance*kappa)
       end if
    end do
    end do
    end do

    call io_write(t,pflags,dis)
    if(body_slip) call slip_sort

  end subroutine body_update
!
! -------------------------------------------------------
!
! -- Integrate force on the body
!
  subroutine body_int_force(time,p)
    use mympi, only: mympi_sum
    implicit none
    real(8),intent(in) :: time,p(ni,nj,nk)
    real(8) :: force(ndims),force_all(ndims),l
    integer :: d

    if(.not.body_on) return
    do d=1,ndims
       force(d) = sum(kernel(d,is:ie,js:je,ks:ke) &
                            *p(is:ie,js:je,ks:ke))
       call mympi_sum(force(d),force_all(d))
    end do

    write(io_log,'(4e15.7)') time,force_all
    call log_print(force_num)  ! set in init_fino

  end subroutine body_int_force
!
! -------------------------------------------------------
!
! -- Implement the smoothed switch function
!
  pure subroutine smooth_switch(distance,flag,delta)
    implicit none
    real(8),intent(in)  :: distance
    integer,intent(in)  :: flag
    real(8),intent(out) :: delta
    if(flag.eq.0) then             ! nothing
       delta = 0.
    else if(abs(flag).eq.1) then   ! solid
       if(flag.gt.0) then             ! no slip
          call solid_body_del(distance,delta)
       else if(distance.gt.zero) then ! slip region
          call surf_del(distance,delta)
       else                           ! internal
          delta = one
       end if
    else if(abs(flag).eq.2) then   ! surface
       call surf_del(distance,delta)
    end if
  end subroutine smooth_switch
!
! -- soild body ...
  pure subroutine solid_body_del(distance,delta)
    implicit none
    real(8),intent(in) :: distance
    real(8),intent(out) :: delta
    if(distance.gt.one) then
       delta = zero
    else if(distance.lt.-one) then
       delta = one
    else
       delta = (1.-sin(distance*pi*0.5))*.5
    end if
  end subroutine solid_body_del
!
! -- surface ...
  pure subroutine surf_del(distance,delta)
    implicit none
    real(8),intent(in) :: distance
    real(8),intent(out) :: delta
    if(abs(distance).gt.one) then
       delta = zero
    else
       delta = (1.+cos(distance*pi))*.5
       delta = delta/(delta+eps)
    end if
  end subroutine surf_del
!
! -------------------------------------------------------
end module body
