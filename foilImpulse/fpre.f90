!------------------------------------------------------------!
!------------------- Pressure Routines ----------------------!
!------------------------------------------------------------!
module pressure
  use solver   ! matrix inversion variables and routines
  use global
  use inout
  private
  public init_pressure,pressure_update
  type(print_flags) :: pflags
contains
!-------------------------------------------------------------
  subroutine init_pressure(p)
    implicit none
    real(8) :: p(ni,nj,nk)

    call io_read
    call io_read(f=pflags)
    call io_write(0,pflags,p)
    call init_solver
    
  end subroutine init_pressure
!
! ------------------------------------------------------------
!
  subroutine pressure_update(t,m,alpha,beta,R,p)  ! temp change for direct-forcing test
    use mympi
    use grid
    use freeint
    use body
    implicit none
    integer,intent(in)    :: t,m
    real(8),intent(in)    :: alpha(ndims,ni,nj,nk),beta(ndims,ni,nj,nk)
    real(8),intent(inout) :: R(ndims,ni,nj,nk),p(ni,nj,nk)
    real(8) :: f(ni,nj,nk),v(ndims,ni,nj,nk)
!
! -- Get source
    call mympi_vector(R)
    call grid_divergence(R,f)
!
! -- Get coefficients and solve
    if(body_unsteady().or.freeint_on.or.t.eq.0) &
         call solver_update_coeffs(alpha*beta)  ! alpha*beta
!!$         call solver_update_coeffs(beta)  ! only beta
    call solver_invert(t,p,f)
!
! -- Apply pressure force to RHS
    call grid_gradient(p,v)
    R = R-alpha*beta*v  ! alpha*beta
!
! -- Pressure output
    if(m.eq.1) then
       call body_int_force(t*dt+time0,p)
       call io_write(t,pflags,p)
    end if

  end subroutine pressure_update
end module pressure
