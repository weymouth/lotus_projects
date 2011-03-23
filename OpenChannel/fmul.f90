!------------------------------------------------------------!
!------------- MIT Multigrid Solver Wrappers ----------------!
!------------------------------------------------------------!
!
! -- Solve the system Ax=b using mit_mgsolver library
!
module solver
  use mit_mgsolver_matrices, only: matrix
  type(matrix),target,private :: A
contains
!-------------------------------------------------------------
  subroutine init_solver
    use global,                only: ndims,thk,ni,nj,nk
    use mympi,                 only: mympi_comm
    use domain,                only: btypes4mgsolver
    use inout,                 only: press_num
    use mit_mgsolver_dims,     only: init_dims
    use mit_mgsolver_utility,  only: log_cutoff,log_file_num
    use mit_mgsolver_mympi,    only: init_mympi
    use mit_mgsolver_matrices, only: init_matrix
    use mit_mgsolver_PCG,      only: init_PCG
    use mit_mgsolver_MG,       only: init_MG
    implicit none
    real(8) :: tol = 1e-6
    integer :: mx = 50

    call log_file_num(press_num) ! set in init_inout
    call log_cutoff(2)   ! don't print when indent>2
    call init_dims(set_ndims=ndims,set_thk=thk)
    call init_mympi(set_btypes=btypes4mgsolver,set_COMM=mympi_COMM())
    call init_matrix(A,ni,nj,nk)
    call init_PCG(set_it_mx=3,set_tol=tol,debug=.true.)
    call init_MG(set_it_mx=mx,set_tol=tol,set_smoother='PCG  ',&
         set_prolongate='cell ',set_adjust=.true.,debug=.true.)

  end subroutine init_solver
!
! ---------------------------------------------------------
!
  subroutine solver_update_coeffs(beta)
    use global,                only: ndims,ni,nj,nk,thk
    use grid,                  only: grid_poisson
    use mympi,                 only: mympi_vector
    use mit_mgsolver_matrices, only: matrix_make_poisson
    implicit none
    integer :: d
    real(8),intent(in),dimension(ndims,ni,nj,nk) :: beta
    real(8),dimension(ndims,ni,nj,nk) :: coeff

    do d=1,ndims
       call grid_poisson(d,coeff(d,:,:,:))
    end do
    coeff = coeff*beta
    call mympi_vector(coeff)
    call matrix_make_poisson(A,coeff)

  end subroutine solver_update_coeffs
!
! ---------------------------------------------------------
!
  subroutine solver_invert(t,x,b)
    use global,                only: ni,nj,nk
    use mit_mgsolver_MG,       only: MG_solver
    use mit_mgsolver_PCG,      only: PCG_solver
    use mit_mgsolver_mympi,    only: mympi_scalar
    implicit none
    integer,intent(in)                        :: t
    real(8),intent(in),dimension(ni,nj,nk)    :: b
    real(8),intent(inout),dimension(ni,nj,nk) :: x    
    if(t.eq.0) then
       call mympi_scalar(x,pflag=1)
       !call pcg_solver(A,x,b,set_nonlin=.true.)
       !call pcg_solver(A,x,b,set_nonlin=.true.)
       !call pcg_solver(A,x,b,set_nonlin=.true.)
       call pcg_solver(A,x,b)
       call pcg_solver(A,x,b)
       call pcg_solver(A,x,b)
    end if
    call MG_solver(A,x,b)
  end subroutine solver_invert
end module solver
