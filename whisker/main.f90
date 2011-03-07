! -------------------------------------------------------
! ---------- Boundary Data Immersion Method -------------
! -------------------------------------------------------
!
program BDIM
!
! -- Solve for the BDIM unsteady velocity and pressure field
!
  use global   ! holds global variables
  use mympi    ! my MPI variables and functions
  use inout    ! I/O variables and routines
  use geom     ! problem geometry variables and functions
  use domain   ! numerical domain variables and functions
  use body     ! immersed body variables and functions
  use slip     ! immersed surface slip conditions
  use freeint  ! free interface variables and functions
  use velo     ! velocity specific functions
  use pressure ! pressure specific functions
  use apriori  ! writes out data for apriori analysis
  implicit none
  integer             :: res_num=1000,t=0
  real(8),allocatable :: u(:,:,:,:),p(:,:,:),del(:,:,:,:),ub(:,:,:,:)

  call init_BDIM
  time_loop: do t=0,tmx-1
     call BDIM_update
  end do time_loop
  deallocate(u,p)
  call mympi_end

contains
!
! -------------------------------------------------------
!
! -- Initialize the program
! Note: the primary input is (num=input_num , file='inp.txt')
! and these routines read from it sequentially. No shuffling!
!
  subroutine init_BDIM
    implicit none
!
! -- read and init global,mympi,grid
    call init_io
!
! -- read and init problem geometry
    call init_geom
!
! -- init domain variables
    call init_domain
!
! -- read and init velocity and pressure
    allocate(u(ndims,ni,nj,nk),p(ni,nj,nk))
    time0 = 0; u = 0; p = 0
    res_num = res_num+mympi_id()
    if(lres) read(res_num) time0,u,p
    write(io_log,'("  time0=",e12.4)') time0
    call log_print
    call init_velo(u)
    call init_pressure(p)
!
! -- read and init immersed body
    allocate(del(ndims,ni,nj,nk),ub(ndims,ni,nj,nk))
    call geom_update(time0)
    call init_body(del,ub)
!
! -- read and init free interface
    call init_freeint
!
! -- finish up initilization
    call init_apriori(u,p)
    call io_profile('init ')
 
  end subroutine init_BDIM
!
! -------------------------------------------------------
!
! Integrate the BDIM governing equation in time
!
!   u+L{u} = u0+del(ub-u0)+(1-del)int{R-grad{p}/rho,t=0,dt}
!
!  where u is the velocity, u0=u{t=0}, R is the right hand side
!  term containing the fluid and gravitational accelerations,
!  p is the total mechanical pressure, rho is the fluid density,
!  del and ub are the indicator function and velocity of
!  the immersed body, and L{u} is a slip model active on the 
!  solid/fluid interface.
!
  subroutine BDIM_update
    implicit none
    real(8),dimension(ndims,ni,nj,nk) :: u0,R,rho
    integer :: m
!
! -- Set variables independant of u,p
    u0 = u
    call geom_update(time0+(t+1)*dt)
    call body_update(t+1,del,ub)
    call domain_update(t+1,u0,del,ub)
    call io_profile('body ')
!
! -- Integrate over dt: int{f,t=0,dt}=(f{0}+f{dt})dt/2
    Huen_loop: do m=1,2
!
! -- Construct the fluid RHS
       call domain_acceleration(R)
       call velo_RHS(t,m,u,u0,R)
!
! -- Adjust RHS for BDIM and inertia
       R = (one-del)*R+del*(ub-u0)*dti
       if(m.eq.1) then
          R = R+dti*u0
       else
          call io_profile('velo ')
          call slip_operator(u,Lu=rho)
          call io_profile('slip ')
          R = R+dti*(u0+u+rho)
       end if
       call io_profile('velo ')
!
! -- Update the pressure and add to RHS: R<=R-grad{p}(1-del)/rho
       call freeint_get(density=rho)
       call io_profile('freei')
       call pressure_update(t,m,one-del,one/rho,R,p)
       call io_profile('press')
!
! -- Update the free interface
       call freeint_update(t+1,m,u,u0)
       call io_profile('freei')
!
! -- Update the velocity: u=R/(I+L)
       R = R*dt
       if(m.eq.2) R = R*0.5
       call slip_inv_operator(R,u)
       call io_profile('slip ')
       call velo_update(t+1,m,u)
       call io_profile('velo ')
    end do Huen_loop
!
! -- Write outputs and log
    call apriori_print(t+1,u,p)
    call io_profile('aprio')
    if(mod(t+1,tmd).eq.0) then
       rewind(res_num)
       write(res_num) (t+1)*dt+time0,u,p
       call flush(res_num)
       call io_profile_report(t+1)
       call io_profile('log  ')
    end if

  end subroutine BDIM_update
!
! -------------------------------------------------------
end program BDIM
