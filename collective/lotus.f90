!
! -- use must check field.distribute before running.
!    This code requires flux be distributed over
!         if(d==3 .and. mympi_wall(d))
!
program squeeze
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  use geom_shape
  use bodyMod,    only: body
  use fluidMod,   only: fluid
  use gridMod,    only: xg
  use imageMod,   only: display
  implicit none
!
! -- Define parameters, declare variables
  real,parameter    :: L = 64         ! major axis size in cells
  real,parameter    :: D = L/2.2  	   ! minor axis size in cells
  real,parameter    :: H = D/4.        ! head radius
  real,parameter    :: R = L/12.       ! jet radius
  real,parameter    :: Re_D = 2500     ! Reynolds number based on D
						   ! The value 2250 comes from a squid travelling at
						   ! 5 L/s with L = 0.03m and nu = 1*10^-6 m^2/s
  real,parameter    :: T = 2*L	       ! size-change timescale
               ! in 1*T the squid swims 1 L and oscillates once
  real,parameter    :: g = 1.0         ! center-to-center gap
  real,parameter    :: phi = 0.0       ! phase difference

  integer           :: b(3) = (/1,1,1/)            ! MPI domain cuts in ijk
  real,parameter    :: lim = 10./D                 ! vorticity level in image
  real,parameter    :: dprnt = 0.04                ! how often to print?
  real,parameter    :: Tend  = 7*T                 ! when should we stop?

  real,parameter    :: nu = D/Re_D                 ! kinematic viscosity
  real,parameter    :: s  = H/(L/2.)               ! head scaling factor
  real,parameter    :: A0 = pi*(D/2.)**2/4.        ! mean frontal area
  real,parameter    :: V0 = 2./3.*pi*L*(D/2.)**2   ! mean frontal area
  real,parameter    :: Aj = pi*R**2                ! jet area
  real,parameter    :: Ai = pi*((H+2*R)**2-H**2)   ! inlet area
  integer           :: n(3)                        ! number of cells in ijk
  real              :: force(3)                    ! force
  real              :: rate_over_length

  logical           :: root,there = .FALSE.        ! flag for stopping
  real(8)           :: t1=0.,dt
  type(fluid)       :: flow
  type(body)        :: geom
!
! -- Initialize mpi, grid, body and fluid
  call init_mympi(ndims=3,set_blocks=b)
  root = mympi_rank()==0

  n = composite(L*(/2.,g,1./),prnt=root)                     ! set n
  call xg(1)%stretch(n(1),-2.*L,-0.6*L,L,6.*L, &
                     h_min=2.,h_max=6.,prnt=root)            ! x
  call xg(3)%stretch(n(3),0.,0.,0.4*L,4*L,prnt=root)         ! z

  if(root) print '("n = ",i4,", center-to-center gap = ",f6.3)',n(2),g

  geom = (sphere(2, 1, radius=L/2., center=0) &   ! sphere
          .map.init_scale(2,length,rate)   &      ! scale in y dimension
       	  .map.init_scale(3,length,rate))  &      ! scale in z dimension
     .or.(sphere(2, 0, radius=L/2., center=0) &   ! dont measure force
          .map.init_scale(2,length2,rate2) &      ! scale in y dimension
          .map.init_scale(3,length2,rate2) &      ! scale in z dimension
          .map.init_rigid(2,fixed,zero))          ! shift

  ! geom = (sphere(2, 1, radius=L/2., center=0) &   ! sphere
  !         .map.init_scale(2,length,rate) &      ! scale in y dimension
  !      	  .map.init_scale(3,length,rate) &      ! scale in z dimension
  !         .map.velo_fld(jet)) &
  !    .or.(sphere(2, 1, radius=L/2., center=(/L/2.,0.,0./)) &
  !         .map.init_scale(2,fixed,zero)  &      ! scale in y dimension
  !         .map.init_scale(3,fixed,zero))        ! scale in z dimension

  call flow%init(n/b, geom, V=(/U(t1),0.,0./), nu=nu)
  flow%velocity%e(2)%flux_correct = .false.

! -- Time update loop
!  do while(flow%time<Tend .and. .not.there)
    dt = flow%dt                                	 		! time step
    t1 = (flow%time+dt)/T                        	 		! time at end of step
    rate_over_length = V0*rate(t1)/length(t1)
    call geom%update(real(t1))                    		! apply mapping to geom
    call flow%update(geom,V=(/U(t1),0.,0./))          ! update the flow
    force = -2.*geom%pforce(flow%pressure)/A0         ! compute the force coefficient

! -- write to file
    if(root) write(9,'(f10.4,f8.4,3e16.8)') t1,dt,force
    if(root) flush(9)
    call display(flow%velocity%vorticity_Z(),name='01_out',lim=lim)
    if(mod(t1,dprnt)<dt/T) then
      if(root) print '(f10.4,7f9.4)',t1,dt,length(t1),rate(t1),&
                      U(t1),jet(real((/L/2.,0.,-H-R/),8))
      if(root) flush(9)
      call display(flow%velocity%vorticity_Z(),name='01_out',lim=lim)
    end if
    call flow%write(geom)
    inquire(file='.kill', exist=there)
!    if(mod(t1,1.)<dt/T) call flow%write(geom,lambda=.TRUE.)
!  end do
!  call flow%write()
  call mympi_end
contains
 real(8) pure function length(t1) 	   ! length scale
    real(8),intent(in) :: t1
    length = (0.925+0.225*cos(2*pi*t1))*D/L ! length scale is divided by squid aspect ratio: 2.2
					   ! and oscillation must go from 1.15 to 0.7 r/r0
					   ! according to literature
  end function length
  real(8) pure function rate(t1)	   ! rate of change of length
    real(8),intent(in) :: t1
    rate = (length(t1+1e-6)-length(t1-1e-6))/(T*2e-6)
  end function rate
  real(8) pure function length2(t1) 	   ! length scale
     real(8),intent(in) :: t1
     length2 = length(t1+phi)
   end function length2
   real(8) pure function rate2(t1)	   ! rate of change of length
     real(8),intent(in) :: t1
     rate2 = rate(t1+phi)
   end function rate2
  pure function jet(x) result(v)
    real(8),intent(in) :: x(3)
    real(8)            :: v(3)
    real :: x2
    v = 0
    if(mod(t1,1.)<0.5) then
      x2 = (x(3)+H+R)**2+x(2)**2
      if(x2<R**2 .and. x(1)>L/3.) v(1) = -rate_over_length/Aj
    else
      x2 = x(3)**2+x(2)**2
      if(x2<(H+2*R)**2 .and. x(1)>L/3.) v(1) = -rate_over_length/Ai
    end if
  end function jet

  real(8) pure function fixed(t1)	   ! length scale
     real(8),intent(in) :: t1
     fixed = g*L
   end function fixed
   real(8) pure function zero(t1)	   ! rate of change of length
     real(8),intent(in) :: t1
     zero = 0.
   end function zero

  real function U(t1)		        ! inflow velocity (or swimming speed)
    real(8),intent(in) :: t1
    real :: ts
    U = 1
    return
    ts = mod(t1,1.)
    if(ts<0.5) then
      U = 0.5+(ts/4.-sin(4*pi*ts)/8./pi)*8.
    else
      U = 0.5+(1./4.-ts/4.)*8.
    end if
  end function U
end program squeeze
