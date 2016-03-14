program squeeze
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg
  use geom_shape, only: sphere,operator(.map.),init_scale,pi
  use bodyMod,    only: body
  use fluidMod,   only: fluid
  use gridMod,    only: xg
  use imageMod,   only: display
  implicit none
!
! -- Define parameters, declare variables
  real,parameter    :: L = 32                      ! minor axis size in cells
  real,parameter    :: D = L/2			   ! major axis size in cells
  real,parameter    :: Re_D = 2250                 ! Reynolds number based on D
						   ! The value 2250 comes from a squid travelling at
						   ! 5 L/s with L = 0.03m and nu = 1*10^-6 m^2/s


  real,parameter    :: avgU = 1.		   ! average speed
  real,parameter    :: nu = Re_D/(avgU*D)
  integer,parameter :: n(3) = L*(/5.,1.5,1.5/)     ! number of cells in ijk 
  integer           :: b(3) = (/4,1,1/)            ! MPI domain cuts in ijk
  real,parameter    :: T = L		           ! flow speed and size-change timescale
						   ! in 1*T the squid swims 1 L and oscillates once
  real              :: force(3)                    ! force at t0

  real,parameter    :: lim = 10./D                 ! vorticity level in image
  real,parameter    :: dprnt =0.02                 ! how often to print? every 1 100ieth of period
  logical     :: root,there = .FALSE.              ! flag for stopping
  real(8)     :: ts
  type(fluid) :: flow
  type(body)  :: geom
!
! -- Initialize mpi, grid, body and fluid
  call init_mympi(ndims=3,set_blocks=b)
  root = mympi_rank()==0

  call xg(1)%stretch(n(1),-2.*L,-L,L,5.*L) 	        ! x (wider for motion)
  call xg(2)%stretch(n(2),0.,0.,L,2.*L)		        ! y 
  call xg(3)%stretch(n(3),0.,0.,L,2.*L)                 ! z

  geom = sphere(2, 1, radius=L/2, center=0)&   ! sphere
          .map.init_scale(2,length,rate)&      ! scale in y dimension
	  .map.init_scale(3,length,rate)       ! scale in z dimension

  call flow%init(n/b, geom, V=(/U(0.),0.,0./), nu=nu)
  flow%dt = min(flow%dt,1.)                     ! limit the time step

! -- Time update loop
  do while(flow%time<2.*T .and. .not.there)
    ts = (flow%time+flow%dt)/T              	 		! get non-dimensional time
    call geom%update(real(ts))                  		! apply mapping to geom
    call flow%update(geom,V=(/U(flow%time),0.,0./))     	! update the flow
    flow%dt = min(flow%dt,1.)                   		! limit the time step
    force = -2.*geom%pforce(flow%pressure)/(3.14159*D**2/8.)	! compute the force over half body

! -- write to file
    if(root) write(9,'(f10.4,f8.4,3e16.8)') ts,flow%dt,force
    if(mod(ts,dprnt)<flow%dt/T) then
!    call flow%write(geom)
      if(root) print '(f10.4,5f9.4)',ts,flow%dt, &
            length(ts)/0.455,rate(ts),U(real(ts))
!        call display(flow%velocity%vorticity_Z(),'out_vort', lim=lim)
    end if
    inquire(file='.kill', exist=there)
  end do
  call flow%write()
contains
 real(8) pure function length(ts) 	   ! length scale
    real(8),intent(in) :: ts
    real(8) :: t,c
    t = ts                                 ! pump at once the frequency
    length = (0.925+0.225*cos(2*pi*t))/2.2 ! length scale is divided by squid aspect ratio: 2.2
					   ! and oscillation must go from 1.15 to 0.7 r/r0
					   ! according to literature
  end function length
  real(8) pure function rate(ts)	   ! rate of change of length
    real(8),intent(in) :: ts
    rate = (length(ts+1e-6)-length(ts-1e-6))/(T*2e-6)
  end function rate
  real(4) pure function U(ts)		   ! inflow velocity (or swimming speed)
    real(4),intent(in) :: ts
    real(4) :: t
    t = ts 
!    U = 0.5+(1.-sin(2*pi*t-3.*pi/2.))/2.  ! inflow speed oscillates around 1.5 and 0.5, avgU=1.
    U = 1.
!     U = 0.				   ! test for measuring forces on pulsating body with no flow
  end function U
end program squeeze
