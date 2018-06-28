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
  real,parameter    :: L = 128          ! major semi-axis size in cells
  real,parameter    :: beta0 = 0.25	   ! aspect ratio
  real,parameter    :: beta1 = 0.0375	 ! pulse amplitude
  real,parameter    :: A0_Ae = 4    	 ! area ratio
  real,parameter    :: per = 3         ! periods of motion
  real,parameter    :: thk = 2         ! membrane half thickness in cells
  real,parameter    :: Re = 25000      ! Approx reynolds number
						   ! Uj=10x(2L)/s with 2L = 0.05m and nu = 1*10^-6 m^2/s

  real,parameter    :: dm = 4*beta1/beta0       ! fraction of mass expelled
  real,parameter    :: xe = sqrt(1-1/A0_Ae)     ! exit location
  real,parameter    :: Vf = 2./3.+xe-xe**3/3.   ! volume factor V= Vf pi L**3
  real,parameter    :: T = Vf*dm/2.*A0_Ae*L     ! size-change timescale
  real,parameter    :: A0 = pi*(beta0*L)**2     ! mean frontal area
  real,parameter    :: nu = 2*L/Re              ! kinematic viscosity

  real,parameter    :: dprnt = pi*T/12.            ! how often to print?
  real,parameter    :: Tend  = per*2*pi*T          ! when should we stop?
  real,parameter    :: f(3)  = (/3.,0.75,0.75/)    ! approx grid factor
  integer           :: b(3)  = (/8,2,1/)           ! MPI domain cuts in ijk
  integer           :: n(3)                        ! number of cells in ijk
  real              :: force(3),dt                 ! force
  logical           :: root,there = .FALSE.        ! flag for stopping
  type(fluid)       :: flow
  type(body)        :: geom
!
! -- Initialize mpi, grid, body and fluid
  call init_mympi(ndims=3,set_blocks=b)
  root = mympi_rank()==0

  if(root) print *,'T',T,'xe',xe,'dm',dm,'Vf',Vf,'nu',nu

  n = composite(L*f,prnt=root)                             ! n
  call xg(1)%stretch(n(1),-4.*L,-1.2*L,2*L,6.*L, &
                     h_min=2.,h_max=6.,prnt=root)          ! x
  call xg(2)%stretch(n(2),0.,0.,0.4*L,L,prnt=root)         ! y
  call xg(3)%stretch(n(3),0.,0.,0.4*L,L,prnt=root)         ! z

  geom = ((sphere(radius=L+thk, center=0) &
            .and.plane(norm=(/1,0,0/), center=(/L*xe,0.,0./))) &
            .map.init_scale(2,betao,dbeta) &
            .map.init_scale(3,betao,dbeta)) &
        -(sphere(radius=L-thk, center=0) &
            .map.init_scale(2,betai,dbeta) &
            .map.init_scale(3,betai,dbeta))
  geom%dis_wall(2:3) = .true. ! ok to adjust velocity on +y,+z planes

  call flow%init(n/b, geom, V=(/0.,0.,0./), nu=nu)
  call flow%write(geom)

! -- Time update loop
  do while(flow%time<Tend .and. .not.there)
    flow%dt = min(flow%dt,1.)
    dt = flow%dt
    call geom%update(flow%time+dt)            		! apply mapping to geom
    call flow%update(geom)                        ! update the flow
    force = -2.*geom%pforce(flow%pressure)/A0     ! compute the force coefficient

! -- write to file
    if(root) write(9,'(f10.4,f8.4,3e16.8)') flow%time/T,flow%dt,force
    if(root) flush(9)
    if(mod(flow%time,dprnt)<dt) then
      if(root) print '(f10.4,4f8.4)',flow%time/T,flow%dt, &
            betao(real(flow%time,8)),betai(real(flow%time,8)), &
            dbeta(real(flow%time,8))*T
      call display(flow%velocity%vorticity_Z(),name='01_out')
      call flow%write(geom,lambda=.TRUE.)
    end if
    inquire(file='.kill', exist=there)
  end do
  call mympi_end
contains
 real(8) pure function betao(t1)  ! length scale of membrane outside
    real(8),intent(in) :: t1
    betao=(L*beta0+thk)/(L+thk)+beta1*cos(t1/T)
  end function betao
  real(8) pure function betai(t1)  ! length scale of membrane inside
     real(8),intent(in) :: t1
     betai=(L*beta0-thk)/(L-thk)+beta1*cos(t1/T)
   end function betai
  real(8) pure function dbeta(t1)	   ! rate of change of length
    real(8),intent(in) :: t1
    dbeta=-beta1*sin(t1/T)/T
  end function dbeta

end program squeeze
