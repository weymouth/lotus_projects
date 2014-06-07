!-------------------------------------------------------!
!--------------- Foil Impulse Test case ----------------!
!-------------------------------------------------------!
program foil_impulse
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use geom_shape  ! to create geometry
  implicit none
  integer,parameter  :: f=2**6             ! resolution  
  real,parameter     :: Re = 1000          ! Reynolds number
  real(8),parameter  :: alpha = 10         ! AOA
  integer,parameter  :: b(3) = (/2,4,2/)   ! blocks
  integer,parameter  :: d(3) = (/6,4,6/)   ! domain size
  real,parameter     :: T = 0.5            ! yank period
  real,parameter     :: U = -3             ! yank speed
!
  integer,parameter  :: ndims = 2   ! dimensions
  real(8),parameter  :: L = f       ! length
  integer,parameter  :: m(3) = f*d  ! points
  real,parameter     :: nu = L/Re   ! viscosity
  real(8),parameter  :: yc = m(2)/2 ! location
  integer            :: n(3)
  real               :: force(3),area
!
  type(fluid)        :: flow
  type(body)         :: foil
  type(set)          :: geom
  type(model_info)   :: info
!
! -- Initialize MPI (if MPI is ON)
#if MPION
  call init_mympi(ndims,set_blocks=b(:ndims))
!
! -- Get grid size
  n = m/b
#else
  n = m
#endif
  if(ndims==2) n(3) = 1
  if(mympi_rank()==0) print *, '-- Foil Impulse --'
  if(mympi_rank()==0) print '("   L=",i0," nu=",f0.4)', f,nu
!
! -- Initialize the foil geometry
!  surface_debug = .true.
  model_fill = .false.
  info%file = 'naca_square.IGS'
  info%x = (/-4.219,-10.271,0./)
  info%s = 0.36626*L
  geom = (model_init(info) &
       .map.((init_affn()**(/alpha,0.D0,0.D0/))+(/yc,yc,0.D0/)))
  if(ndims==3) geom = geom.and.plane(4,1,(/0,0,-1/),yc,0,0)
  foil = geom
!!$  foil =    plane(4,1,(/-1,0,0/),(/yc-0.5*L,yc,yc/),0,0) &
!!$       .and.plane(4,1,(/ 1,0,0/),(/yc+0.5*L,yc,yc/),0,0) &
!!$       .and.plane(4,1,(/0,-1,0/),yc,0,0)
  area = L
  if(ndims==3) area = L*(m(3)-yc)
!
! -- Initialize fluid
!!$  call flow%init(n,foil,V=(/1.,0.,0./),nu=nu,NN=10)
  call flow%init(n,foil,V=(/0.,0.,0./),nu=nu)
  call flow%write
  if(mympi_rank()==0) print *, '-- init complete --'
!
! -- Run it
  do while (flow%time<0*L)
!
! -- parameters
     
!!$     flow%g(1) = 0.5*sin(flow%time/L)/L
!!$     flow%velocity%e(1)%bound_val = 0.5*(1.-cos(flow%time/L))
!!$     if(flow%time>pi*L) then
!!$        flow%g(1) = 0; flow%velocity%e(1)%bound_val = 1.0
!!$     end if
!!$     if(mympi_rank()==0) print '("   t=",f0.4," g=",f0.4," u=",f0.4)', &
!!$          flow%time/L,flow%g(1),flow%velocity%e(1)%bound_val

!!$     flow%g(ndims) = U*sin(flow%time/L/T*2.*pi)/L/T*2.*pi
!!$     flow%velocity%e(ndims)%bound_val = U*(1.-cos(flow%time/L/T*2.*pi))
!!$     if(flow%time>T*L) then
!!$        flow%g(ndims) = 0; flow%velocity%e(ndims)%bound_val = 0
!!$     end if
!!$     if(mympi_rank()==0) print '("   t=",f0.4," g=",f0.4," u=",f0.4)', &
!!$          flow%time/L,flow%g(ndims),flow%velocity%e(ndims)%bound_val
!
!-- run it
     call flow%update
     if(mod(flow%time,L)<flow%dt) call flow%write
!
! -- Print force on the square
     force = foil%pforce(flow%pressure)
     write(9,'(f10.4,f7.4,3e16.8)') flow%time/L,flow%dt,2.*force/area
     flush(9)
  end do
  if(mympi_rank()==0) write(6,*) '--- complete --- '
!
! -- Finalize MPI
#if MPION
  call mympi_end
#endif
end program foil_impulse
!
