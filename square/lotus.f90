!-------------------------------------------------------!
!-------------------- square cylinder ------------------!
!-------------------------------------------------------!
program square_cyl
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  use imageMod,   only: display
  use geom_shape  ! to define geom (set,eps,plane, etc)
  implicit none
  real,parameter     :: L = 50             ! length scale
  real,parameter     :: Re = 5e3           ! Reynolds number
  real,parameter     :: Omega=0.02, Per=2  ! motion frequency, # periods
!
  integer,parameter  :: ndims = 2          ! dimensions
  real,parameter     :: nu = L/Re          ! viscosity
  integer            :: b(3) = (/2,2,1/)   ! blocks
  real               :: m = 4.5*L          ! approximate # points
  integer            :: n(3)               ! # points
  real               :: force(3),moment(3) ! pressure integrals
  logical            :: root,there=.false. ! logical flags
!
  type(fluid)        :: flow
  type(body)         :: square
!
! -- Initialize MPI (if MPI is ON)
#if MPION
  call init_mympi(ndims,set_blocks=b(:ndims))
#else
  b=1
#endif
  root = mympi_rank()==0
!
! -- Initialize array size
  n = composite((/m,m,1./),prnt=root)
!
! -- Initialize and print grid
  call xg(1)%stretch(n(1),-10*L,-L,L,10*L,prnt=root)
  call xg(2)%stretch(n(2),-10*L,-L,L,10*L,prnt=root)
!
! -- Initialize the square geometry
  square = make_bildge_geom().map.init_rigid(6,phi)
!
! -- Initialize fluid
  call flow%init(n/b,square,nu=nu)
  if(root) print *, '-- init complete --'
!
! -- Time update loop
  do while (flow%time*Omega/(2*pi)<Per.and..not.there)
     call square%update(flow%time+flow%dt)
     call flow%update(square)
!
! -- Print fields and force
     force = -square%pforce(flow%pressure)
     moment = -square%pmoment(flow%pressure)
     write(9,1) flow%time*Omega/(2*pi),flow%dt,phi(real(flow%time,8)),&
                    2.*force(:2)/L,2.*moment(3)/L**2
     flush(9)
     if(mod(flow%time,0.25*pi/Omega)<flow%dt) then
       call display(flow%velocity%vorticity_Z(),'vort',lim=0.25)
       print 1,flow%time*Omega/(2*pi),flow%dt,phi(real(flow%time,8)),&
                      2.*force(:2)/L,2.*moment(3)/L**2
     end if
     inquire(file='.kill', exist=there)
1    format(f10.4,f8.4,4e14.6)
  end do
  call flow%write(square)
  if(root) print *, '--- complete ---'
!
! -- Finalize MPI
#if MPION
  call mympi_end
#endif
contains
!
! -- motion definition
  real(8) pure function phi(t)
    real(8),intent(in) :: t
    phi = pi*sin(Omega*t)
  end function phi
!
! -- shape definition
  type(set) function make_bildge_geom() result(geom)
    real,parameter :: s2 = sqrt(2.)/2.
    geom = plane((/-s2,s2,0./),(/-0.5*L,0.,0./)) &
         .and.plane((/s2,s2,0./),(/0.5*L,0.,0./)) &
         .and.plane((/s2,-s2,0./),(/0.5*L,0.,0./)) &
         .and.plane((/-s2,-s2,0./),(/-0.5*L,0.,0./))
  end function
end program square_cyl
