!-------------------------------------------------------!
!---------------- Static square cylinder ---------------!
!-------------------------------------------------------!
program square_cyl
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  use geom_shape  ! to define geom (set,eps,plane, etc)
  use ioMod
  use imageMod
  use fieldMod
  use vectorMod
  implicit none
  real,parameter     :: L = 32             ! length scale
  real,parameter     :: Re = 16e2          ! Reynolds number
!
  integer,parameter  :: ndims = 2          ! dimensions
  real,parameter     :: nu = L/Re          ! viscosity
  logical,parameter  :: p(3) = (/.false.,.false.,.true./)  ! periodic BCs
  integer            :: b(3) = (/4,4,1/)   ! blocks
  integer            :: n(3)
  real               :: area, force(3), s2 = sqrt(2.)/2.
  real               :: Ufric, yp, t0
  logical            :: root
!
  type(fluid)        :: flow
  type(body)         :: square
  type(rgbimage)     :: img
!
! -- Initialize MPI (if MPI is ON)
#if MPION
  call init_mympi(ndims,set_blocks=b(:ndims),set_periodic=p(:ndims))
#else
  b=1
#endif
  root = mympi_rank()==0
!
! -- Print run info
  Ufric = sqrt(0.026/Re**(1./7.)/2.)
  yp = Ufric/nu
  if(root) print *, '-- Square Cylinder --'
  if(root) print '("   L=",f0.4,", nu=",f0.4,", y+=",f0.4)',L,nu,yp
!
! -- Initialize array size
!  n = composite((/4*L,2.5*L,2.5*L/),prnt=root)
  n = composite((/6*L,4*L,2.5*L/),prnt=root)
  if(ndims==2) n(3) = 1
!
! -- Initialize and print grid
  call xg(1)%stretch(n(1),-10*L,-s2*L,s2*L,10*L,h_max=10.,prnt=root)
  call xg(2)%stretch(n(2),-10*L,-s2*L,s2*L,10*L,prnt=root)
  if(ndims==3) xg(3)%h = 10*L/n(3)
!
! -- Initialize the square geometry
  square = plane(4,1,(/-s2,s2,0./),(/-s2*L,0.,0./),0,0) &
       .and.plane(4,1,(/s2,s2,0./),(/s2*L,0.,0./),0,0) &
       .and.plane(4,1,(/s2,-s2,0./),(/s2*L,0.,0./),0,0) &
       .and.plane(4,1,(/-s2,-s2,0./),(/-s2*L,0.,0./),0,0)
  area = L*n(3)*xg(3)%h
!
! -- Initialize fluid
  call flow%init(n/b,square,V=(/1.,0.,0./),nu=nu)
  t0 = flow%time
  if(root) print *, '-- init complete --',t0
!
! -- Time update loop
  do while (flow%time<t0+25*L)
     call flow%update
!
! -- Print fields and force
     force = square%pforce(flow%pressure)
     write(9,'(f10.4,f8.4,3e16.8)') flow%time/L,flow%dt,2.*force/area
     flush(9)
     if(mod(flow%time,0.05*L)<flow%dt) then
       img = flow%velocity%vortZrender(15./L)
       call img%resample((/-300,-320,1280,640/),smooth=.true.)
       call img%write('vortZ',flow%time)
!       call mympi_end
!       stop
     end if
     if(mod(flow%time,5*L)<flow%dt) then
        call flow%write(average=.true.)
        if(root) print '(f10.4,f8.4,3e16.8)',flow%time/L,flow%dt,2.*force/area
     end if
  end do
  if(root) print *, '--- complete ---'
!
! -- Finalize MPI
#if MPION
  call mympi_end
#endif
end program square_cyl
