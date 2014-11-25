!-------------------------------------------------------!
!---------------- Static square cylinder ---------------!
!-------------------------------------------------------!
program square_cyl
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  use geom_shape  ! to define geom (set,eps,plane, etc)
  implicit none
  real,parameter     :: f = 1.0            ! scaling factor
  real,parameter     :: L = 200/f          ! length scale
  real,parameter     :: Re = 400           ! Reynolds number
!
  integer,parameter  :: ndims = 2          ! dimensions
  real,parameter     :: nu = L/Re          ! viscosity
  logical,parameter  :: p(3) = (/.false.,.false.,.true./)  ! periodic BCs
  integer,parameter  :: b(3) = (/4,4,1/)   ! blocks
  integer            :: n(3)
  real               :: area,force(3)
  real               :: Ufric,yp,t0
!
  type(fluid)        :: flow
  type(body)         :: square
!
! -- Initialize MPI (if MPI is ON)
#if MPION
  call init_mympi(ndims,set_blocks=b(:ndims),set_periodic=p(:ndims))
#else
  b=1
#endif
!
! -- Print run info
  Ufric = sqrt(0.026/Re**(1./7.)/2.)
  yp = Ufric/nu
  if(mympi_rank()==0) print *, '-- Square Cylinder --'
  if(mympi_rank()==0) print '("   L=",f0.4,", nu=",f0.4,", y+=",f0.4)',L,nu,yp
!
! -- Initialize array size
  n = composite((/169+L+20,121+L,1./)/b)
  if(ndims==2) n(3) = 1
!
! -- Initialize and print grid
  call xg(1)%init(n(1)*b(1),.5*L,.5*L,0.69,c=L*.15,f=f,r=1.12)
  call xg(2)%init(n(2)*b(2),.5*L,.5*L,1.00,c=L*.25,f=f,r=1.12)
  if(ndims==3) xg(3)%h = 6*b(3)*L/n(3)
  if(mympi_rank()==0) then
     call xg(1)%write
     call xg(2)%write
     call xg(3)%write
     print '("   total points=",i0)', product(n*b)
  end if
!!$  call mympi_end()
!!$  stop
!
! -- Initialize the square geometry
  square = plane(4,1,(/-1,0,0/),(/-L/2.,0.,0./),0,0) &
       .and.plane(4,1,(/ 1,0,0/),(/+L/2.,0.,0./),0,0) &
       .and.plane(4,1,(/0,-1,0/),(/0.,-L/2.,0./),0,0) &
       .and.plane(4,1,(/0, 1,0/),(/0.,+L/2.,0./),0,0)
  area = L*n(3)*b(3)*xg(3)%h
!
! -- Initialize fluid
  call flow%init(n,square,V=(/1.,0.,0./),nu=nu)
  t0 = flow%time
  if(ndims==3) call flow%velocity%e(3)%perturb(0.02)
  if(mympi_rank()==0) print *, '-- init complete --'
!
! -- Time update loop
  do while (flow%time<t0+500*L)
     call flow%update
!
! -- Print fields and force
     force = square%pforce(flow%pressure)
     write(9,'(f10.4,f8.4,3e16.8)') flow%time/L,flow%dt,2.*force/area
     flush(9)
     if(mod(flow%time,5*L)<flow%dt) then
        call flow%write
        if(mympi_rank()==0) print '(f10.4,f8.4,3e16.8)',flow%time/L,flow%dt,2.*force/area
     end if
  end do
  if(mympi_rank()==0) write(6,*) '--- complete --- '
!
! -- Finalize MPI
#if MPION
  call mympi_end
#endif
end program square_cyl
