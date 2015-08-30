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
  real,parameter     :: L = 50             ! length scale
  real,parameter     :: Re = 16e2          ! Reynolds number
!
  integer,parameter  :: ndims = 2          ! dimensions
  real,parameter     :: nu = L/Re          ! viscosity
  real,parameter     :: Us = 5.5           ! reduce velocity
  logical,parameter  :: p(3) = (/.false.,.false.,.true./)  ! periodic BCs
  integer            :: b(3) = (/4,4,1/)   ! blocks
  integer            :: n(3)
  real               :: area, force(3), s2 = sqrt(2.)/2.
  real               :: Ufric, yp, t0, k, pos, velo, mass
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
  n = composite((/5*L,5*L,2.5*L/),prnt=root)
  if(ndims==2) n(3) = 1
!
! -- Initialize and print grid
  call xg(1)%stretch(n(1),-10*L,-0.5*L,L,10*L,h_max=10.,prnt=root)
  call xg(2)%stretch(n(2),-10*L,-1.6*L,1.6*L,10*L,prnt=root)
  if(ndims==3) xg(3)%h = 10*L/n(3)
!
! -- Initialize the square geometry
  square = plane(4,1,(/-s2,s2,0./),(/-0.5*L,0.,0./),0,0) &
       .and.plane(4,1,(/s2,s2,0./),(/0.5*L,0.,0./),0,0) &
       .and.plane(4,1,(/s2,-s2,0./),(/0.5*L,0.,0./),0,0) &
       .and.plane(4,1,(/-s2,-s2,0./),(/-0.5*L,0.,0./),0,0) &
       .map.init_rigid(2,y,v)
  area = L*n(3)*xg(3)%h
  mass = 2.64*(0.5*L**2)*n(3)*xg(3)%h
  k = (1.38*mass)*(2*3.14159/Us/L)**2
!
! -- Initialize fluid
  call flow%init(n/b,square,V=(/1.,0.,0./),nu=nu)
  t0 = flow%time
  force = -square%pforce(flow%pressure)
  if(root) print *, '-- init complete --',t0
!
! -- Time update loop
  do while (flow%time<t0+250*L)
     call motion_update()
     call square%update(flow%time+flow%dt)
     call flow%update(square)
!
! -- Print fields and force
     force = -square%pforce(flow%pressure)
     write(9,1) flow%time/L,flow%dt,2.*force/area,pos/L,velo
     flush(9)
     if(mod(flow%time,0.5*L)<flow%dt) then
       img = flow%velocity%vortZrender(15./L)
       call img%resample((/-300,-320,1280,640/),smooth=.true.)
       call img%write('vortZ',flow%time)
     end if
     if(mod(flow%time,5*L)<flow%dt) then
        call flow%write(average=.true.)
        if(root) print 1,flow%time/L,flow%dt,2.*force/area,pos/L,velo
     end if
1    format(f10.4,f8.4,5e14.6)
  end do
  if(root) print *, '--- complete ---'
!
! -- Finalize MPI
#if MPION
  call mympi_end
#endif
contains
!
! -- motion definitions
  subroutine motion_update()
    logical :: first=.true.
    real :: dt, dt0=0, am1 = 0, a0, v0 = 0, y0 = 0
    dt = flow%dt
    a0 = (force(2)-y0*k)/mass
    if(first) then
      dt0 = dt; am1 = a0; first = .false.
    end if
    velo = v0+dt*(a0+0.5*dt*(a0-am1)/dt0)
    pos = y0+dt*(v0+0.5*dt*a0)
    dt0 = dt; am1 = a0
    v0 = velo; y0 = pos
  end subroutine motion_update
  real(8) pure function y(ts)
    real(8),intent(in) :: ts
    y = pos
  end function y
  real(8) pure function v(ts)  ! rotation velocity
    real(8),intent(in) :: ts
    v = velo
  end function v
end program square_cyl
