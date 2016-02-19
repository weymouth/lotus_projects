!-------------------------------------------------------!
!---------------- Static square cylinder ---------------!
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
!
  integer,parameter  :: ndims = 3          ! dimensions
  real,parameter     :: nu = L/Re          ! viscosity
  real,parameter     :: Us = 7             ! reduced velocity
  logical,parameter  :: p(3) = (/.false.,.false.,.true./)  ! periodic BCs
  integer            :: b(3) = (/2,2,4/)   ! blocks
  integer            :: n(3)
  real,parameter     :: s2 = sqrt(2.)/2., D = L*s2
  real               :: area, force(3), t0
  logical            :: root,there=.false.
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
  root = mympi_rank()==0
!
! -- Print run info
  if(root) print *, '-- Square Cylinder --'
!
! -- Initialize array size
  n = composite((/4.5*L,4.5*L,9*L/),prnt=root)
  if(ndims==2) n(3) = 1
!
! -- Initialize and print grid
  call xg(1)%stretch(n(1),-10*L,-0.5*L,L,10*L,h_max=10.,prnt=root)
  call xg(2)%stretch(n(2),-10*L,-1.6*L,1.6*L,10*L,prnt=root)
  if(ndims==3) xg(3)%h = 2
  if(ndims==3.and.root) print *,'length',n(3)*xg(3)%h/D
!
! -- Initialize the square geometry
  square = plane(4,1,(/-s2,s2,0./),(/-0.5*L,0.,0./),0,0) &
       .and.plane(4,1,(/s2,s2,0./),(/0.5*L,0.,0./),0,0) &
       .and.plane(4,1,(/s2,-s2,0./),(/0.5*L,0.,0./),0,0) &
       .and.plane(4,1,(/-s2,-s2,0./),(/-0.5*L,0.,0./),0,0) &
       .map.init_rigid(2,y,v)
  area = L*n(3)*xg(3)%h
!
! -- Initialize fluid
  call flow%init(n/b,square,V=(/1.,0.,0./),nu=nu)
  t0 = flow%time
  if(t0>0) call flow%write(lambda=.true.)
  if(t0==0.and.ndims==3) call flow%velocity%e(3)%perturb(0.05)
  call flow%reset_u0
  force = -square%pforce(flow%pressure)
  if(root) print *, '-- init complete --',t0
!
! -- Time update loop
  do while (flow%time<t0+3*Us*L.and..not.there)
     call square%update(flow%time+flow%dt)
     call flow%update(square)
!
! -- Print fields and force
     force = -square%pforce(flow%pressure)
     write(9,1) flow%time/L,flow%dt,2.*force/area,&
                y(real(flow%time,8))/L,v(real(flow%time,8))
     flush(9)
     if(mod(flow%time,Us*L/24.)<flow%dt) then
       call display(flow%velocity%vorticity_Z(),'01_vort', &
                    box=(/-2,-3,12,6/)*int(L),lim=0.25)
       if(root) print 1,flow%time/L,flow%dt,2.*force/area,&
                  y(real(flow%time,8))/L,v(real(flow%time,8))
     end if
     if(mod(flow%time,Us*L/24.)<flow%dt) call flow%write(lambda=.true.)
     inquire(file='.kill', exist=there)
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
  real(8) pure function y(ts)
    real(8),intent(in) :: ts
    y = 0.5*L*cos(2*pi*ts/L/Us)
  end function y
  real(8) pure function v(ts)  ! rotation velocity
    real(8),intent(in) :: ts
    v = -0.5*sin(2*pi*ts/L/Us)*2*pi/Us
  end function v
end program square_cyl
