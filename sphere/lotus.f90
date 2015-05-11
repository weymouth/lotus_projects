program sphere_flow
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use gridMod,    only: xg,composite
  use geom_shape  ! sphere geometry
  implicit none

  real,parameter     :: f = 2, D = 64/f, Re = 1000
  integer            :: n(3)
  type(fluid)        :: flow
  type(body)         :: geom
!
! -- Initialize
  n = composite(D*(/6,4,4/))
  call xg(1)%init( n(1), 0.75, 2.25, 1.0, L=D, r=1.037, f=f, d=4.)
  call xg(2)%init( n(2), 0.75, 0.75, 1.0, L=D, r=1.037, f=f)
  call xg(3)%init( n(3), 0.75, 0.75, 1.0, L=D, r=1.037, f=f)

  stop

  geom = sphere(2,1,radius=0.5*D,center=0).map.init_rigid(4,alpha,omega)

  call flow%init(n,geom,V=(/1.,0.,0./),nu=D/Re)
!
! -- Time update loop
  do while (flow%time/D<10)
     call flow%update
!
! -- print
     if(mod(flow%time,D)<flow%dt) then
        call flow%write
        print '("time = ",f0.3)',flow%time/D
     end if
     write(9,'(f10.4,f8.4,3e16.8)') flow%time/D,flow%dt, &
          2.*geom%pforce(flow%pressure)/(3.14159/4.*D**2)
  end do
contains
!
! -- alpha can stay 0 since the body is axisymmetric
  real(8) pure function alpha(t)
    real(8),intent(in) :: t
    alpha = 0
  end function alpha
!
! -- omega*R is 10% of the forward velocity
  real(8) pure function omega(t)
    real(8),intent(in) :: t
    omega = 0.1/(0.5*D)
  end function omega
end program sphere_flow
