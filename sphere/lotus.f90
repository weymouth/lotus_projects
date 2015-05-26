program sphere_flow
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use gridMod,    only: xg,composite
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use geom_shape  ! sphere geometry, rigid mapping
  implicit none

  real,parameter     :: f = 1, D = 64/f, Re = 2100
  integer,parameter  :: b(3) = (/4,4,1/)
  integer            :: n(3)
  real               :: t0=0,dt,t1,u0,u1
  type(fluid)        :: flow
  type(body)         :: geom
!
! -- Initialize
  call init_mympi(ndims=3,set_blocks=b)

  geom = sphere(2,1,radius=0.5*D,center=0).map.init_rigid(3,z,w).map.init_rigid(6,phi,r)

  n = composite( D*(/4,4,6/), prnt=mympi_rank()==0 )
  call xg(1)%init( n(1), 1.0, 1.0,  1., L=D, f=f, prnt=mympi_rank()==0 )
  call xg(2)%init( n(2), 1.0, 1.0,  1., L=D, f=f, prnt=mympi_rank()==0 )
  call xg(3)%init( n(3), 4.5, 0.6, 15., L=D, f=f, prnt=mympi_rank()==0 )

  call flow%init(n/b,geom,nu=D/Re)
!
! -- Time update loop
  do while (t0<4.5)
     t0 = flow%time/D
     dt = flow%dt
     t1 = t0+dt/D
     call geom%update(t1)
     call flow%update(geom)
!
! -- print
     if(mod(t1,0.1)<dt/D) then
        call flow%write(geom)
        if(mympi_rank()==0) print '("time = ",f0.2,", dt = ",f0.3, &
             ", z = ",f0.2,", w = ",f0.2)',t1,dt,z(real(t1,8)),w(real(t1,8))
     end if
     write(9,'(f10.4,f8.4,3e16.8)') t1,dt, &
          2.*geom%pforce(flow%pressure)/(3.14159/4.*D**2)
  end do
  call mympi_end
contains
!
! -- position and velocity
  real(8) pure function z(ts)
    real(8),intent(in) :: ts
    z = -D*log(cosh(ts))
  end function z

  real(8) pure function w(ts)
    real(8),intent(in) :: ts
    w = -tanh(ts)
  end function w
!
! -- angle and spin
  real(8) pure function phi(ts)
    real(8),intent(in) :: ts
    phi = 0 ! not important
  end function phi

  real(8) pure function r(ts)
    real(8),intent(in) :: ts
    r = -2.4/D
  end function r
end program sphere_flow
