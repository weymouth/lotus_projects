program sphere_flow
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use gridMod,    only: xg,composite
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use geom_shape  ! sphere geometry, rigid mapping
  implicit none

  real,parameter     :: f = 1, D = 100/f, Re = 30e3, T = 0.2
  real,parameter     :: dprnt = 0.02
  integer,parameter  :: b(3) = (/4,2,2/)
  logical,parameter  :: gust = .false.
  integer            :: n(3),p
  real               :: t0=0,dt,t1=0,u0,u1,V(3)
  logical            :: root
  type(fluid)        :: flow
  type(body)         :: geom
!
! -- Initialize
  call init_mympi(ndims=3,set_blocks=b)
  root = mympi_rank()==0

  if(gust) then
    geom = sphere(2, 1, radius=0.5*D, center=0)
  else
    geom = sphere(2, 1, radius=0.5*D, center=0).map.init_rigid(1, x, u)
  end if

  n = composite( D*(/5, 2, 2/), prnt=root)
  call xg(1)%stretch(n(1), -10*D, -0.6*D, 1.5*D, 10*D, h_max=5., prnt=root)
  call xg(2)%stretch(n(2), 0., 0., D, 10*D, prnt=root)
  call xg(3)%stretch(n(3), 0., 0., D, 10*D, prnt=root)

  V = 0; V(1) = merge(0.5,1.,gust)
  call flow%init(n/b,geom,V,nu=D/Re)
  if(.not.gust) then
    flow%velocity%e(1)%p = flow%velocity%e(1)%p+0.5
    call flow%velocity%applyBC(flow%mu0)
    flow%u0%e(1)%p = flow%velocity%e(1)%p
  end if
  flow%dt = 0.5

  p = floor(flow%time/D/dprnt)
  if(root) print *,'Init complete'
!
! -- Time update loop
  do while (t1<3.0)
     t0 = flow%time/D
     dt = flow%dt
     t1 = t0+dt/D
     if(t1<2+T+dt/D) then
       if(gust) then
         V(1) = 1-u(real(t1,8))
       else
         call geom%update(t1)
       end if
     end if
     call flow%update(geom=geom,V=V)
     flow%dt = 0.5
!
! -- print
     if(floor(t1/dprnt)>p) then
        call flow%write(geom)
        if(root) print '("time = ",f0.2,", dt = ",f0.3, &
             ", x = ",f0.2,", u = ",f0.2)',t1,dt,x(real(t1,8)),u(real(t1,8))
        p = p+1
     end if
     write(9,'(f10.4,f8.4,3e16.8)') t1,dt, &
          2.*geom%pforce(flow%pressure)/(3.14159/4.*D**2)
  end do
  call mympi_end
contains
!
! -- position and velocity
  real(8) pure function x(ts)
    real(8),intent(in) :: ts
    if(ts<2) then
      x = 0.5*D*ts
    else if(ts<2+T) then
      x = 0.5*D*(ts-0.5*(ts-2)**2/T) !*V(1)
    else
      x = 0.5*D*(2+0.5*T)
    end if
    x = x-D
  end function x

  real(8) pure function u(ts)
    real(8),intent(in) :: ts
    if(ts<2) then
      u = 0.5
    else if(ts<2+T) then
      u = 0.5*(1-(ts-2)/T)
    else
      u = 0
    end if
  end function u
end program sphere_flow
