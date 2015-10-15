!-------------------------------------------------------!
!---------------- Turbulent Channel Flow ---------------!
!-------------------------------------------------------!
program channel
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use gridMod,    only: xg
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use imageMod,   only: display
  use geom_shape
  implicit none
  integer,parameter :: f = 128, s=3, w=2
  real,parameter    :: L = w*f-2*s, nu = L/5600
  integer,parameter :: n(3) = (/128,f,64/)
  integer     :: b(3) = (/4,4,1/)
  real        :: z(3) = (/0.,0.,3.14159*L/2./)
  logical     :: root,there=.false.
  real        :: Ub,Utau
  type(fluid) :: flow
  type(body)  :: wall
!
! -- Initialize MPI
  call init_mympi(3,set_blocks=b,set_periodic=(/.true.,.false.,.true./))
  root = mympi_rank()==0
!
! -- Init grid
  xg%h = L*3.14159/n(3)
  call xg(2)%squeeze(n(2),0.,real(w*f),prnt=root)
!
! -- Init channel walls
  wall = plane(4,1,(/0,1,0/),(/0.,real(s),0./),0,0).or. &
!         sphere(2,1,0.1*L,0.2*L,0.,0.).or.sphere(2,1,0.1*L,z+0.2*L,0.,0.).or.&
!         sphere(2,1,0.1*L,0.4*L,0.,0.).or.sphere(2,1,0.1*L,z+0.4*L,0.,0.).or.&
!         sphere(2,1,0.1*L,0.6*L,0.,0.).or.sphere(2,1,0.1*L,z+0.6*L,0.,0.).or.&
!         sphere(2,1,0.1*L,0.8*L,0.,0.).or.sphere(2,1,0.1*L,z+0.8*L,0.,0.).or.&
         plane(4,1,(/0,-1,0/),(/0.,L+real(s),0./),0,0)
!
! -- Initialize fluid
  call flow%init(n/b,wall,nu=nu)
  if(flow%time==0) call flow%u0%e(1)%eval(parabola)
  if(root) print *, '-- init complete --'
!
! -- Time update loop
  do while (.not.there)
     call flow%update(b=nu*(/12.,0.,0./)/L**2)
     Ub = flow%velocity%e(1)%mean()
     Utau = fric(flow%velocity%e(1))
!
! -- Print fields and force
     write(9,'(f10.4,f8.4,3e14.6)') flow%time*Ub/L,flow%dt*Ub,Ub*L/nu,Utau*L/2./nu
     flush(9)
     if(mod(flow%time,0.1*L/Ub)<flow%dt) &
        call display(flow%velocity%vorticity_Z(),'01vort',lim=14.*Ub/L)
     inquire(file='.kill', exist=there)
  end do
  call flow%write()
  if(root) print *, '--- complete ---'
  call mympi_end()
contains
  real pure function parabola(x) result(u)
    real,intent(in) :: x(3)
    real :: y
    y = (x(2)-s)/L
    u = 6.*y*(1.-y)
  end function parabola
  real function fric(u)
    use gridMod,  only: dv
    use fieldMod, only: field
    use mympiMod, only: mympi_sum,mympi_wall
    type(field),intent(in) :: u
    real,pointer :: up(:,:,:)
    real :: area=0
    if(area==0) then
      area = sum(real(dv(:,s+1,:),8))
      area = merge(area,0.,mympi_wall(-2))
      call mympi_sum(area)
    end if
    up => u%point()
    fric = sum(real(up(:,s+1,:)*dv(:,s+1,:),8))/area
    fric = merge(fric,0.,mympi_wall(-2))
    call mympi_sum(fric)
    fric = sqrt(nu*fric/0.5)
  end function fric
end program channel
