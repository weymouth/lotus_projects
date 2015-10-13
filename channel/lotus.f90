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
  real,parameter    :: L = 128-6, nu = 9e-3
  integer,parameter :: n(3) = (/128,128,64/)
  integer     :: b(3) = (/2,1,1/)
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
  xg(2)%h = 1.
!
! -- Init channel walls
  wall = plane(4,1,(/0,1,0/),(/0.,3.,0./),0,0).or. &
!         sphere(2,1,5.,0.2*L,0.,0.).or.sphere(2,1,5.,z+0.2*L,0.,0.).or.&
!         sphere(2,1,5.,0.4*L,0.,0.).or.sphere(2,1,5.,z+0.4*L,0.,0.).or.&
!         sphere(2,1,5.,0.6*L,0.,0.).or.sphere(2,1,5.,z+0.6*L,0.,0.).or.&
!         sphere(2,1,5.,0.8*L,0.,0.).or.sphere(2,1,5.,z+0.8*L,0.,0.).or.&
         plane(4,1,(/0,-1,0/),(/0.,L+3.,0./),0,0)
!
! -- Initialize fluid
  call flow%init(n/b,wall,nu=nu)
!  call IC(flow%u0%e(1))
  if(root) print *, '-- init complete --'
!
! -- Time update loop
  do while (.not.there)
     call flow%update(b=nu*(/8.,0.,0./)/L**2)
     Ub = flow%velocity%e(1)%mean()
     Utau = fric(flow%velocity%e(1))
!
! -- Print fields and force
     write(9,'(f10.4,f8.4,3e14.6)') flow%time*Ub/L,flow%dt*Ub,Ub*L/nu,Utau*L/2./nu
     flush(9)
     if(mod(flow%time,L/Ub)<flow%dt) &
        call display(flow%velocity%vorticity_Z(),'vort',lim=14.*Ub/L)
     inquire(file='.kill', exist=there)
  end do
  call flow%write()
  if(root) print *, '--- complete ---'
  call mympi_end()
contains
  subroutine IC(u)
    use fieldMod, only: field
    type(field),intent(inout) :: u
    real :: x(3),y
    integer :: j
    do concurrent(j=u%js:u%je)
        x = u%pos(1,j,1); y = (x(2)-2)/L
        u%p(:,j,:) = 4.*y*(1.-y)
    end do
    call u%applyBC()
  end subroutine IC
  real function fric(u)
    use gridMod,  only: dv
    use fieldMod, only: field
    use mympiMod, only: mympi_sum
    type(field),intent(in) :: u
    real,pointer :: up(:,:,:)
    real :: area=0
    if(area==0) then
      area = sum(real(dv(:,4,:),8))
      call mympi_sum(area)
    end if
    up => u%point()
    fric = sum(real(up(:,4,:)*dv(:,4,:),8))/area
    call mympi_sum(fric)
    fric = sqrt(nu*fric/0.5)
  end function fric
end program channel
