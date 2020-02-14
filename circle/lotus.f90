program circle_flow
  use bodyMod,    only: body
  use fluidMod,   only: fluid
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  use imageMod,   only: display
  use geom_shape
  implicit none
!
! -- Define parameter, declare variables
  real,parameter :: D = 32, Re = 100, nu = D/Re ! diameter and Reynolds number
  integer :: n(3), b(3) = (/2,1,1/), box(4) = (/-2*D,-3*D,20*D,6*D/)   ! MPI blocks (product must equal n_procs)
  logical :: root, there = .false.     ! root processor
  type(fluid) :: flow
  type(body) :: geom
  integer :: i
!
! -- Initialize MPI (if MPI is OFF, b is set to 1)
  call init_mympi(2,set_blocks=b)
  root = mympi_rank()==0
!
  if(root) print *,'Setting up the grid, body and fluid'
  if(root) print *,'-----------------------------------'
  n = composite((/8.*D,4.*D,0./),prnt=root)
  call xg(1)%stretch(n(1), -10*D, -0.6*D, 3*D, 20*D, h_max=8.,prnt=root)
  call xg(2)%stretch(n(2), -15*D, -0.9*D, 0.9*D, 15*D, prnt=root)
  geom = cylinder(axis=3, radius=0.5*D, center=0)!.and.plane(norm=(/1,0,0/),center=0.)
  call flow%init(n/b,geom,V=(/1.,0.,0./),nu=nu,exit=.true.)

  call display(flow%velocity%vorticity_Z(), 'vort', lim = 0.25, box=box)
  call SGS(flow%velocity)
!
  ! if(root) print *,'Starting time update loop'
  ! if(root) print *,'-----------------------------------'
  ! if(root) print *,' -t- , -dt- '
  ! flow%time = 0
  ! do while(flow%time<1*D.and..not.there)
  !   call flow%update()
  !   write(9,'(f10.4,f8.4,6e16.8)') flow%time/D,flow%dt, &
  !      2./D*nu*geom%vforce(flow%velocity), &
  !     -2./D*geom%pforce(flow%pressure)
  !     ! (flow%velocity%e(i)%max(),flow%velocity%e(i)%max(neg=.TRUE.),0.,i=1,2)
  !   flush(9)
  !   if(mod(flow%time,1*D)<flow%dt) then
  !     if(root) print '(f6.1,",",f6.3)',flow%time/D,flow%dt
  !   end if
  !   inquire(file='.kill', exist=there)
  ! end do
  ! if(root) print *,'Loop complete: writing restart files and exiting'
  ! if(root) print *,'-----------------------------------'
  ! call flow%write()
  call mympi_end

contains
  subroutine SGS(u)
    use vectorMod, only: vfield
    use fieldMod,  only: field
    type(vfield),intent(in) :: u
    type(vfield) :: uC,uCh,f,uh,fh,fuh
    type(field)  :: a,b
    integer :: d
    real,pointer :: p(:,:,:)
    real :: l
    call uC%init(u%e(1)%size(),l=0)
    call a%init(u%e(1)%size(),P=0.)
    do d=1,u%ndims
      p => uC%e(d)%point()
      p = p+0.5*(u%e(d)%point()*u%e(d)%point(d)) ! u_C
    end do
    uCh = uC%halve()             !<u_C>

    do d=1,u%ndims
      a%p = uC%e(d)%p**2         ! u_C^2
      b = a%halve()              !<u_C^2>
      b%p = b%p-uCh%e(d)%p**2    ! <u_C^2>-<u_C>^2
      call display(b,'SGS',limout=l,box=box)
      print *,d,l
    end do

    a%p = uC%e(1)%p*uC%e(2)%p
    b = a%halve()
    b%p = b%p-uCh%e(1)%p*uCh%e(2)%p
    call display(b,'SGS',limout=l,box=box)
    print *,d,l

    call f%convect_diffuse(u,nu=0.)   ! conv(u)
    fh = f%halve()                    ! <conv(u)>
    uh = u%halve()                    ! <u>
    call fuh%convect_diffuse(uh,nu=0.) ! conv(<u>)
    do d=1,u%ndims
      b%p = fh%e(d)%p-fuh%e(d)%p
      call display(b,'SGS',limout=l,box=box)
      print *,d,l
    end do
  end subroutine SGS

end program circle_flow
