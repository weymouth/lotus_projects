!-------------------------------------------------------!
!---------------- Abreast square cylinders -------------!
!-------------------------------------------------------!
program square_abreast
  use fluidMod,   only: fluid
  use fieldMod,   only: field
  use bodyMod,    only: body
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank,mympi_coords
  use gridMod,    only: xg,composite
  use imageMod,   only: display
  use geom_shape
  implicit none
  real,parameter :: L = 32                   ! length scale
  real,parameter :: Re = 100                 ! Reynolds number
  real,parameter :: g = 4                    ! gap
!
  real,parameter :: nu = L/Re                ! viscosity
  real,parameter :: xf = 4, yf = g+1
  real,parameter :: m(3) = (/(xf+6)*L,(yf+2)*L,1./)  ! points (approx)
  integer        :: b(3) = (/1,1,1/),box(4)=(/-L,1.5*L,2*L,2*L/)         ! blocks
  real           :: t0=0,Tau=5E2,g_var=2*L
!
  integer        :: n(3)
  logical        :: root,there=.false.,process
  type(fluid)    :: flow
  type(set)      :: square
  type(body)     :: bodies
  type(field)    :: omega
!
! -- Initialize MPI (if MPI is OFF, b is set to 1)
  call init_mympi(2,set_blocks=b(:2))
  root = mympi_rank()==0
!
! -- Print run info
  if(root) print *, '-- Abreast Square Cylinders --'
!
! -- Initialize array size
  n = composite(m,prnt=root)
!
! -- Initialize and print grid
  call xg(1)%stretch(n(1),-15*L,-L,xf*L,15*L,h_max=5.,prnt=root)
  call xg(2)%stretch(n(2),0.,0.,yf*L,(yf+15)*L,prnt=root)
!
! -- Initialize the square geometry and maps
  ! square = plane(4,1,(/0,1,0/),(/0.,0.5*L,0./),0,0) &
  !      .and.plane(4,1,(/0,-1,0/),(/0.,-0.5*L,0./),0,0) &
  !      .and.plane(4,1,(/-1,0,0/),(/-0.5*L,0.,0./),0,0) &
  !      .and.plane(4,1,(/1,0,0/),(/0.5*L,0.,0./),0,0)
  square = cylinder(1,1,3,L/2.,0,0,0)
  bodies = (square.map.init_rigid(2,top,zip)) &
       .or.(square.map.init_rigid(2,bot,zip))

  ! print *,square%at(real((/0,0,0/),8))
  ! print *,square%at(real((/0.5*L,0.,0./),8))
  ! print *,square%at(real((/0.5*L,0.5*L,0./),8))
  ! print *,square%at(real((/-0.5*L,0.5*L,0./),8))
  ! print *,square%at(real((/L,0.,0./),8))
!
! -- Initialize fluid
  call flow%init(n/b,bodies,V=(/1.,0.,0./),nu=nu)
  if(root) print *,'gap=',2*g_var/L
  if(root) print *,'time=',flow%time/L,': resetting to zero'
  flow%time = 0

  ! call display(flow%mu0%e(1),'02xmu0',box=box)
  ! call display(flow%mu0%e(2),'03ymu0',box=box)
  ! call display(flow%mu1(1)%e(1),'04xmu1',box=box)
  ! call display(flow%mu1(2)%e(2),'05ymu1',box=box)
  ! call flow%write(bodies)
!
! -- Time update loop
  if(root) print *, '-- Init complete  --'
  do while (.not.there .and. flow%time<100*L+flow%dt)
!
! -- measure stuff
     if(mod(flow%time,0.1*L)<flow%dt) then
       write(9,'(f10.4,f8.4,3e14.6)') flow%time/L,flow%dt,-2./L*bodies%pforce(flow%pressure)
       flush(9)
      !  omega = flow%velocity%vorticity_Z()
      !  call measure(omega)
       if(mod(flow%time,50*L)<flow%dt) then
         call flow%write()
       omega = flow%velocity%vorticity_Z()
         call display(omega,'03vort',lim=4./L)
       end if
!
! -- update the body and flow
       g_var = 2*L-flow%time/Tau
       if(root) print *, flow%time/L,2*g_var/L
       call bodies%update(flow%time)
     end if
     call flow%update(bodies)
     inquire(file='.kill', exist=there)
  end do
  if(root) print *, '--- complete ---'
  call mympi_end()
contains
!
! -- motion definitions
  real(8) pure function top(t)
    implicit none
    real(8),intent(in) :: t
    top = g_var+L/2.
  end function top
  real(8) pure function bot(t)
    implicit none
    real(8),intent(in) :: t
    bot = -(g_var+L/2.)
  end function bot
  real(8) pure function zip(t)
    implicit none
    real(8),intent(in) :: t
    zip = 0.
  end function zip
!
! -- measurement routines
  subroutine measure(omega)
    use gridMod, only: dv
    type(field),intent(in) :: omega
    type(field) :: hold,x,y
    real :: neg(n(1)),pos(n(1)),circ,xcen,ycen,coord(3)
    integer :: i,ns,ps
    real,pointer :: h(:,:,:) => null()
! -- positive vorticity only
    hold = omega; h => hold%point()
    h = h*dv
    where(h<0) h = 0
! -- coords
    y = omega; x = omega
    do i = y%js-1,y%je+1
      coord = y%pos(1,i,1); y%p(:,i,:) = coord(2)
    end do
    do i = y%is-1,y%ie+1
      coord = x%pos(i,1,1); x%p(i,:,:) = coord(1)
    end do
! -- integrate vorticity over y
    rewind(10)
    do i=1,n(1)
        neg(i) = integral(hold,i,i,coord=0)
        pos(i) = integral(hold,i,i,coord=1)
        write(10,'(3f8.4)') xg(1)%x(i)/L,neg(i),pos(i)
    end do
! -- find size and centriod of vorticies
    ns = 1; ps = 1
    do i=2,n(1)-1
        if(neg(i)<neg(i-1).and.neg(i)<neg(i+1)) then
          circ = sum(neg(ns:i))
          if(circ>0.1*L) then
            ycen = integral(hold,ns,i,coord=0,weight=y)/circ
            xcen = integral(hold,ns,i,coord=0,weight=x)/circ
            write(11,'(4f8.4)') xcen/L,ycen/L,circ/L,(xg(1)%x(i)-xg(1)%x(ns))/L
            ns = i
          end if
        end if
        if(pos(i)<pos(i-1).and.pos(i)<pos(i+1)) then
          circ = sum(pos(ps:i))
          if(circ>0.1*L) then
            ycen = integral(hold,ps,i,coord=1,weight=y)/circ
            xcen = integral(hold,ps,i,coord=1,weight=x)/circ
            write(12,'(4f8.4)') xcen/L,ycen/L,circ/L,(xg(1)%x(i)-xg(1)%x(ps))/L
            ps = i
          end if
        end if
    end do
    flush(10);flush(11);flush(12)
  end subroutine measure
  real function integral(a,gs,ge,coord,weight) result(ai)
    use mympiMod, only: mympi_coords,mympi_sum
    type(field),intent(in) :: a
    integer,intent(in) :: gs,ge,coord
    type(field),intent(in),optional :: weight
    integer :: ijk(3),ls,le

    ijk = a%local((/gs,1,1/)); ls = min(max(ijk(1),a%is),a%ie+1)
    ijk = a%local((/ge+1,1,1/)); le = min(max(ijk(1),a%is),a%ie+1)-1
    if(present(weight)) then
      ai = merge(sum(a%p(ls:le,a%js:a%je,1)*weight%p(ls:le,a%js:a%je,1)), &
            0.,mympi_coords(2)==coord)
    else
      ai = merge(sum(a%p(ls:le,a%js:a%je,1)),0.,mympi_coords(2)==coord)
    end if
    call mympi_sum(ai)
  end function integral
end program square_abreast
