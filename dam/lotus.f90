program dambreak
  use gridMod, only : xg
  use twophaseMod, only: twophase
  use imageMod, only: display
  use bodyMod, only: body
  use geom_shape
  use mympiMod
  implicit none
  integer,parameter :: n(3) = (/32,12,5/)*16
  real,parameter :: l = 160, w = 1.228*l, h=0.55*l, tau = 500, g = l/tau**2
  type(twophase) :: a
  type(body) :: box
  integer :: b(3) = (/8,2,1/)
  integer :: cnt=0
  logical :: root,there=.false.
  real    :: tot

  call init_mympi(3,set_blocks=b)
  root = mympi_rank()==0

  xg(1)%h = 3.22*l/real(n(1))
  call xg(2)%stretch(n(2),0.,0.,0.8*l,2*l,prnt=root)

  if(root) print *,'--- 2D dam break --- ',tau
  box = plane(norm=(/ 1,0,0/),center=l*(/0.825,0.,0./)).and. &
        plane(norm=(/-1,0,0/),center=l*(/0.664,0.,0./)).and. &
        plane(norm=(/0,0,-1/),center=l*(/0.,0.,0.295/)).and. &
        plane(norm=(/0, 1,0/),center=l*(/0.,0.161,0./))

  call a%init(n/b,geom=box,nu=1e-2)
  a%g = (/0.,-g,0./)

  call a%f%eval(initf)
  call display(a%f%field,'vof')
  call a%write(box)
!
! -- Update
  do while (a%time<10*tau .and..not.there)
    a%dt = min(1.,a%dt)
    call a%update
    write(9,'(f10.4,f8.4,3e16.8)') a%time/tau,a%dt, &
          box%pforce(a%pressure)/(g*l**3)
    if(floor(a%time/(0.1*tau))>cnt) then
      cnt = cnt+1
      tot = a%f%integral(pow=1)/(w*h*n(3))
      if(root) print '(3f8.4)',a%time/tau,a%dt,tot
      call display(a%f%field,'vof')
      call a%write
    end if
    inquire(file='.kill', exist=there)
  end do
  call mympi_end
contains
  real pure function initf(x) result(p)
    real,intent(in) :: x(3)
    p = merge(0, 1, x(1)<n(1)-w .or. x(2)>h)
  end function initf
end program dambreak
