program sphere_flow
  use mympiMod,   only: init_mympi
  use gridMod,    only: xg
  use geom_shape, only: sphere,operator(.map.),init_scale,pi
  use bodyMod,    only: body
  use fluidMod,   only: fluid
  use imageMod,   only: display
  implicit none
!
! -- Define parameters, declare variables
  real,parameter    :: D = 64                           ! sphere size in cells
  real,parameter    :: Re_D = 30e3                      ! Reynolds number based on D
  integer,parameter :: n(3) = D*(/4,2,2/)               ! number of cells in ijk
  integer           :: b(3) = (/2,1,1/)                 ! MPI domain cuts in ijk
  real,parameter    :: lim = 10./D                      ! vorticity level in image
  real,parameter    :: dprnt = 0.05                     ! how often to print
  integer,parameter :: box(4) = D*(/-1,-1,4,2/)         ! image size in pixels
  logical     :: there = .false.                        ! flag for stopping
  type(fluid) :: flow
  type(body)  :: geom
!
! -- Initialize grid, body and fluid
  call init_mympi(ndims=3,set_blocks=b)
  call xg(1)%stretch(n(1),-2*D,-D/2,D/2,5*D,h_max=8.)
  call xg(2)%stretch(n(2),-2*D,-D/2,D/2,2*D)
  call xg(3)%stretch(n(3),-2*D,-D/2,D/2,2*D)
  geom = sphere(2, 1, radius=0.5*D, center=0).map.init_scale(1,length,rate)
  call flow%init(n/b, geom, V=(/1.,0.,0./), nu=D/Re_D)
!
! -- Time update loop
  do while(flow%time<2*D .and. .not.there)
    if(flow%time<D+flow%dt) call geom%update(flow%time/D)
    call flow%update(geom)
    write(9,'(f10.4,f8.4,3e16.8)') &
        flow%time/D,flow%dt,-2.*geom%pforce(flow%pressure)/(3.14159*D**2/4)
    flush(9)
    if(mod(flow%time, D*dprnt)<flow%dt) then
      call display(flow%velocity%vorticity_Z(), 'out_vort', lim=lim, box=box)
    end if
    inquire(file='.kill', exist=there)
  end do
  call flow%write()

contains
  real(8) pure function length(ts)  ! length scale
    real(8),intent(in) :: ts
    length = 1-log(cosh(pi))/pi
    if(ts<1) length = length+log(cosh(pi*(ts-1)))/pi
  end function length
  real(8) pure function rate(ts)  ! rate of change of length
    real(8),intent(in) :: ts
    rate = tanh(pi*(ts-1))/D
    if(ts>1) rate = 0
  end function rate
end program sphere_flow
