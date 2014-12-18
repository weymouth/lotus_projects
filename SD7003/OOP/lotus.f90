!-------------------------------------------------------!
!------------------- SD7003 Test case ------------------!
!-------------------------------------------------------!
program foil_impulse
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use geom_shape  ! to create geometry
  use gridMod,    only: xg
  implicit none
  real, parameter    :: L = 100            ! resolution
  integer,parameter  :: n = 9*2**6         ! number of points  
  real,parameter     :: Re = 1e4           ! Reynolds number
  real,parameter     :: tStop=30,dtPrint=1
!
  real               :: nu = L/Re,t=0,force(3)=0,dt
  type(fluid)        :: flow
  type(body)         :: foil
  type(model_info)   :: top,bot
!
  call xg(1)%init(n,0.7*L,1.0*L,0.8,r=1.04)
  call xg(1)%write
  call xg(2)%init(n,0.1*L,0.2*L,1.0,r=1.04,h=0.25,c=16.)
  call xg(2)%write
!
! -- Initialize the foil geometry
!  model_fill = .false.
  top%file = 'SD7003_top.stl'
  top%x = (/-0.5,0.,0./)
  top%r = (/-4.0,0.,0./)
  top%s = L
  top%n = (/128,32,10/)
  bot = top
  bot%file = 'SD7003_bottom.stl'
  eps = 2.0
  foil = model_init(top).and.model_init(bot)
!
! -- Initialize fluid
  call flow%init((/n,n,1/),foil,V=(/1.,0.,0./),nu=nu)
!  flow%body = .false.
  print *, '-- init complete --'
!
! -- Time update loop
  do while (t/L<tStop)
!
!-- update and write fluid
     dt = flow%dt
     call flow%update
     t = flow%time
     if(mod(t,dtPrint*L)<dt) call flow%write
!
! -- print force
     force = foil%pforce(flow%pressure)
     write(9,1) t/L,dt,2.*force/L
     print 1,t/L,dt,2.*force/L
1    format(f10.4,f8.4,3e16.8)
     flush(9)
  end do
  write(6,*) '--- complete --- '
end program foil_impulse
!
