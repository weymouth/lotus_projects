!-------------------------------------------------------!
!------------------- SD7003 Test case ------------------!
!-------------------------------------------------------!
program SD7003Test
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use gridMod,    only: xg
  use imageMod,   only: display
  implicit none
  real,parameter     :: c = 200            ! resolution
  integer,parameter  :: n = 9*2**6         ! number of points
  real,parameter     :: Re = 1e4           ! Reynolds number
  real,parameter     :: tStop=30,dtPrint=1
!
  real               :: nu = c/Re,t=0,force(3),dt,u
  type(fluid)        :: flow
  type(body)         :: foil
!
! -- Initialize grid and print stats
  print *, '-- Foil test case --'
  print 2,int(c),nu,sqrt(0.026/Re**(1./7.)/2.)/nu
2 format("    c/h=",i0,", nu=",f0.4,", y+=",f0.4)
  call xg(1)%stretch(n,-5*c,-0.7*c,1.0*c,10*c,h_min=2.,h_max=6.,prnt=.true.)
  call xg(2)%stretch(n,-10*c,-0.2*c,0.2*c,10*c,prnt=.true.)
!
! -- Initialize the foil and flow
  foil = SD(c,4.)
  call flow%init((/n,n,1/),foil,nu=nu)
  call flow%write(foil)
  print *, '-- init complete --'
!
! -- Time update loop
  do while (t<tStop)
     dt = flow%dt
     u = min(1.,(t+dt/c)/5.)
     call flow%update(V=(/u,0./))
!
! -- write to file
     t = flow%time/c
     force = 2.*foil%pforce(flow%pressure)/c
     write(9,1) t,dt,force
1    format(f10.4,f8.4,3e16.8)
     flush(9)
!
! -- print & display to screen
     if(mod(t,dtPrint)<dt/c) then
        call flow%write
        print 1,t,dt,force
        call display(flow%velocity%vorticity_Z(),'vort',&
          lim=15./c,box=(/-200,-320,2080,640/))
     end if
  end do
  print *,'--- complete --- '
contains
  type(set) function SD(chord,alpha)
    use geom_shape
    real,intent(in) :: chord,alpha
    type(model_info) :: top,bot
    top%file = 'SD7003_top.stl'
    top%x = (/-0.5,0.,0./)
    top%r = (/-alpha,0.,0./)
    top%s = chord
    top%n = (/256,64,1/)
    top%xmax(1) = 2.0*chord
    top%xmin(2) = -0.25*chord
    bot = top
    bot%file = 'SD7003_bottom.stl'
    ! model_fill = .false.
    eps = 2.0
    SD = model_init(top).and.model_init(bot)
  end function SD
end program SD7003Test
!
