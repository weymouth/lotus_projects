!-------------------------------------------------------!
!------------------- SD7003 Test case ------------------!
!-------------------------------------------------------!
program SD7003Test
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use gridMod,    only: xg
  use mympiMod
  implicit none
  real, parameter    :: c = 100            ! resolution
  integer,parameter  :: n = 9*2**6         ! number of points  
  real,parameter     :: Re = 1e4           ! Reynolds number
  logical,parameter  :: O2 = .false.       ! O(2) flag
  real,parameter     :: tStop=100,dtPrint=1
!
  integer            :: b = 4
  real               :: nu = c/Re,t=0,force(3)=0,dt
  type(fluid)        :: flow
  type(body)         :: foil
#if MPION
  call init_mympi(2,set_blocks=(/b,b/))
#else
  b=1
#endif
!
! -- Initialize grid and print stats
  call xg(1)%init(n,0.7*c,1.0*c,0.8,f=2.)
  call xg(2)%init(n,0.2*c,0.2*c,1.0,f=2.)
  if(mympi_rank()==0) then
     print *, '-- Foil test case --'
     print '("    c/h=",i0,", nu=",f0.4,", y+=",f0.4)', int(c),nu,sqrt(0.026/Re**(1./7.)/2.)/nu
     call xg(1)%write()
     call xg(2)%write()
  end if
!
! -- Initialize the foil and flow
  foil = SD(c,4.)
  call flow%init((/n/b,n/b,1/),foil,V=(/1.,0.,0./),nu=nu)
  call flow%write(foil)
  if(mympi_rank()==0) print *, '-- init complete --'
!
! -- Turn on/off O(2) !!!
  flow%body = O2
!
! -- Time update loop
  do while (t<tStop)
     dt = flow%dt
     call flow%update
!
     t = flow%time/c
     force = 2.*foil%pforce(flow%pressure)/c
     write(9,1) t,dt,force
     flush(9)
!
     if(mod(t,dtPrint)<dt/c) then
        call flow%write
        if(mympi_rank()==0) print 1,t,dt,force
     end if
  end do

  if(mympi_rank()==0) print *,'--- complete --- '
  call mympi_end
1 format(f10.4,f8.4,3e16.8)
contains
  type(set) function SD(chord,alpha)
    use geom_shape
    real,intent(in) :: chord,alpha
    type(model_info) :: top,bot
    top%file = 'SD7003_top.stl'
    top%x = (/-0.5,0.,0./)
    top%r = (/-alpha,0.,0./)
    top%s = chord
    top%n = (/1.5*chord,0.5*chord,10./)
    top%xmax(1) = 1.2*chord
    top%xmin(2) = -0.2*chord
    bot = top
    bot%file = 'SD7003_bottom.stl'
    eps = 2.0
    SD = model_init(top).and.model_init(bot)
  end function SD
end program SD7003Test
!
