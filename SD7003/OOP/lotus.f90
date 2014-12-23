!-------------------------------------------------------!
!------------------- SD7003 Test case ------------------!
!-------------------------------------------------------!
program SD7003Test
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use gridMod,    only: xg
  use mympiMod
  implicit none
  real,parameter     :: f = 2              ! resolution
  real,parameter     :: c = 200/f          ! resolution
  integer,parameter  :: n = 9*2**7/f       ! number of points  
  real,parameter     :: Re = 1e4           ! Reynolds number
  logical,parameter  :: O2 = .true.       ! O(2) flag
  real,parameter     :: tStop=30,dtPrint=1
!
  integer            :: b = 4
  real               :: nu = c/Re,t=0,force(3)=0,dt,u,u0
  type(fluid)        :: flow
  type(body)         :: foil
#if MPION
  call init_mympi(2,set_blocks=(/b,b/))
#else
  b=1
#endif
!
! -- Initialize grid and print stats
  call xg(1)%init(n,0.7*c,1.0*c,0.8,f=f)
  call xg(2)%init(n,0.2*c,0.2*c,1.0,f=f)
  if(mympi_rank()==0) then
     print *, '-- Foil test case --'
     print '("    c/h=",i0,", nu=",f0.4,", y+=",f0.4)', int(c),nu,sqrt(0.026/Re**(1./7.)/2.)/nu
     call xg(1)%write()
     call xg(2)%write()
  end if
!
! -- Initialize the foil and flow
  foil = SD(c,4.)
  call flow%init((/n/b,n/b,1/),foil,nu=nu)
  call flow%write(foil)
  if(mympi_rank()==0) print *, '-- init complete --'
!
! -- Turn on/off O(2) !!!
  flow%body = O2
!
! -- Time update loop
  do while (t<tStop)
     dt = flow%dt
     u  = min(1.,(t+dt/c)/5.)
     u0 = min(1.,t/5.)
     flow%velocity%e(1)%bound_val = u
     flow%g(1) = (u-u0)/dt
!
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
    top%n = (/256,64,1/)
    top%xmax(1) = 2.0*chord
    top%xmin(2) = -0.25*chord
    bot = top
    bot%file = 'SD7003_bottom.stl'
    eps = 2.0
    SD = model_init(top).and.model_init(bot)
  end function SD
end program SD7003Test
!
