program gust_model
  use fluidMod, only: fluid
  use bodyMod,  only: body,bodyUnion
  use imageMod, only: display
  use geom_shape
  use gridMod
  implicit none
  type(fluid)        :: flow
  type(set)          :: foil,walls
  type(bodyUnion)    :: geom
  real,parameter     :: c = 64
  real,parameter     :: Re = 1e3
  logical,parameter  :: kinematic = .KINEMATIC_FLAG.
  real,parameter     :: delay = 2
  real,parameter     :: k = K_VAL, f = k/(pi*c)
  real,parameter     :: v = tan(ALPHA_VAL*pi/180.)
  integer,parameter  :: n(3) = (/8.*c,8.*c,1./)
  real               :: x0=0, y0=0, p0=0, dotx=1, doty=0, dotp=0, force(3), moment(3)
  logical            :: there = .false., ext(-3:3) = .true.
!
! -- Set up grid geom and flow
  xg(1)%left = -n(1)/3.
  xg(2)%left = -n(2)/2.; xg(2)%right = xg(2)%left+n(2)

  foil = cylinder(axis=3, radius=c/2, center=0.).map.init_scale(2,w)
  if(kinematic) foil = foil.map.init_rigid(6,p,dp)
  call geom%add(foil)

  walls = (plane(norm=(/0,1,0/),center=(/0.,xg(2)%left+2,0./)).or. &
           plane(norm=(/0,-1,0/),center=(/0.,xg(2)%right-2,0./)).or. &
           plane(norm=(/1,0,0/),center=(/xg(1)%left+2,0.,0./))) &
           .map.init_velocity(gust_velo)
  if(.not.kinematic) call geom%add(walls)

  ext((/-2,-1,2/)) = kinematic
  call flow%init(n,geom,V=(/1.,0.,0./),nu=c/Re,external=ext)
!
! -- Initialize the velocity
  flow%time = -delay*c
  if(.not.kinematic) then
    call flow%velocity%eval(gust_velo)
    call flow%reset_u0()
  end if
!
! -- Time update loop
  do while (flow%time<1./f+c.and..not.there)  ! run 3 cycles
     if(kinematic) call gust_kinematics(flow%time+flow%dt)
     call geom%update(flow%time+flow%dt)
     call flow%update(geom,V=(/dotx,-doty,0./))
     force = -2.*geom%bodies(1)%pforce(flow%pressure)/c
     moment = -2.*geom%bodies(1)%pmoment(flow%pressure)/c**2
     write(9,'(f10.4,f8.4,6e16.8)') flow%time*f,flow%dt,&
        force(:2),moment(3),x0/c,y0/c,p0
     flush(9)
     if(mod(abs(flow%time),0.25*c)<flow%dt) &
        print *,flow%time/c,flow%time*f,x0/c,y0/c,p0
     if(flow%time>0 .and. mod(flow%time,0.125/f)<flow%dt) &
        call display(flow%velocity%vorticity_Z(),'flow',lim=0.25)
     inquire(file='.kill', exist=there)
  end do
  call flow%write()
contains

   real pure function gust(t)
     real,intent(in) :: t
     gust = 0.5*v*(1.-cos(2.*pi*f*t))
     if(t<0 .or. t*f>1) gust = 0.
   end function gust
   subroutine gust_kinematics(t)
     real,intent(in) :: t
     real :: up,vp
     dotp = (gust(flow%time-c)-gust(flow%time))/c  ! average dv/dx
     up = 1+dotp; vp = gust(flow%time-c/2.)-dotp
     dotx = cos(p0)*up+sin(p0)*vp
     doty = cos(p0)*vp-sin(p0)*up
     x0 = x0+flow%dt*dotx
     y0 = y0+flow%dt*doty
     p0 = p0+flow%dt*dotp
   end subroutine gust_kinematics
   pure function gust_velo(x) result(g)
     real,intent(in) :: x(3)
     real            :: g(3)
     g = (/1.,-gust(flow%time-(x(1)+0.5*c)),0./)
   end function gust_velo
   real(8) pure function y(t)
     real(8),intent(in) :: t
     y = y0
   end function y
   real(8) pure function p(t)
     real(8),intent(in) :: t
     p = p0
   end function p
   real(8) pure function dy(t)
     real(8),intent(in) :: t
     dy = doty
   end function dy
   real(8) pure function dp(t)
     real(8),intent(in) :: t
     dp = dotp
   end function dp
   real(8) pure function w(t)
     real(8),intent(in) :: t
     w = 0.12
   end function w

end program gust_model
