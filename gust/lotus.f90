program gust_model
  use fluidMod, only: fluid
  use bodyMod,  only: body,bodyUnion
  use imageMod, only: display
  use geom_shape
  use gridMod
  implicit none
  type(fluid)        :: flow
  type(body)         :: foil,walls
  type(bodyUnion)    :: geom
  real,parameter     :: c = 64
  real,parameter     :: Re = 1e4
  real,parameter     :: k = 0.25, f = k/(pi*c)
  real,parameter     :: v = tan(15.*pi/180.)
  integer,parameter  :: n(3) = (/8.*c,8.*c,1./)
  real               :: y0=0, p0=0, doty=0, dotp=0, force(3), moment(3)
  logical            :: there = .FALSE., wall_flag(-3:3) = .FALSE.
!
! -- Set up geom
  xg(1)%left = -n(1)/3.;
  xg(2)%left = -4*c; xg(2)%right = 4*c
  foil = cylinder(3,c/2.,center=0.).map.init_scale(2,w)! &
        ! .map.init_rigid(6,p,dp).map.init_rigid(2,y,dy)
  call geom%add(foil)
  walls = (plane(norm=(/0,1,0/),center=(/0.,xg(2)%left+2,0./)).or. &
           plane(norm=(/0,-1,0/),center=(/0.,xg(2)%right-2,0./)).or. &
           plane(norm=(/1,0,0/),center=(/xg(1)%left+2,0.,0./))) &
           .map.init_velocity(gust_velo)
  call geom%add(walls)
  wall_flag((/-2,-1,2/)) = .TRUE.
  call flow%init(n,geom,V=(/1.,0.,0./),nu=c/Re,external=.not.wall_flag)
!
! -- Initialize the velocity
  flow%time = -2.*c
  call flow%velocity%eval(gust_velo)
  call flow%reset_u0()
!
! -- Time update loop
  do while (flow%time<1./f+c.and..not.there)  ! run 3 cycles
    !  call gust_kinematics(flow%time+flow%dt)
     call geom%update(flow%time+flow%dt)
     call flow%update(geom)
     force = -2.*geom%bodies(1)%pforce(flow%pressure)/c
     moment = -2.*geom%bodies(1)%pmoment(flow%pressure)/c**2
     write(9,'(f10.4,f8.4,5e16.8)') flow%time*f,flow%dt,&
        force(:2),moment(3),y0/c,p0
     flush(9)
     if(flow%time>0 .and. mod(flow%time,0.25*c)<flow%dt) then
        print *,flow%time/c,flow%time*f
        call display(flow%velocity%vorticity_Z(),'flow',lim=0.25)
     end if
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
     real :: t_LE,t_TE,at_LE,at_TE
     t_LE = tan(atan(gust(t))-p0)
     t_TE = tan(atan(gust(t-c))-p0)
     at_LE = sin(p0)*t_LE-cos(p0); at_TE = sin(p0)*t_TE-cos(p0)
     dotp = (t_LE-t_TE)/(at_LE+at_TE)/(0.5*c)
     doty = (t_TE*at_LE+t_LE*at_TE)/(at_LE+at_TE)
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
