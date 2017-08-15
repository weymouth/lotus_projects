program gust_model
  use fluidMod, only: fluid
  use bodyMod,  only: body,bodyUnion
  use imageMod, only: display
  use mympiMod
  use geom_shape
  use gridMod
  implicit none
  type(fluid)        :: flow
  type(set)          :: foil,walls
  type(bodyUnion)    :: geom
  real,parameter     :: c = 90
  real,parameter     :: Re = 1e3
  character(8),parameter  :: run = 'RUN_FLAG'
  real,parameter     :: delay = 2
  real,parameter     :: k = K_VAL, f = k/(pi*c)
  real,parameter     :: v = tan(ALPHA_VAL*pi/180.)
  integer,parameter  :: n(3) = (/8.*c,8.*c,1./)
  real               :: x0=0, y0=0, p0=0, dotx=1, doty=0, dotp=0, force(3), moment(3)
  logical            :: there = .false., ext(-3:3) = .true., root
  integer            :: bl(3) = (/4,4,1/), box(4) = (/-2.67*c,-4*c,8*c,8*c/)
!
! -- Set up grid geom and flow
  call init_mympi(2,set_blocks=bl(:2))
  root = mympi_rank()==0

  call xg(1)%stretch(n(1),-2.67*c,-0.6*c,5.3*c,10*c,prnt=there.and.root)
  call xg(2)%stretch(n(2),-10*c,-3*c,3*c,10*c,prnt=there.and.root)

  foil = plane(norm=(/0,1,0/),center=(/0.,2.5,0./)).and. &
         plane(norm=(/0,-1,0/),center=(/0.,-2.5,0./)).and. &
         plane(norm=(/1,0,0/),center=(/c/2.,0.,0./)).and. &
         plane(norm=(/-1,0,0/),center=(/-c/2.,0.,0./))
  if(run.ne.'gust') foil = foil.map.init_rigid(6,p,dp)
  call geom%add(foil)

  walls = (plane(norm=(/0,1,0/),center=(/0.,xg(2)%x(2),0./)).or. &
           plane(norm=(/0,-1,0/),center=(/0.,xg(2)%x(n(2)),0./)).or. &
           plane(norm=(/1,0,0/),center=(/xg(1)%x(2),0.,0./))) &
           .map.init_velocity(gust_velo)
  if(run.eq.'gust') call geom%add(walls)
  ext((/-2,-1,2/)) = run.ne.'gust'

  call flow%init(n/bl,geom,V=(/1.,0.,0./),nu=c/Re,external=ext)
!
! -- Initialize the velocity
  flow%time = -delay*c
  if(run.eq.'gust') then
    call flow%velocity%eval(gust_velo)
    call flow%reset_u0()
  end if
!
! -- Time update loop
  do while (flow%time<1./f+c.and..not.there)  ! run 3 cycles
     call gust_kinematics(flow%time+flow%dt)
     call geom%update(flow%time+flow%dt)
     call flow%update(geom,V=(/dotx,-doty,0./))
     force = -2.*geom%bodies(1)%pforce(flow%pressure)/c
     moment = -2.*geom%bodies(1)%pmoment(flow%pressure)/c**2
     write(9,'(f10.4,f8.4,7e16.8)') flow%time*f,flow%dt,&
        force(:2),moment(3),flow%time/c,x0/c,y0/c,p0
     flush(9)
     if(root.and.mod(abs(flow%time),c)<flow%dt) &
        print *,flow%time/c,flow%time*f
     if(flow%time>0 .and. mod(flow%time,0.125/f)<flow%dt) &
        call display(flow%velocity%vorticity_Z(),'flow',lim=0.25,box=box)
     inquire(file='.kill', exist=there)
  end do
  call display(flow%velocity%vorticity_Z(),'flow',lim=0.25,box=box)
  call flow%write()
  call mympi_end
contains

   real pure function gust(t)
     real,intent(in) :: t
     gust = 0.5*v*(1.-cos(2.*pi*f*t))
     if(t<0 .or. t*f>1) gust = 0.
   end function gust
   subroutine gust_kinematics(t)
     real,intent(in) :: t
     real :: vp, g = sqrt(1./3.)
     select case(run)
        case('gust')
          return
        case('edges') ! leading and trailing edge matching
          dotp = (gust(t-c)-gust(t))/c
          vp = (gust(t-c)+gust(t))/2.
        case('center') ! center matching
          dotp = (gust(t-c/2.-1.)-gust(t-c/2.+1.))/2.
          vp = gust(t-c/2.)
        case('gauss') ! center matching
          dotp = (gust(t-c/2.*(1.+g))-gust(t-c/2.*(1.-g)))/(c*g)
          vp = (gust(t-c/2.*(1.+0.))*0.5688888888888889 &
               +gust(t-c/2.*(1.-0.5384693101056831))*0.4786286704993665 &
               +gust(t-c/2.*(1.+0.5384693101056831))*0.4786286704993665 &
               +gust(t-c/2.*(1.-0.9061798459386640))*0.2369268850561891 &
               +gust(t-c/2.*(1.+0.9061798459386640))*0.2369268850561891)/2.
        case('heave') ! heave only
          dotp = 0
          vp = gust(flow%time-c/2.)
        case default
          if(root) print *,'bad RUN_FLAG:',run
          there=.true.
     end select
     dotx = cos(p0)+sin(p0)*vp
     doty = cos(p0)*vp-sin(p0)
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
