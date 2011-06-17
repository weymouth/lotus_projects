!-------------------------------------------------------!
!-------------- Box yank geometry routines -------------!
!-------------------------------------------------------!
module geom
  use analytic
  implicit none
  private
  public init_geom,geom_update,geom_body,geom_fint,geom_velo
  public geom_domain_set,geom_domain
  integer    :: d_save
  type(set)  :: body,fint

  real(8)    :: alpha=0.174,w=0.12,dh,dt,t0
contains
!
!-------------------------------------------------------!
!
  subroutine init_geom
    integer :: i
    fint = .set.plane(1,(/0,0,1/),0,0,0);
    dh = 100.D0/70.D0
    dt = 0.5
    t0 = 7.0
  end subroutine init_geom
!
!-------------------------------------------------------!
!
  subroutine geom_update(time)
    real(8),intent(in) :: time
    real(8),dimension(3) :: fcorner,bcorner,speed
    real(8) :: ca,sa,h,s,pi
    ca = cos(alpha); sa = sin(alpha); pi = acos(-1.D0)

    if(time.le.t0) then
       h = 0
       s = 0
    else if(time.lt.t0+dt) then
       h = dh*(1-cos((time-t0)/dt*pi))*0.5
       s = dh/dt*pi*sin((time-t0)/dt*pi)*0.5
    else
       h = dh
       s = 0
    end if

    speed = (/0.D0,0.D0,s/)
    fcorner = (/-ca*0.5+w*sa*0.5,-sa*0.5-w*ca*0.5,h/)
    bcorner = (/ ca*0.5-w*sa*0.5, sa*0.5+w*ca*0.5,h/)

    body = .set.plane(1,(/ sa,-ca,0.D0/),fcorner,0,speed) &
      .and..set.plane(1,(/-ca,-sa,0.D0/),fcorner,0,speed) &
      .and..set.plane(1,(/0,0,-1/),fcorner,0,speed) &
      .and..set.plane(1,(/-sa, ca,0.D0/),bcorner,0,speed) &
      .and..set.plane(1,(/ ca, sa,0.D0/),bcorner,0,speed)
    return
  end subroutine geom_update
!
!-------------------------------------------------------!
!
  subroutine geom_body(point,distance,normal,velocity,flag,kappa)
    use global, only: ndims
    real(8),intent(in)  :: point(3)
    integer,intent(out) :: flag
    real(8),intent(out) :: distance,normal(ndims),velocity(ndims)
    real(8),intent(out),optional :: kappa
    type(prop) :: a
!
! -- get
    a = body.at.point
!
! -- copy
    distance = a%distance
    normal   = a%normal(:ndims)
    velocity = a%velocity(:ndims)
    if(present(kappa)) &
       kappa = a%kappa(ndims-1)
    flag = a%flag

  end subroutine geom_body
!
!-------------------------------------------------------!
!
  subroutine geom_fint(point,distance,normal)
    real(8),intent(in)  :: point(3)
    real(8),intent(out) :: distance,normal(3)
    type(prop) :: a
!
! -- get
    a = fint.at.point
!
! -- copy
    distance = a%distance
    normal   = a%normal

  end subroutine geom_fint
!
!-------------------------------------------------------!
!
  pure real(8) function geom_velo(d,point)
    integer,intent(in)  :: d
    real(8),intent(in)  :: point(3)
    type(prop) :: a

    geom_velo = 0

  end function geom_velo
!
!-------------------------------------------------------!
!
  subroutine geom_domain_set(d)
    integer,intent(in) :: d
    d_save = abs(d)
  end subroutine geom_domain_set
!
  elemental subroutine geom_domain(x,y,z,delta,velocity)
    real(8),intent(in)    :: x,y,z
    real(8),intent(inout) :: delta,velocity
    type(prop) :: a
!
! -- copy
    delta = 1.
    velocity = 0.

  end subroutine geom_domain
!
! -------------------------------------------------------
end module geom
