!-------------------------------------------------------!
!---------- Geometry Routines using geom.txt -----------!
!-------------------------------------------------------!
module geom
  use analytic
  use nurbs
  implicit none
  private
  public init_geom,geom_update,geom_body,geom_fint,geom_velo
  public geom_domain_set,geom_domain
  integer    :: d_save
  type(set)  :: body,fint,velo,domain
  type(nurbsSet)  :: bodyIGS
  logical :: useIGS
contains
!
!-------------------------------------------------------!
!
  subroutine init_geom
    use inout
    use mympi
    integer :: i
!
! -- open geometry files
    open(7,file='inp.ana',IOSTAT=i)
    if(i.ne.0) stop 'geom: missing analytic input, inp.ana'
    inquire(file='inp.IGS',exist=useIGS)
    if(useIGS) inquire(file='inp.IGS2',exist=useIGS)
    if(useIGS) then
       open(18,file='inp.IGS')
       open(19,file='inp.IGS2')
    end if
!
! -- read into sets on all processors
    do i=0,mympi_idmx()
       if(mympi_id().eq.i) then
          if(useIGS) then
             call init_nurbs(bodyIGS,18,19)
             rewind(6)
          end if
          call set_read(7,body)
          call set_read(7,fint)
          call set_read(7,velo)
          call set_read(7,domain)
          rewind(7)
       end if
    end do
    close(7);close(6,iostat=i)
!
! -- write to log
    write(io_log,'("  geometry files read")')
    call log_print
  end subroutine init_geom
!
!-------------------------------------------------------!
!
  subroutine geom_update(time)
    real(8),intent(in) :: time
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
    if(useIGS) then
       a = bodyIGS.at.point
    else
       a = body.at.point
    end if
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
!
! -- get
    a = velo.at.point
!
! -- copy
    geom_velo = merge(0,1,a%distance.gt.0)*a%velocity(d)

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
! -- get
    a = domain.at.(/x,y,z/)
!
! -- copy
    delta = 1.
    velocity = merge(0,1,a%distance.gt.0)*a%velocity(d_save)

  end subroutine geom_domain
!
! -------------------------------------------------------
end module geom
