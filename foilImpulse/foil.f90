!-------------------------------------------------------!
!------------------ Foil Yank Geom ---------------------!
!-------------------------------------------------------!
!!$module mympi
!!$contains
!!$  integer function mympi_id() result(y)
!!$    y=0
!!$  end function mympi_id
!!$end module mympi

module geom
  use nurbs
  type(nurbs_trim),allocatable,dimension(:) :: mytrim
  real(8),dimension(3) :: rmin,rmax
  integer :: ntrim
  real(8) :: h,s,eps
contains
!
  subroutine init_geom
    use iges
    use mympi, only : mympi_id
    implicit none
    integer :: i,s,s2,myid,ierr
    real(8) :: val(3),cent(3)
    logical :: debug=.false.
    character(100) :: filename
!
! -- open IGS file and convert to scratch files
    myid = mympi_id()
    open(7,file='inp.IGS',IOSTAT=ierr)
    if(ierr.ne.0) stop 'geom: missing inp.IGS'
    call iges_to_srf(7,20,ntrim)
    if(ntrim.ne.144) stop 'no 144 in IGS'
    close(7)
!
! -- read pointer file and allocate
    read(20,*) ntrim
    allocate(mytrim(ntrim))
!
! -- init every surface on every proc
    do s=1,ntrim
       read(20,*) s2
       call init_trim(s2,mytrim(s))
    end do  
!
! -- normalize the shape (scale, shift, rotate, flip)
!
! -- scale
    open(7,file='inp.geo',IOSTAT=ierr)
    if(ierr.ne.0) stop 'geom: missing inp.geo'
    read(7,*,END=100) eps
    if(debug.and.myid.eq.0) write(*,'("eps",3e12.4)') eps
    read(7,*,END=100) val
    if(debug.and.myid.eq.0) write(*,'("inp scale",3e12.4)') val
    call geom_extents
    if(debug.and.myid.eq.0) write(*,'("inp extents",3e12.4)') rmax-rmin
    if(all(val.ne.0.D0)) then       ! scale each dimension
       val = val/(rmax-rmin)
    else if(all(val.eq.0.D0)) then  ! don't scale
       val = 1.D0
    else                            ! scale uniformly
       i = sum(maxloc(val))
       val = val(i)/(rmax(i)-rmin(i))
    end if
    if(debug.and.myid.eq.0) write(*,'("scale",3e12.4)') val
    do s=1,ntrim
       call trim_scale(val,mytrim(s))
    end do    
    call geom_extents
    if(debug.and.myid.eq.0) write(*,'("extents",3e12.4)') rmax-rmin
!
! -- shift
    read(7,*,END=100) cent
    if(debug.and.myid.eq.0) write(*,'("inp cent",3e12.4)') cent
    val = cent-(rmin+rmax)*0.5
    if(debug.and.myid.eq.0) write(*,'("shift",3e12.4)') val
    do s=1,ntrim
       call trim_shift(val,mytrim(s))
    end do    
!
! -- rotate
    read(7,*,END=100) val
    if(debug.and.myid.eq.0) write(*,'("inp rotate",3e12.4)') val
    val = val/180*acos(-1.D0)
    close(7)
    if(debug.and.myid.eq.0) write(*,'("rotate",3e12.4)') val
    do s=1,ntrim
       call trim_rotate(cent,val,mytrim(s))
    end do
    call geom_extents ! get new extents!
    if(debug.and.myid.eq.0) write(*,'("extents",3e12.4)') rmax-rmin
!
! -- flip
    do s=1,ntrim
       call trim_init_grid(6,mytrim(s))
    end do
    if(debug.and.myid.eq.0) print *,'starting array distance'
    call trim_array_distance(cent,mytrim,val(1))
    if(debug.and.myid.eq.0) write(*,'("dis to cent",e12.4)') val(1)
    if(val(1).gt.0) then
       do s=1,ntrim
          call trim_flip(mytrim(s))
       end do
    end if
!
! -- plot
    open(7,file='surf.dat')
    do s=1,ntrim
       if(myid.eq.0) call trim_plot(7,mytrim(s))
    end do
    close(7)
    if(debug) stop 'fgeo debug stop'
    return
100 stop 'geom: end of inp.geo'
  end subroutine init_geom
!
!-------------------------------------------------------!
!
! -- get extents
  subroutine geom_extents
    implicit none
    integer :: s
    rmin = 1e12; rmax = -rmin
    do s=1,ntrim
       rmax = max(rmax,trim_max(mytrim(s)))
       rmin = min(rmin,trim_min(mytrim(s)))
    end do
  end subroutine geom_extents
!
!-------------------------------------------------------!
!
  subroutine geom_update(time)
    real(8),intent(in) :: time
    real(8) :: dh,dt,t0,pi
    dh = 100.D0/70.D0
    dt = 0.5
    t0 = 7.0
    pi = acos(-1.D0)

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

  end subroutine geom_update
!
!-------------------------------------------------------!
!
  subroutine geom_body(ri,distance,normal,velocity,flag,kappa)
    implicit none
    real(8),intent(in)  :: ri(:)
    real(8),intent(out) :: distance,normal(:),velocity(:)
    integer,intent(out) :: flag
    real(8),intent(out),optional :: kappa
    real(8) :: r(3)
!
! -- Get ~3D point and set knowns
!
    velocity = (/0.D0,0.D0,s/)
    r = ri-(/0.D0,0.D0,h/)
    if(present(kappa))then
       kappa = 0. ! incorrect, but ok for now
    end if
    flag = 1 ! solid
!
! -- Check if the point is in range
!
    if(any(r.gt.rmax+3*eps.or.r.lt.rmin-3*eps)) then
       flag = 0
       distance = 3*eps
       normal = 0
       return
    end if
!
! -- If so, get distance and normal on closest surface
!
    call trim_array_distance(r,mytrim,distance,normal)
  end subroutine geom_body
!
!-------------------------------------------------------!
!
  subroutine geom_fint(point,distance,normal)
    real(8),intent(in)  :: point(3)
    real(8),intent(out) :: distance,normal(3)
    distance = 10
    normal   = 0
  end subroutine geom_fint
!
!-------------------------------------------------------!
!
  pure real(8) function geom_velo(d,point)
    integer,intent(in)  :: d
    real(8),intent(in)  :: point(3)
    geom_velo = 0
  end function geom_velo
!
!-------------------------------------------------------!
!
  subroutine geom_domain_set(d)
    integer,intent(in) :: d
  end subroutine geom_domain_set
!
  elemental subroutine geom_domain(x,y,z,delta,velocity)
    real(8),intent(in)    :: x,y,z
    real(8),intent(inout) :: delta,velocity
    delta = 1
    velocity = 0
  end subroutine geom_domain
!
! -------------------------------------------------------
end module geom
!!$program test
!!$  use geom
!!$  call init_geom
!!$end program test
