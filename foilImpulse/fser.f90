!------------------------------------------------------------!
!----------------- Serial version of MPI Routines -----------!
!------------------------------------------------------------!
module mympi
  integer,private,allocatable,dimension(:) :: coords,points,me,sub_points
  integer,private                          :: myid,myidmx,ios,mpi_grid
  integer,private                          :: gis,gie,gjs,gje,gks,gke
  real(8),private                          :: mpi_time
contains
!----------------------------------------------------------------
!
! -- Initialize a serial run, reading in inputs, computing dims
! setting up the topography
!
  subroutine init_mympi
    use global
    use utility, only: nrerror,log,log_print
    implicit none
    integer :: d
    logical,allocatable,dimension(:) :: periodic
    integer,allocatable,dimension(:) :: blocks
!
! -- Set up default MPI communicator
!
    myid = 0
    call cpu_time(mpi_time)
!
! -- Open input file
!
    open(unit=input_num,file='inp.txt',status='old',iostat=ios)
    if(mympi_any(ios.gt.0)) call nrerror('mympi_init: can not open inp.txt')
!
! -- Read global inputs
!
    call mympi_read(input_num)
    call mympi_read(input_num,i=ndims)
    call mympi_read(input_num,l=lres)
    call mympi_read(input_num,r=dt)
    call mympi_read(input_num,i=tmx)
    call mympi_read(input_num,i=tmd)
    call mympi_read(input_num,i=thk)
    if(ndims.ne.2.and.ndims.ne.3) call nrerror('init_mympi:wrong dims')
    dti = 1./dt
    write(log,'("  ndims,dt=",i2,e12.4)') ndims,dt
    call log_print
!
! -- Read mpi inputs
!
    allocate(blocks(0:ndims-1),points(0:ndims-1),periodic(0:ndims-1))
    call mympi_read(input_num)
    call mympi_read(input_num,ia=points)
    call mympi_read(input_num,la=periodic)
    call mympi_read(input_num,ia=blocks)
    myidmx = 0
    write(log,'("  processors=",i5)') 1
    call log_print
!
! -- Create parrallel topology
!
    allocate(me(-ndims:ndims),coords(0:ndims-1))
    coords = 0
    me = -1 ; me(0) = 0 ; myid = 0; myidmx = 0
    forall(d=1:ndims , periodic(d-1)) me( d) = 0 
    forall(d=1:ndims , periodic(d-1)) me(-d) = 0 
!
! -- Set number of points in each block
!
    ni = points(0)+2*thk; is = 1+thk; ie = ni-thk
    nj = points(1)+2*thk; js = 1+thk; je = nj-thk
    if(ndims.eq.3) then
       nk = points(2)+2*thk; ks = 1+thk; ke = nk-thk
    else if(ndims.eq.2) then
       nk = 1; ks = 1; ke = 1
    end if
    write(log,'("  ni,nj,nk=",3i5)') ni,nj,nk
    call log_print

  end subroutine init_mympi
!
!----------------------------------------------------------------
!
! -- Basic query functions
!
  logical function mympi_domain_bound(d)
    integer,intent(in) :: d
    mympi_domain_bound = me(d).lt.0
  end function mympi_domain_bound
!
  integer function mympi_coords(d)
    integer,intent(in) :: d
    mympi_coords = coords(d-1)
  end function mympi_coords
!
  integer function mympi_points(d)
    integer,intent(in) :: d
    mympi_points = points(d-1)
  end function mympi_points
!
  logical function mympi_out(p)
    use global
    implicit none
    integer,dimension(4),intent(in) :: p
    mympi_out = p(2).lt.is.and.me(-1).lt.0
    mympi_out = p(2).gt.ie.and.me( 1).lt.0.or.mympi_out
    mympi_out = p(3).lt.js.and.me(-2).lt.0.or.mympi_out
    mympi_out = p(3).gt.je.and.me( 2).lt.0.or.mympi_out
    if(ndims.eq.3) then
       mympi_out = p(4).lt.ks.and.me(-3).lt.0.or.mympi_out
       mympi_out = p(4).gt.ke.and.me( 3).lt.0.or.mympi_out
    end if
  end function mympi_out
!
  integer function mympi_id()
    implicit none
    mympi_id = myid
  end function mympi_id
!
  integer function mympi_idmx()
    implicit none
    mympi_idmx = myidmx
  end function mympi_idmx
!
  integer function mympi_comm()
    mympi_comm = mpi_grid
  end function mympi_comm
!
!----------------------------------------------------------------
!
! -- Read and distribute file info (serial version)
!
  subroutine mympi_read(num,l,la,i,ia,r,ra,ru,f)
    use global,  only: print_flags
    use utility, only: nrerror
    implicit none
    integer,intent(in)                  :: num
    logical,intent(inout),optional      :: l,la(:)
    integer,intent(inout),optional      :: i,ia(:)
    real(8),intent(inout),optional      :: r,ra(:),ru
    type(print_flags),intent(out),optional :: f
    character(len=1) :: c
    logical :: debug=.false.
!
! -- read
    ios = 0
    root: if(myid.eq.0) then
       if(present(l)) then
          read(num,*,iostat=ios) l
          if(debug) write(*,*) 'l',l,ios
       else if(present(la)) then
          read(num,*,iostat=ios) la
          if(debug) write(*,*) 'la',la,ios
       else if(present(i)) then
          read(num,*,iostat=ios) i
          if(debug) write(*,*) 'i',i,ios
       else if(present(ia)) then
          read(num,*,iostat=ios) ia
          if(debug) write(*,*) 'ia',ia,ios
       else if(present(r)) then
          read(num,*,iostat=ios) r
          if(debug) write(*,*) 'r',r,ios
       else if(present(ra)) then
          read(num,*,iostat=ios) ra
          if(debug) write(*,*) 'ra',ra,ios
       else if(present(ru)) then
          read(num,iostat=ios) ru
          if(debug) write(*,*) 'ru',ru,ios
       else if(present(f)) then
          read(num,*,iostat=ios) f%prnt,f%rwnd,f%lwrs,&
               f%file,f%tmod,f%skip,f%ghst
          if(debug) write(*,*) 'f',f,ios
       else
          read(num,*,iostat=ios)
          if(debug) write(*,*) 'b ',ios
       end if
!
! -- error check       
       if(ios.ne.0) then
          backspace(num)
          read(num,*,iostat=ios) c
          if(c.eq.'!') then
             backspace(num)
             ios = -1
          else
             ios = 1
          end if
       end if
    end if root
    if(mympi_any(ios.eq.1)) then
       call nrerror('mympi_read: read')
    else if(mympi_any(ios.eq.-1)) then
       return
    end if

  end subroutine mympi_read
!
!----------------------------------------------------------------
!
! -- Sum over proccessors (serial version)
!
  subroutine mympi_sum(each,all,average)
    implicit none
    real(8),intent(in)          :: each
    real(8),intent(out)         :: all
    logical,intent(in),optional :: average
    all = each
    if(present(average)) then
       if(average) all = all/product(points)
    end if
  end subroutine mympi_sum
!
  subroutine mympi_domain_sum(d,each,all)
    implicit none
    integer,intent(in)          :: d
    real(8),intent(in)          :: each
    real(8),intent(out)         :: all
  	all = each
  end subroutine mympi_domain_sum
!
!----------------------------------------------------------------
!
! -- Get the elapsed time
!
  subroutine mympi_timer(elapsed)
    implicit none
    real(8),intent(out) :: elapsed
    real(8) :: time

    call cpu_time(time)
    elapsed = time-mpi_time
    mpi_time = time

  end subroutine mympi_timer
!
!----------------------------------------------------------------
!
! -- Compare across proccessors (serial version)
!
  logical function mympi_any(each)
    implicit none
    logical :: each
    mympi_any = each
  end function mympi_any
!
  logical function mympi_all(each)
    implicit none
    logical :: each
    mympi_all = each
  end function mympi_all
!
  real(8) function mympi_max(each)
    implicit none
    real(8) :: each
    mympi_max = each
  end function mympi_max
!
!----------------------------------------------------------------
!
! -- Fill ghosts of a scalar  (serial version)
! Note: zeroing the unused ghosts, using d/dn=0 on domain bounds
!       'all'   =T means set all ghosts, not only 1 (default F)
!
  subroutine mympi_scalar(p,all)
    use global,  only: ndims,ni,nj,nk
    use utility, only: pbound
    implicit none
    real(8),dimension(ni,nj,nk),intent(inout) :: p
    logical,intent(in),optional               :: all
    real(8),dimension(:,:,:),pointer          :: put,get
    integer :: d
!
! -- loop boundaries
    do d=-ndims,ndims
       if(d.eq.0) cycle
!
! -- zero all ghosts
       get => pbound(p,d,all=.true.)
       get = 0
!
! -- point to desired ghosts
       get    => pbound(p, d,ghost=.true. ,all=all)
       if(me(d).eq.0) then
          put => pbound(p,-d,ghost=.false.,all=all)
       else
          put => pbound(p, d,ghost=.false.,all=all)
       end if
!
! -- set boundary value
       get = put
    end do

  end subroutine mympi_scalar
!
!----------------------------------------------------------------
!
! -- Fill ghosts of a vector (serial version)
! Note: zeroing the unused ghosts, using d/dn=0 on domain bounds
!       'all'=T means set all ghosts, not only 1 (default F)
!
  subroutine mympi_vector(u,all)
    use global,  only: ndims,ni,nj,nk
    use utility, only: pbound
    implicit none
    real(8),dimension(ndims,ni,nj,nk),intent(inout),target :: u
    logical,intent(in),optional                            :: all
    real(8),dimension(:,:,:),pointer                       :: v,put,get
    integer :: d,d2
    logical :: vector
!
! -- loop components and boundaries
    do d2=1,ndims
       v => u(d2,:,:,:)
       do d=-ndims,ndims
          if(d.eq.0) cycle
!
! -- flag to preserve values on domain faces
          vector = me(d).lt.0.and.abs(d).eq.d2
!
! -- zero all ghosts
          get => pbound(   v, d,ghost=.true. ,all=.true.,vector=vector)
          get = 0
!
! -- point to desired ghosts
          get => pbound(   v, d,ghost=.true. ,all=all)
          if(me(d).eq.0) then
             put => pbound(v,-d,ghost=.false.,all=all)
          else
             put => pbound(v, d,ghost=.false.,all=all,vector=vector)
          end if
          get = put
       end do
    end do

  end subroutine mympi_vector
!
!----------------------------------------------------------------
!
  subroutine mympi_end
    implicit none
    if(myid.eq.0) print *,'clean exit'
    return
  end subroutine mympi_end
!
!----------------------------------------------------------------
!
! -- Create datatypes for the full domain and subdomain arrays
!
  subroutine mympi_types(gisi,giei,gjsi,gjei,gksi,gkei)
  	use global, only: ndims
  	implicit none
  	integer,intent(in) :: gisi,giei,gjsi,gjei,gksi,gkei
  	gis=gisi; gie=giei; gjs=gjsi; gje=gjei; gks=gksi; gke=gkei
    allocate(sub_points(ndims))
    sub_points(1) = gie-gis+1
    sub_points(2) = gje-gjs+1
    if(ndims.eq.3) sub_points(3) = gke-gks+1    
  end subroutine mympi_types
!
!----------------------------------------------------------------
!
  subroutine mympi_write_scalar(file,array,lowres)
    implicit none
    integer,intent(in) :: file
    real(8),intent(in) :: array(:,:,:)
    logical,intent(in),optional :: lowres
!
! -- write to file
    if(present(lowres)) then
       if(lowres)      write(file,1) array
       if(.not.lowres) write(file,2) array
    else
       write(file) array
    end if

1   format(15f6.2)
2   format(15e12.4)
  end subroutine mympi_write_scalar
!
  subroutine mympi_write_vector(file,array)
    use global, only: ndims
    implicit none
    integer,intent(in) :: file
    real(8),intent(in) :: array(:,:,:,:)
    integer :: d
    do d=1,ndims
       call mympi_write_scalar(file,array(d,:,:,:))
    end do
  end subroutine mympi_write_vector
!
  subroutine mympi_read_scalar(file,array,bounds)
    implicit none
    integer,intent(in)          :: file
    real(8),intent(out)         :: array(:,:,:)
    logical,intent(in),optional :: bounds
!
! -- root reads from file
    read(file) array
    if(present(bounds)) then
       if(.not.bounds) return
    end if
    call mympi_scalar( array, all=.true. )
  end subroutine mympi_read_scalar
!
  subroutine mympi_read_vector(file,array)
    use global, only: ndims
    implicit none
    integer,intent(in) :: file
    real(8),intent(out) :: array(:,:,:,:)
    integer :: d
    do d=1,ndims
       call mympi_read_scalar(file,array(d,:,:,:),bounds=.false.)
    end do
    call mympi_vector( array, all=.true. )
  end subroutine mympi_read_vector
!
  subroutine mympi_write_subarray(file,array)
    implicit none
    integer,intent(in)   :: file
    real(8),intent(in)   :: array(:,:,:)
    write(file) sub_points
    write(file) array(gis:gie,gjs:gje,gks:gke)
    close(file)
  end subroutine mympi_write_subarray
!
!----------------------------------------------------------------
end module mympi
