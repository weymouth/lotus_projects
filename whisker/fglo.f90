!------------------------------------------------------------!
!--------- Global Variables and Utility Modules -------------!
!------------------------------------------------------------!
module global
  integer :: ndims,tmx,tmd,thk,ni,nj,nk,is,js,ks,ie,je,ke
  logical :: lres
  real(8) :: dt,dti,time0
  integer, allocatable :: points(:)
  real(8),parameter :: pi = 3.14159265359, zero = 0., one = 1.
!
  integer :: input_num = 5
  type print_flags
     logical :: prnt,rwnd,lwrs
     integer :: skip,tmod,file,ghst
  end type print_flags
!
end module global
module utility
  integer,private :: log_num = 9,indent = 0, cutoff = 60
  character(LEN=60),private,parameter :: space = &
       '                                                            '
  character(LEN=60) :: log=space
  interface get
     module procedure get,sget,iget
  end interface
contains
!
! -------------------------------------------------------
!
! -- Find the median of three values
!
  pure elemental subroutine median(fun,a,b,c)
    implicit none
    real(8),intent(in) :: a,b,c
    real(8),intent(out) :: fun
    fun=max(min(a,b),min(max(a,b),c))
  end subroutine median
!
! -------------------------------------------------------
!
! -- Quick utility for testing optionals
!
  pure elemental logical function optest(default,test)
    implicit none
    logical,intent(in) :: default
    logical,intent(in),optional :: test
    optest = default
    if(present(test)) then
       optest = test
    end if
  end function optest
!
! -------------------------------------------------------
!
! -- Array shifting subroutine
! Note: array has 'i' zeroed collums
!
  pure subroutine shift(i,j,a,b)
    use global, only : ndims,ni,nj,nk
    implicit none
    integer,intent(in) :: i,j
    real(8),intent(in),dimension(ni,nj,nk) :: a
    real(8),intent(out),dimension(ni,nj,nk) :: b
    integer :: ip,in
    
    if(j.gt.ndims) return
    ip = max(i,0); in = min(i,0) ! pos,neg

    if(j.eq.1) then
       b(1:1-in,:,:) = 0
       b(ni-ip:ni,:,:) = 0
       b(1-in:ni-ip,:,:) = a(1+ip:ni+in,:,:)
    else if(j.eq.2) then
       b(:,1:1-in,:) = 0
       b(:,nj-ip:nj,:) = 0
       b(:,1-in:nj-ip,:) = a(:,1+ip:nj+in,:)
    else
       b(:,:,1:1-in) = 0
       b(:,:,nk-ip:nk) = 0
       b(:,:,1-in:nk-ip) = a(:,:,1+ip:nk+in)
    end if

  end subroutine shift
!
! -------------------------------------------------------
!
! -- Shifted pointer function
! Note: dimensions reduced to exclude ghosts from vector field
!
  function pshift(a,shft,d) result(b)
    use global, only : ndims,ni,nj,nk,thk,is,ie,js,je,ks,ke
    implicit none
    real(8),dimension(:,:,:),pointer              :: b
    real(8),intent(in),dimension(ni,nj,nk),target :: a
    integer,intent(in)                            :: shft,d

    if(d.gt.ndims.or.d.lt.1) call nrerror('pshift d')
    if(thk+shft.lt.0.or.thk-shft.lt.1) call nrerror('pshift shft')
    if(d.eq.1) then
       b => a(is+shft:ie+1+shft,js:je,ks:ke)
    else if(d.eq.2) then
       b => a(is:ie,js+shft:je+1+shft,ks:ke)
    else
       b => a(is:ie,js:je,ks+shft:ke+1+shft)
    end if

  end function pshift
!
! -------------------------------------------------------
!
! -- Boundary pointer array
! Note: 'ghost' =T means get ghosts, not internals (default T)
!       'all'   =T means get all ghosts, not 1     (default F)
!       'vector'=T preserves domain boundary vals  (default F)
!
  function pbound(a,d,ghost,all,vector) result(b)
    use global, only: ndims,thk
    implicit none
    real(8),dimension(:,:,:),pointer           :: b
    real(8),dimension(:,:,:),target,intent(in) :: a
    integer,intent(in)                         :: d
    logical,intent(in),optional                :: ghost,all,vector
    integer :: first,last,skip,wdth,walk,put,vec

    if(abs(d).gt.ndims.or.abs(d).lt.1) call nrerror('pbound d')

    walk = 1; put = 0; wdth = 0; vec = 0
    if(present(ghost))  then
       if(.not.ghost) walk = -1
       if(.not.ghost) put  =  1
    end if
    if(present(all))  then
       if(all) wdth = thk-1
    end if
    if(present(vector)) then
       if(vector) vec = 1
    end if

    if(d.lt.0) then
       first = thk+put+put*vec
       last = first-walk*wdth
       skip = -walk
    else
       first = size(a,d)+1-thk-put+vec
       last = min(first+walk*wdth,size(a,d))
       skip = walk
    end if

    if(abs(d).eq.1) then
       b => a(first:last:skip,:,:)
    else if(abs(d).eq.2) then
       b => a(:,first:last:skip,:)
    else if(abs(d).eq.3) then
       b => a(:,:,first:last:skip)
    end if

  end function pbound
!
! -------------------------------------------------------
!
! -- Index shift
!
  function ishift(point,shift,d)
    use global, only: ndims,ni,nj,nk
    implicit none
    integer,intent(in) :: point(4),shift,d
    integer :: ishift(4)

    if(d.lt.0.or.d.gt.ndims) call nrerror('ishift error 1')
    ishift = point
    ishift(d+1) = ishift(d+1)+shift
    if(d.eq.0) ishift(1) = shift
    if(any(ishift.gt.(/ndims,ni,nj,nk/).or. &
           ishift.lt.(/1,1,1,1/))) then
       print '(10i4)',point,shift,d,ishift
       call nrerror('ishift error 2')
    end if

  end function ishift
!
! -- Get value from index
!
  real(8) function get(u,point)
    implicit none
    real(8),dimension(:,:,:,:) :: u
    integer,intent(in) :: point(4)
    get = u(point(1),point(2),point(3),point(4))
  end function get
  real(8) function sget(p,point)
    implicit none
    real(8),dimension(:,:,:) :: p
    integer,intent(in) :: point(4)
    sget = p(point(2),point(3),point(4))
  end function sget
  integer function iget(p,point)
    implicit none
    integer,dimension(:,:,:) :: p
    integer,intent(in) :: point(4)
    iget = p(point(2),point(3),point(4))
  end function iget
!
!-------------------------------------------------------
!
! -- Log file printing
!
  subroutine log_print(file_num)
    implicit none
    integer,optional,intent(in) :: file_num
    integer :: i
    i = len_trim(log)
    if(present(file_num)) then
       write(file_num,*) log(:i)
    else if(indent.lt.cutoff) then
       write(log_num,*) space(:indent),log(:i)
    end if
    log = space
  end subroutine log_print
!
  subroutine log_file_num(set_num)
    implicit none
    integer :: set_num
    log_num = set_num
  end subroutine log_file_num
  subroutine log_indent(adjust_indent)
    implicit none
    integer :: adjust_indent
    indent = indent+adjust_indent
  end subroutine log_indent
  subroutine log_cutoff(set_cutoff)
    implicit none
    integer :: set_cutoff
    cutoff = set_cutoff
  end subroutine log_cutoff
!
! -------------------------------------------------------
!
! -- Basic error checking and handling
!
!
  integer function assert_eq(nn,string)
    implicit none
    integer,dimension(:),intent(in) :: nn
    character(LEN=*),intent(in)     :: string
!
! -- report and die if integers are not all identical
    if(all(nn(2:)==nn(1))) then
       assert_eq = nn(1)
    else
       write(log_num,*) 'nrerror: assert_eq failed at:',&
            string
       stop 'program terminated by assert_eq'
    end if
  end function assert_eq
!
  subroutine nrerror(string)
    implicit none
    character(LEN=*),intent(in) :: string
!
! -- report and die
    write(log_num,*) 'nrerror:',string
    stop 'program terminated by nrerror'
  end subroutine nrerror
!
  subroutine rerror(string)
    implicit none
    character(LEN=*),intent(in) :: string
!
! -- report
    write(log_num,*) 'rerror:',string
  end subroutine rerror
!
end module utility
!
!!$program test
!!$  use global
!!$  use utility
!!$  use profile
!!$  implicit none
!!$  real(8),allocatable :: v(:,:,:,:)
!!$  real(8),pointer     :: p(:,:,:)
!!$  integer             :: d,i,j
!!$  logical             :: T=.true.,F=.false.,ghost,domain
!!$
!!$  ndims = 2; ni = 68; nj = 68; nk = 1; thk = 2
!!$  allocate(v(ndims,ni,nj,nk))
!!$
!!$  do i=1,2
!!$  do j=1,2
!!$     ghost  = merge(T,F,j.eq.1)
!!$     domain = merge(T,F,i.eq.2)
!!$     v = 0
!!$     do d=-ndims,ndims
!!$        if(d.eq.0) cycle
!!$        p => pbound(v(abs(d),:,:,:),d,ghost=ghost,domain=domain)
!!$        p = 1
!!$     end do
!!$     print '("ghost=",l," domain=",l)',ghost,domain
!!$     print 1,v(1,1:4,nj/2,1)
!!$     print 1,v(1,ni:ni-3:-1,nj/2,1)
!!$     print 1,v(2,ni/2,1:4,1)
!!$     print 1,v(2,ni/2,nj:nj-3:-1,1)
!!$     v = 0
!!$     do d=-ndims,ndims
!!$        if(d.eq.0) cycle
!!$        p => pbound(v(abs(d),:,:,:),d,ghost=ghost,vector=domain)
!!$        p = 1
!!$     end do
!!$     print '("ghost=",l," vector=",l)',ghost,domain
!!$     print 1,v(1,1:4,nj/2,1)
!!$     print 1,v(1,ni:ni-3:-1,nj/2,1)
!!$     print 1,v(2,ni/2,1:4,1)
!!$     print 1,v(2,ni/2,nj:nj-3:-1,1)
!!$  end do
!!$  end do
!!$  
!!$1 format(4e12.4)
!!$
!!$end program test
