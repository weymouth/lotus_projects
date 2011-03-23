!------------------------------------------------------------!
!                                                            !
!                 Blob Extraction Routines                   !
!                                                            !
!------------------------------------------------------------!
!                                                            !
module blob
  implicit none
  private
  public :: blob_find
!                                                            !
!------------------------------------------------------------!
!                                                            !
!
! -- Define types and operators
!
  type :: set                       ! tree structure
     private
     integer :: rank,parent
     logical :: seed
  end type set
  type :: data                      ! data structure
     real(8) :: vol,mom(3)
     integer :: imin(3),imax(3)
  end type data
contains
  subroutine MakeSet(S,x)           ! initialize tree as root
    type(set),intent(out) :: S(:)
    integer,intent(in)    :: x
    S(x)%parent = x
    S(x)%rank = 0
    S(x)%seed = .false.
  end subroutine MakeSet
  subroutine Union(S,x,y)           ! merge two trees
    type(set),intent(inout) :: S(:)
    integer,intent(in)      :: x,y 
    integer :: xRoot,yRoot         
!
! the tree with higher rank is made root to speed-up Find.
! two trees with equal rank (that aren't already merged)
! are merged arbitrarily and their rank is incremented.
!
    xRoot = Find(S,x)               
    yRoot = Find(S,y)              
    if((S(xRoot)%seed.and.S(yRoot)%seed).or.xRoot.eq.yRoot) then
       return
    else if(S(xRoot)%rank.gt.S(yRoot)%rank.or.S(xRoot)%seed) then 
       S(yRoot)%parent = xRoot
    else if(S(xRoot)%rank.lt.S(yRoot)%rank.or.S(yRoot)%seed) then 
       S(xRoot)%parent = yRoot
    else
       S(yRoot)%parent = xRoot
       S(xRoot)%rank   = S(xRoot)%rank+1
    end if
  end subroutine Union
  integer recursive function Find(S,x)result(y) ! find root
    type(set),intent(inout) :: S(:)
    integer,intent(in)      :: x
    y = S(x)%parent
    if(S(y)%parent.ne.y) then ! if parent isn't a root,
       S(x)%parent = Find(S,y)   ! link directly to parent's root,
       y = S(x)%parent           ! and return it.
    end if
  end function Find
  recursive subroutine PrintSet(S,x) ! print out tree
    type(set),intent(in) :: S(:)
    integer,intent(in)   :: x
    integer :: y
    y = S(x)%parent
    print '("label =",i3,", rank =",i3,", parent =",i3)',x,S(x)%rank,y
    if(x.ne.y) then
       call PrintSet(S,y)
    else
       print '("-------- root ---------")'
    end if
  end subroutine PrintSet
  elemental subroutine InitData(d)   ! init data elements
    type(data),intent(out) :: d
    d%vol  = 0
    d%mom  = 0
    d%imin = 1e6
    d%imax = -1e6
  end subroutine InitData
  subroutine Add(d,vol,cen)          ! update blob data with point info
    type(data),intent(inout) :: d
    real(8),intent(in)       :: vol    ! point 'volume'
    integer,intent(in)       :: cen(3) ! point 'centroid'
    d%vol    = d%vol+vol
    d%mom    = d%mom+vol*dble(cen)     ! moment=centriod*volume
    d%imin   = min(d%imin,cen)
    d%imax   = max(d%imax,cen)
  end subroutine Add
!                                                                !
!----------------------------------------------------------------!
!                                                                !
! -- Find connected "blobs" of nonzero points in the input array
!
! Blobs are found by scanning the array twice. First, nonzero 
!  points are found and given the label of a connected existing
!  blob if present and given a new label otherwise. Additionally,
!  connections between the blobs themselves are tracked using the
!  fast Union/Find algorithms above. In the second pass, the blob
!  fragments are merged and relabeled sequentially.
! This routine can also print lagrangian data for each blob,
!  print a histogram of the blob lengths, return the blob label
!  array, and remove blobs from the input array whose length is 
!  less than 'cutoff'.
! The routine is serial with no dependencies.
!
  subroutine blob_find(r,cutoff,d_num,h_num,R_num,thresh,time,iout,t_num)
    real(8),intent(inout)        :: r(:,:,:)
    integer,intent(in),optional  :: cutoff,d_num,h_num,R_num,t_num
    real(8),intent(in),optional  :: time,thresh(:)
    integer,intent(out),optional :: iout(:,:,:)
    integer,allocatable    :: list(:),a(:,:,:)
    type(set),allocatable  :: blob(:)
    type(data),allocatable :: blob_data(:)
    type(data) :: b
    integer    :: ni,nj,nk,n,i,j,k,d,p(3),p2(3),label,local,tmax,t
    real(8)    :: tval,sumvol,sumf
!
! -- init
    ni = size(r,1); nj = size(r,2); nk = size(r,3)
    allocate(blob(size(r)/2),a(ni,nj,nk))
    a = 0; n = 0
!
! -- loop over threshold levels
    tmax = 1; tval = 0.
    if(present(thresh)) tmax = size(thresh)
    do t=1,tmax
       if(present(thresh)) tval = thresh(t)
!
! -- find points above the threshold
    do k=1,nk
    do j=1,nj
    do i=1,ni
       if(r(i,j,k).le.tval) cycle
       p = (/i,j,k/)
!
! -- check for blob fragments and link them
       label = a(i,j,k)
       do d=3,1,-1
          p2 = p; p2(d) = p(d)-1         ! look in direction d
          if(any(p2.le.0)) cycle         ! loop if out
          local = a(p2(1),p2(2),p2(3))   ! get local value
          if(local.ne.0) then            ! local tree exists
             if(label.ne.0) &               ! labeled tree exists
                call Union(blob,local,label)   ! merge trees (if not seeds)
             if(label.eq.0) label = local   ! label tree
          end if
       end do
!
! -- if none, start new blob
       if(label.eq.0) then
          n = n+1
          label = n
          call MakeSet(blob,n)
       end if
!
! -- label point
       a(i,j,k) = label
    enddo
    enddo    
    enddo
!
! -- fix the roots at this level before looping
    do label = 1,n
       local = Find(blob,label)
       blob(local)%seed = .true.
    enddo
    enddo
    if(n.eq.0) return
!
! -- relabel blobs sequentially and gather data
!
    allocate(list(n),blob_data(n))
    call InitData(blob_data); list = 0; n = 0
    do k=1,nk
    do j=1,nj
    do i=1,ni
       if(a(i,j,k).eq.0) cycle
       a(i,j,k) = Find(blob,a(i,j,k))  ! get tree root
       label = list(a(i,j,k))          ! get root label
       if(label.eq.0) then             ! new root
          n = n+1; label = n              ! get label
          list(a(i,j,k)) = n              ! label root
       end if
       a(i,j,k) = label
!
! note: all of the blob data are in numerical units
! note:  to get physical units, scale by the grid metrics
!        call Add(blob_data(label),r(i,j,k)*dv,(/i,j,k/)*dx+offset)
!
       call Add(blob_data(label),r(i,j,k),(/i,j,k/))
    end do
    end do
    end do
!
! -- print blob data and save approximate lengths
!
! note: box size+1> exact length > box size-1
! note:  where, box size = imax-imin
!
    list = -1
    if(present(d_num).and.present(time)) write(d_num,1) time
    do label=1,n
       b = blob_data(label)
       if(present(d_num)) write(d_num,'(i6,4e12.4,3i6)') &
            label,b%vol,b%mom/b%vol,b%imax-b%imin
       if(present(R_num).and.present(time)) write(R_num, &
            '(e12.4,",",e12.4,",",e12.4,",",e12.4,",",e12.4,",",i6,",",i6,",",i6)') &
            time,b%vol,b%mom/b%vol,b%imax-b%imin
       list(label) = maxval(b%imax-b%imin) ! 'length'=max(box size)
    end do
!
! -- print total volume
    if(present(t_num).and.present(time)) then
	sumvol=0.D0
	sumf=sum(r(1:ni,1:nj,1:nk))
	do label=1,n
	b = blob_data(label)
	sumvol=sumvol+b%vol
	end do
	write(t_num,'(e12.4,2e18.8)') time,sumvol,sumf
    endif
!
! -- print histogram of lengths
    if(present(h_num).and.present(time)) write(h_num,1) time
    do label=0,maxval(list)
       n = count(list.eq.label)
       if(n.eq.0) cycle
       if(present(h_num)) write(h_num,'(2i6)') label,n
    end do
!
! -- remove blobs with length<cutoff
!
! note: set cutoff.le.ghost thickness: erroneous removals
! note:  can be corrected by updating the boundary values
!
    if(present(cutoff).and.cutoff.gt.0) then
       forall(i=1:ni,j=1:nj,k=1:nk,list(a(i,j,k)).lt.cutoff)
          r(i,j,k) = 0
          a(i,j,k) = -a(i,j,k)   ! mark as removed
       end forall
    end if
!
! -- output blob label array
    if(present(iout)) iout = a
    return

1   format('ZONE SOLUTIONTIME = ',e12.4,', F=POINT')
  end subroutine blob_find
end module blob
!!$program test
!!$  use blob
!!$  implicit none
!!$!  integer,parameter :: ni=5,nj=8
!!$  integer,parameter :: ni=22,nj=18
!!$  character(7)  :: fmt
!!$  character(3)  :: c(ni,nj,1)
!!$  real(8)       :: r(ni,nj,1)
!!$  integer       :: a(ni,nj,1)
!!$
!!$  fmt = '(   a3)'
!!$  write(fmt(2:4),'(i3)') ni
!!$  call random_seed
!!$  call random_number(r)
!!$!  r(:,1,1) = (/0,1,0,0,0/)
!!$!  r(:,2,1) = (/1,1,0,0,1/)
!!$!  r(:,3,1) = (/0,1,0,1,1/)
!!$!  r(:,4,1) = (/1,1,0,0,1/)
!!$!  r(:,5,1) = (/0,1,0,1,1/)
!!$!  r(:,6,1) = (/1,1,0,0,1/)
!!$!  r(:,7,1) = (/0,1,0,1,1/)
!!$!  r(:,8,1) = (/1,1,1,1,1/)
!!$  a = merge(1,0,r.ge.0.5)
!!$
!!$  print '("--------- input --------")'
!!$  where(r.lt.0.5) r = 0
!!$  write(c,'(f3.1)') r
!!$  where(c.eq.'0.0') c = '   '
!!$  where(c.eq.'1.0') c = '  1'
!!$  where(c(:,:,:)(1:1).eq.'0') c(:,:,:)(1:1) = ' '
!!$  write(*,fmt) c
!!$
!!$  call blob_find(r,cutoff=2,d_num=6,h_num=6,time=0.,iout=a)
!!$  print '("--------- labeled --------")'
!!$  write(c,'(i3)') a
!!$  where(c.eq.'  0') c = '   '
!!$  write(*,fmt) c
!!$
!!$end program test
!!$program test
!!$  use blob
!!$  use profile
!!$  implicit none
!!$  integer,parameter :: ni=100,nj=100,nk=100,N=10
!!$  real(8)       :: r(ni,nj,nk),c
!!$  integer       :: a(ni,nj,nk),i
!!$
!!$  call init_profile
!!$
!!$  call random_seed
!!$  call random_number(r)
!!$  call profile_time('rand ')
!!$
!!$  do i=1,N
!!$     c = dble(i)/dble(N+1)
!!$     where(r.lt.c) r = 0
!!$     call blob_find(r,h_num=8,time=c)
!!$  end do
!!$  call profile_time('blob ')
!!$  call profile_report(c)
!!$
!!$end program test
