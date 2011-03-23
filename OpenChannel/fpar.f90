!------------------------------------------------------------!
!------------------- MPI Data and Routines ------------------!
!------------------------------------------------------------!
module mympi
  use utility, only: nrerror
  include 'mpif.h'
  integer,private,allocatable,dimension(:) :: coords,me,blocks
  integer,private                          :: myid,myidmx,ios,ierr
  integer,private                          :: mpi_grid,local_type,local_subtype
  integer,private,allocatable,dimension(:) :: mpi_dom_grid,global_type
  integer,private,allocatable,dimension(:) :: global_subtype,sub_points
  logical,private                          :: subflag
  logical,private,allocatable,dimension(:) :: global_subflag,periodic
  real(8),private                          :: mpi_time
contains
!----------------------------------------------------------------
!
! -- Initialize a parrallel run, reading in inputs, computing dims
! setting up the topography
!
  subroutine init_mympi
    use global
    use utility, only: nrerror,log,log_print
    implicit none
    integer :: d
    logical :: reorder=.true.
!
! -- Set up default MPI communicator
!
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    call MPI_COMM_DUP (MPI_COMM_WORLD, mpi_grid, ierr)    
    mpi_time = MPI_WTIME()
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
    allocate(blocks(0:ndims-1),points(ndims),periodic(0:ndims-1))
    call mympi_read(input_num)
    call mympi_read(input_num,la=periodic)
    call mympi_read(input_num,ia=blocks)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, d, ierr)
    if(product(blocks).ne.d) call nrerror('init_mympi: wrong number of processors')
    myidmx = d-1
    write(log,'("  processors=",i5)') d
    call log_print
!
! -- Create parrallel topology
!
    allocate(me(-ndims:ndims),coords(0:ndims-1))
    call MPI_CART_CREATE(MPI_COMM_WORLD, ndims, blocks, &
                         periodic, reorder, mpi_grid, ierr)
    call MPI_COMM_RANK  (mpi_grid, myid, ierr)               ! redefine myid
    call MPI_CART_COORDS(mpi_grid, myid, ndims, coords, ierr)
    do d=1,ndims
       call MPI_CART_SHIFT(mpi_grid, d-1, 1, me(-d), me(d), ierr)
    end do
    if(ierr.ne.0) call nrerror('init_mympi: mpi_grid')
    me(0) = myid
!
! -- Make domain boundary communicators
!    
    allocate(mpi_dom_grid(ndims))
    do d=1,ndims
       periodic = .true.; periodic(d-1) = .false. ! periodic used as 'remain'
       call MPI_Cart_sub(mpi_grid,periodic,mpi_dom_grid(d),ierr)       
    end do
    if(ierr.ne.0) call nrerror('init_mympi: mpi_dom_grid')
  end subroutine init_mympi
!
! -- Basic query functions
!
  logical function mympi_domain_bound(d)
    integer,intent(in) :: d
    mympi_domain_bound = me(d).lt.0
  end function mympi_domain_bound
!
  integer function mympi_blocks(d)
    integer,intent(in) :: d
    mympi_blocks = blocks(d-1)
  end function mympi_blocks
!
  integer function mympi_coords(d)
    integer,intent(in) :: d
    mympi_coords = coords(d-1)
  end function mympi_coords
!
  pure logical function mympi_periodic(d)
    integer,intent(in) :: d
    mympi_periodic = periodic(d-1)
  end function mympi_periodic
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
  integer function mympi_idmx()
    implicit none
    mympi_idmx = myidmx
  end function mympi_idmx
!
  integer function mympi_id()
    implicit none
    mympi_id = myid
  end function mympi_id
!
  integer function mympi_comm()
    mympi_comm = mpi_grid
  end function mympi_comm
!
!----------------------------------------------------------------
!
! -- Read and distribute file info
!
  subroutine mympi_read(num,l,la,i,ia,r,ra,ru,f)
    use global,  only: print_flags
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
!
! -- distribute
    if(present(l)) then
       call MPI_BCAST(l,1,MPI_LOGICAL,0,mpi_grid,ierr)
    else if(present(la)) then
       call MPI_BCAST(la,size(la),MPI_LOGICAL,0,mpi_grid,ierr)
    else if(present(i)) then
       call MPI_BCAST(i,1,MPI_INTEGER,0,mpi_grid,ierr)
    else if(present(ia)) then
       call MPI_BCAST(ia,size(ia),MPI_INTEGER,0,mpi_grid,ierr)
    else if(present(r)) then
       call MPI_BCAST(r,1,MPI_DOUBLE_PRECISION,0,mpi_grid,ierr)
    else if(present(ra)) then
       call MPI_BCAST(ra,size(ra),MPI_DOUBLE_PRECISION,0,mpi_grid,ierr)
    else if(present(ru)) then
       call MPI_BCAST(ru,1,MPI_DOUBLE_PRECISION,0,mpi_grid,ierr)
    else if(present(f)) then
       call MPI_BCAST(f%prnt,1,MPI_LOGICAL,0,mpi_grid,ierr)
       call MPI_BCAST(f%lwrs,1,MPI_LOGICAL,0,mpi_grid,ierr)
       call MPI_BCAST(f%rwnd,1,MPI_LOGICAL,0,mpi_grid,ierr)
       call MPI_BCAST(f%file,1,MPI_INTEGER,0,mpi_grid,ierr)
       call MPI_BCAST(f%skip,1,MPI_INTEGER,0,mpi_grid,ierr)
       call MPI_BCAST(f%tmod,1,MPI_INTEGER,0,mpi_grid,ierr)
       call MPI_BCAST(f%ghst,1,MPI_INTEGER,0,mpi_grid,ierr)
    end if
    if(ierr.ne.0) call nrerror('mympi_read: MPI')

  end subroutine mympi_read
!
!----------------------------------------------------------------
!
! -- Get the elapsed time
!
  subroutine mympi_timer(elapsed)
    implicit none
    real(8),intent(out) :: elapsed
    real(8) :: time

    call MPI_BARRIER(mpi_grid,ierr)
    time = MPI_WTIME()
    elapsed = time-mpi_time
    mpi_time = time
    if(ierr.ne.0) call nrerror('mympi_timer')

  end subroutine mympi_timer
!
!----------------------------------------------------------------
!
! -- Sum over proccessors
!
  subroutine mympi_sum(each,all)
    implicit none
    real(8),intent(in)          :: each
    real(8),intent(out)         :: all
    call MPI_ALLREDUCE(each, all, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, mpi_grid, ierr)
    if(ierr.ne.0) call nrerror('mympi_sum')
  end subroutine mympi_sum
!
  subroutine mympi_domain_sum(d,each,all)
    implicit none
    integer,intent(in)          :: d
    real(8),intent(in)          :: each
    real(8),intent(out)         :: all
    call MPI_ALLREDUCE(each, all, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, mpi_dom_grid(d), ierr)
    if(ierr.ne.0) call nrerror('mympi_domain_sum')
  end subroutine mympi_domain_sum
!
!----------------------------------------------------------------
!
! -- Compare across proccessors
!
  logical function mympi_any(each)
    implicit none
    logical :: each
    call MPI_ALLREDUCE(each, mympi_any, 1, &
         MPI_LOGICAL, MPI_LOR, mpi_grid, ierr)
    if(ierr.ne.0) call nrerror('mympi_any')
  end function mympi_any
!
  logical function mympi_all(each)
    implicit none
    logical :: each
    call MPI_ALLREDUCE(each, mympi_all, 1, &
         MPI_LOGICAL, MPI_LAND, mpi_grid, ierr)
    if(ierr.ne.0) call nrerror('mympi_all')
  end function mympi_all
!
  real(8) function mympi_max(each)
    implicit none
    real(8) :: each
    call MPI_ALLREDUCE(each, mympi_max, 1, &
         MPI_DOUBLE_PRECISION, MPI_MAX, mpi_grid, ierr)
    if(ierr.ne.0) call nrerror('mympi_max')
  end function mympi_max
!
!----------------------------------------------------------------
!
! -- Fill ghosts of a scalar
! Note: zeroing the unused ghosts, using d/dn=0 on domain bounds
!       'all'   =T means set all ghosts, not only 1 (default F)
!
  subroutine mympi_scalar(p,all)
    use global,  only: ndims,ni,nj,nk
    use utility, only: pbound,nrerror
    implicit none
    real(8),dimension(ni,nj,nk),intent(inout) :: p
    logical,intent(in),optional               :: all
    real(8),dimension(:,:,:),pointer          :: put,get
    integer :: d,num,status(MPI_STATUS_SIZE)
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
       get => pbound(p, d,ghost=.true. ,all=all)
       put => pbound(p,-d,ghost=.false.,all=all)
!
! -- communicate with neighbor
       num = size(put)
       CALL MPI_SENDRECV(&
            put, num, MPI_DOUBLE_PRECISION, me(-d), d+ndims,  &
            get, num, MPI_DOUBLE_PRECISION, me( d), d+ndims,  &
            mpi_grid, status, ierr )
       if(ierr.ne.0) call nrerror('mympi_scalar')
!
! -- if domain bound, get info from this boundary
       if(me(d).lt.0) then
          put => pbound(p,d,ghost=.false.,all=all)
          get = put
       end if
    end do

  end subroutine mympi_scalar
!
!----------------------------------------------------------------
!
! -- Fill ghosts of a vector
! Note: zeroing the unused ghosts
!       'all'=T means set all ghosts, not only 1 (default F)
!
  subroutine mympi_vector(u,all)
    use global,  only: ndims,ni,nj,nk
    use utility, only: pbound,nrerror
    implicit none
    real(8),dimension(ndims,ni,nj,nk),intent(inout),target :: u
    logical,intent(in),optional                            :: all
    real(8),dimension(:,:,:),pointer                       :: v,put,get
    integer :: d,d2,num,status(MPI_STATUS_SIZE)
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
          get => pbound(v, d,ghost=.true. ,all=.true.,vector=vector)
          get = 0
!
! -- point to desired ghosts
          get => pbound(v, d,ghost=.true. ,all=all)
          put => pbound(v,-d,ghost=.false.,all=all)
!
! -- communicate with neighbor
          num = size(put)
          CALL MPI_SENDRECV(&
               put, num, MPI_DOUBLE_PRECISION, me(-d), 6*d2+d,  &
               get, num, MPI_DOUBLE_PRECISION, me( d), 6*d2+d,  &
               mpi_grid, status, ierr )
          if(ierr.ne.0) call nrerror('mympi_vector')
!
! -- if domain boundary, use local info for ghosts
          if(me(d).lt.0) then
             put => pbound(v,d,ghost=.false.,all=all,vector=vector)
             get = put
          end if
       end do
    end do

  end subroutine mympi_vector
!
!----------------------------------------------------------------
!
  subroutine mympi_end
    implicit none
    if(myid.eq.0) print *,'clean exit'
    call MPI_FINALIZE(ierr)
  end subroutine mympi_end
!
!----------------------------------------------------------------
!
! -- Create datatypes for the full domain and subdomain arrays
!
  subroutine mympi_types(gis,gie,gjs,gje,gks,gke)
    use global
    implicit none
    integer,intent(in)           :: gis,gie,gjs,gje,gks,gke
    integer                      :: d,mpi_subgrid,id
    integer,dimension(0:ndims-1) :: sizes,subsizes,starts
    integer,allocatable          :: temp(:,:)
    logical,dimension(0:ndims-1) :: remain
!
! -- local array type which includes domain ghosts
!
    sizes(0:1) = (/ni,nj/); subsizes(0:1) = (/ie-is+1,je-js+1/)
    if(ndims.eq.3) sizes(2) = nk
    if(ndims.eq.3) subsizes(2) = ke-ks+1
    starts = 0
    do d=1,ndims
       if(coords(d-1).eq.0) subsizes(d-1) = subsizes(d-1)+thk
       if(coords(d-1).ne.0) starts(d-1) = starts(d-1)+thk
       if(coords(d-1).eq.blocks(d-1)-1) &
            subsizes(d-1) = subsizes(d-1)+thk
    end do
    call MPI_Type_create_subarray(ndims,sizes,subsizes,starts,&
         MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, local_type, ierr)
    call MPI_Type_commit(local_type,ierr)
    if(ierr.ne.0) call nrerror('mympi_types: local_type')
!
! -- global domain array type
!
    allocate(global_type(0:myidmx))
    sizes = points+2*thk
!
! -- create an image of each array on root
    if(myid.eq.0) then
       do id=myidmx,0,-1
          call MPI_Cart_coords(mpi_grid,id,ndims,coords,ierr)
          subsizes(0:1) = (/ie-is+1,je-js+1/)
          if(ndims.eq.3) subsizes(2) = ke-ks+1
          starts = coords*subsizes
          do d=1,ndims
             if(coords(d-1).eq.0) subsizes(d-1) = subsizes(d-1)+thk
             if(coords(d-1).ne.0) starts(d-1) = starts(d-1)+thk
             if(coords(d-1).eq.blocks(d-1)-1) &
                                 subsizes(d-1) = subsizes(d-1)+thk
          end do
          call MPI_Type_create_subarray(ndims,sizes,subsizes,starts,&
               MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
               global_type(id), ierr)
          call MPI_Type_commit(global_type(id),ierr)
       end do
    end if
    if(ierr.ne.0) call nrerror('mympi_types: global_type')
!
! -- local subdomain array type using input positions
!
    sizes(0:1) = (/ni,nj/); starts(0:1) = (/gis-1,gjs-1/)
    subsizes(0:1) = (/max(gie-gis+1,0),max(gje-gjs+1,0)/)
    if(ndims.eq.3) then
       sizes(2) = nk; subsizes(2) = max(gke-gks+1,0); starts(2) = gks-1
    end if
    subflag = .not.any(subsizes.eq.0)
    if(subflag) then
       call MPI_Type_create_subarray(ndims,sizes,subsizes,starts,&
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, local_subtype, ierr)
       call MPI_Type_commit(local_subtype,ierr)
    end if
    if(ierr.ne.0) call nrerror('mympi_types: local_subtype')
!
! -- global subdomain array type
!
    allocate(temp(ndims,maxval(blocks)),sub_points(ndims),&
         global_subflag(0:myidmx),global_subtype(0:myidmx))
    do d=1,ndims
       remain = .false.; remain(d-1) = .true.
       call MPI_Cart_sub(mpi_grid,remain,mpi_subgrid,ierr)       
       call MPI_Allgather(subsizes(d-1),1,MPI_INTEGER,&
            temp(d,:),1,MPI_INTEGER,mpi_subgrid,ierr)
       sizes(d-1) = sum(temp(d,1:blocks(d-1)))
    end do
    sub_points = sizes(0:ndims-1)
    if(ierr.ne.0) call nrerror('mympi_types: mpi_subgrid')
!
! -- create an image of each subarray on root
    if(myid.eq.0) then
       do id=myidmx,0,-1
          call MPI_Cart_coords(mpi_grid,id,ndims,coords,ierr)
          do d=1,ndims
             subsizes(d-1) = temp(d,coords(d-1)+1)
             starts(d-1) = sum(temp(d,1:coords(d-1)))
          end do
          subflag = .not.any(subsizes.eq.0)
          global_subflag(id) = subflag
          if(subflag) then
             call MPI_Type_create_subarray(ndims,sizes,subsizes,starts,&
                  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
                  global_subtype(id), ierr)
             call MPI_Type_commit(global_subtype(id),ierr)
          end if
       end do
    end if
    if(ierr.ne.0) call nrerror('mympi_types: global_subtype')
  end subroutine mympi_types
!
!----------------------------------------------------------------
!
  subroutine mympi_write_subarray(file,array)
    implicit none
    integer,intent(in)   :: file
    real(8),intent(in)   :: array(:,:,:)
    integer :: status(MPI_STATUS_SIZE),request,id
    real(8) :: rarray(sub_points(1),sub_points(2),sub_points(3))
!
! -- everyone sends local_subtype
    if(subflag) call MPI_Isend( array, 1, local_subtype, &
         0, myid, MPI_COMM_WORLD, request, ierr)
!
! -- root receives each global_subtype
    if(myid.eq.0) then
       do id=0,myidmx
          if(global_subflag(id)) call MPI_Recv( rarray, 1, &
               global_subtype(id), id, id, &
               MPI_COMM_WORLD, status, ierr)
       end do
!
! -- write to file and close
!!$       write(file,*) sub_points
!!$       write(file,'(15e12.4)') rarray
       write(file) sub_points
       write(file) rarray
       close(file)
    end if

    if(ierr.ne.0) call nrerror('mympi_write_subarray')
  end subroutine mympi_write_subarray
!
!----------------------------------------------------------------
end module mympi
