!------------------------------------------------------------!
!------------------- I/O Data and Routines ------------------!
!------------------------------------------------------------!
module inout
!
! Note: Much of this is wrapper. The fact is that the low-level
!       routines need I/O but the I/0 requires low-level routines
!       and data. There isn't a clean place to cut.
!
!       My solution is to make this wrapper I/O module so that the
!       high-level modules are unaware of the tangled confusion.
!
  use global,  only: print_flags
  use utility, only: io_log=>log,log_print,io_error=>nrerror
  interface io_log_field
      module procedure io_log_scalar,io_log_vector
   end interface
   integer              :: force_num = 10, press_num = 8, CFL_num = 11
   integer,private      :: io_log_num = 9, profile_num = 0
   real(8),private      :: profile_dt(100)
   character(5),private :: profile_name(100)
contains
!
! -- Initialize low-level stuff
!
  subroutine init_io
    use global
    use mympi
    use utility, only: log_file_num
    use grid,    only: init_grid     !<--- Content here
    implicit none
    integer :: gis,gie,gjs,gje,gks,gke
    
    call log_file_num(io_log_num)
    write(io_log,'("Simulation Initiated")')
    call log_print
    call init_profile

    call init_mympi
    call init_grid(gis,gie,gjs,gje,gks,gke)
    call mympi_types(gis,gie,gjs,gje,gks,gke)

  end subroutine init_io
!
! -- Read and distribute file info
!
  subroutine io_read(file_num,l,la,i,ia,r,ra,ru,f)
    use global, only: input_num
    use mympi,  only: mympi_read     !<--- Content here
    implicit none
    integer,intent(in),optional            :: file_num
    logical,intent(inout),optional         :: l,la(:)
    integer,intent(inout),optional         :: i,ia(:)
    real(8),intent(inout),optional         :: r,ra(:),ru
    type(print_flags),intent(out),optional :: f
    integer          :: num

    num = input_num
    if(present(file_num)) then
       num = file_num
    end if
    call mympi_read(num,l,la,i,ia,r,ra,ru,f)

  end subroutine io_read
!
! -------------------------------------------------------
!
! -- Profile routines
!
  subroutine init_profile
    implicit none
    profile_dt = 0
    profile_name = '     '
  end subroutine init_profile
!
  subroutine io_profile(name)
    use mympi, only: mympi_timer
    implicit none
    character(5),intent(in) :: name
    logical :: debug=.false.,found
    real(8) :: time
    integer :: i

    call mympi_timer(time)
    found = .false.
    do i=1,profile_num                      ! check old profiles
       if(name.eq.profile_name(i)) then     ! if old, add dt
          profile_dt(i) = profile_dt(i)+time
          found = .true.
       end if
    end do
    if(.not.found) then                     ! set up new profile
       profile_num = profile_num+1
       profile_dt(profile_num) = time
       profile_name(profile_num) = name
    end if
    if(debug) then                          ! write to log file
       write(io_log,'("profile: ",a5)') name
       call log_print
       call flush(io_log_num)
    end if
    
  end subroutine io_profile
!
  subroutine io_profile_report(t)
    use mympi
    use global,  only: tmx
    implicit none 
    integer,intent(in) :: t
    real(8)            :: time
    logical            :: unsorted(profile_num)
    integer            :: i,j

    time = sum(profile_dt)
    if(t.eq.tmx) then                        ! say goodbye
       write(io_log,'("Time step =",i5,", Clean exit")') t
    else                                     ! estimate remaining time
       write(io_log,2) t,time*(tmx-t)/t
    end if
    call log_print
       
    if(time.gt.10) then
       write(io_log,'(" Profile Report. Total Time=",e12.4)') time
       call log_print
       unsorted = .true.
       sort_loop: do i=1,profile_num
          j = sum(maxloc(profile_dt(1:profile_num),unsorted))
          unsorted(j) = .false.
          if(profile_dt(j)/time.lt.0.005.or.profile_dt(j).lt.1) &
               exit sort_loop
          write(io_log,1) i,profile_name(j),100*profile_dt(j)/time &
               ,profile_dt(j)
          call log_print
       end do sort_loop
    end if
    call flush(io_log_num)
    call flush(force_num)
    call flush(press_num)
    call flush(CFL_num)

1   format(i3,':  ',a5,f9.4,'%,',e12.4)
2   format('Time step =',i5,', Estimated time remaining =',e12.4)
  end subroutine io_profile_report
!
! -- General print routine for scalar/vector fields
! 
  subroutine io_write(t,pflags,p,u,error)
    use global
    use mympi
    use grid,    only: grid_position_array
    use utility, only: shift
    implicit none
    integer,intent(in)                                    :: t
    type(print_flags),intent(in),optional                 :: pflags
    real(8),intent(in),optional,dimension(ni,nj,nk)       :: p
    real(8),intent(in),optional,dimension(ndims,ni,nj,nk) :: u
    logical,intent(in),optional                           :: error
    real(8),dimension(ni,nj,nk) :: v
    integer :: d,file,skip,tmod,ghst,iss,iee,jss,jee,kss,kee,sizes(ndims)
    logical :: rwnd,lwrs,scalar,debug=.false.
    real(8) :: time

    scalar = present(p)
    if(scalar.and.present(u)) call io_error('io_write: two fields')
    if(.not.scalar.and..not.present(u)) call io_error('io_write: no field')

    if(present(pflags)) then
       if(.not.pflags%prnt) return
       file = pflags%file
       skip = pflags%skip
       tmod = pflags%tmod
       rwnd = pflags%rwnd
       lwrs = pflags%lwrs
       ghst = pflags%ghst
    else if(present(error)) then
       file = 9000
       skip = 1
       tmod = 1
       rwnd = .true.
       lwrs = .false.
       ghst = 0
    else
       call io_error('io_write: no flags')
    end if
    sizes = (/ni,nj,nk/)
    if(debug) write(9,'("pflags=",i4,3i2,3l)') file,skip,tmod,ghst, &
         scalar,rwnd,lwrs
    if(mod(t,tmod).eq.0) then
       file = file+mympi_id()
       time = t*dt+time0
       if(rwnd.or.t.eq.0) rewind(file)
       if(ndims.eq.3) then
          if(scalar) then
             write(file,*)'VARIABLES=x,y,z,p'
          else
             write(file,*)'VARIABLES=x,y,z,u,v,w'
          end if
          if(rwnd) then
             write(file,5) sizes
          else
             write(file,1) time,sizes
          end if
       else
          if(scalar) then
             write(file,*)'VARIABLES=x,y,p'
          else
             write(file,*)'VARIABLES=x,y,u,v'
          end if
          if(rwnd) then
             write(file,6) sizes
          else
             write(file,2) time,sizes
          end if
       end if
       if(rwnd.or.t.eq.0) then
          do d=1,ndims
             v = grid_position_array(d,0)
             write(file,3) v
          end do
       else
          if(ndims.eq.2) write(file,*)'VARSHARELIST=([1-2]=1)'
          if(ndims.eq.3) write(file,*)'VARSHARELIST=([1-3]=1)'
       end if
       if(scalar) then
          if(lwrs) write(file,3) p
          if(.not.lwrs) write(file,4) p
       else
          do d=1,ndims
             call shift(1,d,u(d,:,:,:),v)
             v = (v+u(d,:,:,:))*0.5
             if(lwrs) write(file,3) v
             if(.not.lwrs) write(file,4) v
          end do
       end if
       call flush(file)
    end if

1   format('ZONE SOLUTIONTIME = ',e12.4,', I =',i5, &
         ', J = ',i5,', K = ',i5,', F=BLOCK')
2   format('ZONE SOLUTIONTIME = ',e12.4,', I =',i5, &
         ', J = ',i5,', F=BLOCK')
3   format(15f6.2)
4   format(15e12.4)
5   format('ZONE I =',i5,', J = ',i5,', K = ',i5,', F=BLOCK')
6   format('ZONE I =',i5,', J = ',i5,', F=BLOCK')
  end subroutine io_write
!
!----------------------------------------------------------------
!
! -- Log the name and L_1{s,v}
!
  subroutine io_log_vector(name,v,ghosts)
    use global
    use mympi
    implicit none
    character(5),intent(in)     :: name
    real(8),intent(in)          :: v(ndims,ni,nj,nk)
    logical,intent(in),optional :: ghosts
    real(8),dimension(ndims)    :: each,all
    logical                     :: flag = .false.
    integer                     :: d

    if(present(ghosts)) flag = .not.ghosts
    do d=1,ndims
       if(flag) then
          each(d) = sum(abs(v(d,is:ie,js:je,ks:ke)))
       else
          each(d) = sum(abs(v(d,:,:,:)))
       end if
       call mympi_sum(each(d),all(d))
    end do

    write(io_log,'(a5,3e16.8)') name,all
    call log_print

  end subroutine io_log_vector
!
  subroutine io_log_scalar(name,s,ghosts)
    use global
    use mympi
    implicit none
    character(5),intent(in)     :: name
    real(8),intent(in)          :: s(ni,nj,nk)
    logical,intent(in),optional :: ghosts
    real(8)                     :: each,all
    logical                     :: flag = .false.
    integer                     :: d

    if(present(ghosts)) flag = .not.ghosts
    if(flag) then
       each = sum(abs(s(is:ie,js:je,ks:ke)))
    else
       each = sum(abs(s))
    end if
    call mympi_sum(each,all)

    write(io_log,'(a5,e16.8)') name,all
    call log_print

  end subroutine io_log_scalar
!
!----------------------------------------------------------------
end module inout
