program measure_blobs
  use blob
  implicit none
  character(20) :: string
  character(52) :: buff
  character :: name
  integer   :: ni,nj,nk,logfile=3000,file,unit=10,ios
  real(8)   :: time,time0,dt=0.05
  real(8),allocatable :: f(:,:,:),thresh(:)

! read arguments
  ni = command_argument_count()
  print *,"number of thresholds=",ni
  if(ni.ne.0) then
     allocate(thresh(ni))
     do nj=1,ni
        call getarg(nj,string)
        read(string,*,IOSTAT=ios) thresh(nj)
        if(ios.ne.0) stop 'bad arg'
     end do
  else
     allocate(thresh(1))
     thresh = 0.D0
  end if

! read log, size array
  read(logfile,'(A52)',iostat=ios) buff
  if(ios.ne.0) stop 'no apriori logfile'
  read(buff(13:17),*) file
  read(file) ni,nj,nk
  allocate(f(ni,nj,nk))
  open(unit+1,file='blob.dat')
  write(unit+1,*)'VARIABLES=index,vol,x,y,z,lx,ly,lz'
  open(unit+2,file='blob.csv')
  write(unit+2,*)'time,vol,x,y,z,lx,ly,lz'

! loop in time
  do
  ! read log
     read(logfile,'(A52)',iostat=ios) buff
     if(ios.ne.0) stop 'EOF'
     read(buff(13:17),*) file
     read(buff(27:38),*) time
     read(buff(52:52),*) name
     print *,file,time,name

  ! read data
     if(name.ne.'f') cycle
!     if(mod(time,dt).ne.0) cycle
     read(file) ni,nj,nk
     read(file) f
     close(file)

  ! find air blobs
     f = 1.D0-f
     call blob_find(f,thresh=thresh,time=time,d_num=unit+1,R_num=unit+2)
  end do
911 stop 'bye'
end program measure_blobs
