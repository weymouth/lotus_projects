program maketec
  use blob
  implicit none
  character :: name
  integer   :: ni,nj,nk,logfile=4000,file,unit=10,ios
  real(8)   :: time,time0
  real(8),allocatable,dimension(:,:,:) :: f,u,v
  integer, allocatable, dimension(:,:,:) :: ilabel
  real(8),allocatable :: z(:),fr(:,:)
  logical :: bubbles,honly,full
  character(8) :: rtype,otype

  call getarg(1,rtype)
  call getarg(2,otype)

  bubbles=.false.
  honly=.false.
  full=.true.
  if(trim(rtype).eq.'bubbles') bubbles=.true.
  if(trim(otype).eq.'honly') honly=.true.
  if(honly) full=.false.

! init
  read(logfile,1) file,time0,name
  read(file) ni,nj,nk
  allocate(f(ni,nj,nk))
  allocate(ilabel(ni,nj,nk))

  if(honly) then
  	open(unit+1,file='apriori_blob.dat')
	open(unit+5,file='apriori_vol.dat')
  else
  open(unit,file='apriori.dat')
  open(unit+1,file='apriori_blob.dat')
	open(unit+3,file='apriori_label.dat')
  open(unit+5,file='apriori_vol.dat')
  write(unit,*) 'VARIABLES=x,y,z,f'
  write(unit,2) time0,ni,nj,nk
  endif

! x,y,z
	if(full) then
  do ios=1,3
     if(ios.ne.1) then
        read(logfile,1) file,time,name
        read(file) ni,nj,nk
     end if
     read(file) f
     write(unit,3) f
     close(file)
  end do
	endif

! loop in time
  ios = 0
  do while (ios.eq.0)

  ! read log
     read(logfile,1,iostat=ios) file,time,name
     if(ios.ne.0) stop 'EOF'

  ! read data
     if(name.eq.'p') cycle
     read(file) ni,nj,nk
     read(file) f
     close(file)

  ! write header
     if(time.ne.time0) then
	if(full) then
        write(unit,*) 'VARIABLES=x,y,z,f'
        write(unit,2) time,ni,nj,nk
        write(unit,*)'VARSHARELIST=([1-3]=1)'
	endif
	

        time0 = time
     end if


  ! get/write data
     if(name.eq.'f') then
        if(full) write(unit,4) f
     if(bubbles) then
        f = 1-f
     endif
        call blob_find(f,thresh=(/0.5,0.25,0.0/),time=time,d_num=unit+1,t_num=unit+5,iout=ilabel)
	write(unit+3,7) ilabel

     else if(name.eq.'u') then
     else if(name.eq.'v') then
     else if(name.eq.'w') then
     end if

  end do

1 format("File = fort.",i5,", Time = ",e12.4,", Variable = ",a)
2 format('ZONE SOLUTIONTIME = ',e12.4,', I =',i5, &
       ', J = ',i5,', K = ',i5,', F=BLOCK')
3   format(15e12.4)
4   format(15f6.2)

5 format('ZONE SOLUTIONTIME = ',e12.4,', I =',i5, &
       ', J = ',i5,', F=BLOCK')
6 format('fmax(t=',f5.1,', z=',f5.1,') =',e12.4)
7   format(15i6)
end program maketec
