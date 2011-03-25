
  program tec2vtk

  implicit none
  real(8), dimension(:,:,:),allocatable :: x,y,z,f
  real(8) :: x0,y0,z0,dx,dy,dz
  integer :: ni,nj,nk
  real(8) :: time

  character(120) :: line
  character(80) :: fin,fout
  character(10) :: timestr,timestr1
  integer :: i,ii

  read(4001) ni,nj,nk
  read(4002)
  read(4003)

  allocate(x(ni,nj,nk),y(ni,nj,nk),z(ni,nj,nk),f(ni,nj,nk))
  read(4001) x
  read(4002) y
  read(4003) z

  x0=x(1,1,1); y0=y(1,1,1); z0=z(1,1,1)
  dx=x(2,1,1)-x0; dy=y(1,2,1)-y0; dz=z(1,1,2)-z0

  call getarg(1,fin)
  call getarg(2,timestr1)
  read(timestr1,'(i4)') ii

  print*,'num steps',ii
  i=0
  do i=0,ii

  if(i.lt.10) write(timestr,'("000",i1)') i
  if(i.lt.100.and.i.ge.10) write(timestr,'("00",i2)') i
  if(i.lt.1000.and.i.ge.100) write(timestr,'("0",i3)') i
  if(i.ge.1000) write(timestr,'(i4)') i

  fout=trim(fin)//'.t.'//trim(timestr)
  write(*,'(a20)') fout
  call writevtkpoints(x0,y0,z0,ni,nj,nk,dx,dy,dz,fout)
  !call writevtkgrid(x,y,z,fout)

  read(5000+i)
  read(5000+i) f
  call writevtkscalar(f,'fraction',fout)
  close(2)

  end do

  deallocate(x,y,z,f)
 
100 format("ZONE SOLUTIONTIME = ",e12.4,", I = ",i4,", J =  ",i4,", K =  ",i4)
101 format("ZONE SOLUTIONTIME = ",e12.4)
 9 continue

contains
  subroutine writevtkpoints(x0,y0,z0,ni,nj,nk,dx,dy,dz,fout)

  implicit none
  real(8), intent(in) :: x0,y0,z0,dx,dy,dz
  character(*) :: fout
  integer :: ni,nj,nk
  logical :: op

  inquire(2,opened=op)
  if(.NOT. op) then
    open(2,file=fout,status='UNKNOWN')
    write(2,1) 
    write(2,*) 'Volume fraction'
    write(2,*) 'ASCII'
    write(2,*) 'DATASET STRUCTURED_POINTS'
    write(2,*) 'DIMENSIONS ',ni,nj,nk
    write(2,2) x0,y0,z0
    write(2,3) dx,dy,dz
  endif
    write(2,*) 'POINT_DATA ',ni*nj*nk
1 format('# vtk DataFile Version 3.0')
2 format(' ORIGIN ',3e12.4)
3 format('SPACING ',3e12.4)

  end subroutine writevtkpoints

  subroutine writevtkgrid(x,y,z,fout)

  implicit none
  real(8), dimension(:,:,:) :: x,y,z
  character(*) :: fout
  integer :: ni,nj,nk,i,j,k
  logical :: op

  ni=size(x,1); nj=size(x,2); nk=size(x,3)
  inquire(2,opened=op)
  if(.NOT. op) then
    open(2,file=fout,status='UNKNOWN')
    write(2,1) 
    write(2,*) 'Volume fraction'
    write(2,*) 'ASCII'
    write(2,*) 'DATASET STRUCTURED_GRID'
    write(2,*) 'DIMENSIONS ',ni,nj,nk
    write(2,*) 'POINTS ',ni*nj*nk,' float'
    do k=1,nk
	do j=1,nj
	do i=1,ni
	write(2,'(3f9.5)') x(i,j,k),y(i,j,k),z(i,j,k)
	end do
	end do
    end do
  endif
    write(2,*) 'POINT_DATA ',ni*nj*nk
1 format('# vtk DataFile Version 3.0')

  end subroutine writevtkgrid

  subroutine writevtkscalar(f,name,fout)

  implicit none
  real(8), dimension(:,:,:) :: f
  character(*) :: fout
  character(*) :: name
  integer :: ni,nj,nk,i,j,k
  logical :: op

  ni=size(f,1); nj=size(f,2); nk=size(f,3)
  inquire(2,opened=op)
  if(op) then
    write(2,*) 'SCALARS ',trim(name),' float 1'
    write(2,*) 'LOOKUP_TABLE default'
    write(2,'(f5.2)') f

!    do k=1,nk
!	do j=1,nj
!	do i=1,ni
!	write(2,'(f3.2)') f(i,j,k)
!	end do
!	end do
 !   end do
  else
	stop 'no grid data written'
  endif
    
  end subroutine writevtkscalar
  end program

! VARIABLES=x,y,z,p
!ZONE SOLUTIONTIME =   0.0000E+00, I =   36, J =    36, K =    36, F=BLOCK
! -0.2284E+01 -0.2224E+01 -0.2164E+01 -0.2105E+01 -0.2048E+01 -0.1992E+01 -0.1937E+01 -0.1883E+01 -0.1831E+01 -0.1780E+01 -0.1731E+01 -0.1682E+01 -0.1635E+01 -0.1589E+01 -0.1544E+01
! -0.1500E+01 -0.1457E+01 -0.1415E+01 -0.1374E+01 -0.1334E+01 -0.1295E+01 -0.1257E+01 -0.1220E+01 -0.1184E+01 -0.1149E+01 -0.1114E+01 -0.1081E+01 -0.1048E+01 -0.1016E+01 -0.9844E+00
! -0.9531E+00 -0.9219E+00 -0.8906E+00 -0.8594E+00 -0.8281E+00 -0.7969E+00 -0.2284E+01 -0.2224E+01 -0.2164E+01 -0.2105E+01 -0.2048E+01 -0.1992E+01 -0.1937E+01 -0.1883E+01 -0.1831E+01
