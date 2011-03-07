!-------------------------------------------------------!
!----------------- IGES Routines -----------------------!
!-------------------------------------------------------!
! -- Gabriel Weymouth 
! -- MIT, 07/05/07
! -- DSS, 10/22/08  Added read for 144,142,102,126
!-------------------------------------------------------!
!
! -- Convert NURBS surfaces from IGES to SRF format
!
! Note: IGES files are designed to store very general types of geometric data. 
! Note: The Parameters of a specific entity are given in the 'P' section. The
! Note: file is fixed width (80), but the data from positions 65 to 80 is not
! Note: parameter data. Additionally, the parameters are not written in any
! Note: fixed format and will run over many lines. These factors make typical
! Note: file reading impossible. 
!
! Note: To overcome this, the routine below takes advantage a few factors.
! Note: (1) Each line is marked with it's section character in the 73rd
! Note:     position.
! Note: (2) A Directory of parameter information is written in the 'D' section,
! Note:     inlcuding the number of lines each entity is written over.
! Note: (3) The NURBS surface entities are marked by the string '128' in the
! Note:     'D' and 'P' sections.
! Note: The final trick used below is that a character string may be read as an
! Note: internal file to translate it from text to numeric data.
!
! Note: The routine below checks for surfaces in the directory and gets the
! Note: number of lines. Then it reads surface data from section 'P' as text,
! Note: concatenating the first 64 characters in each line together to form one
! Note: long text string. Then the program reads from that text string as an
! Note: internal file to get the numeric parameter data for each surface.
! Note: Finally, this data is written in a scratch file in SRF format, allowing
! Note: it to be read by the NURBS module.
!
module iges
contains
  subroutine iges_to_srf(inum,onum,n)
    implicit none
    integer,parameter :: maxn=1000,maxlines=2000
    integer,intent(in) :: inum,onum
    integer,intent(out):: n
    integer :: ind(3,maxn),ierr,junk,nu,nv,ordu,ordv,i,j,num
    character(80) :: line
    character(64*maxlines) :: record
    real(8),allocatable :: knotu(:),knotv(:)
    real(8),allocatable :: wght(:,:),xcp(:,:,:)
    integer,allocatable :: work(:)
    logical :: debug = .false.
    if(debug) print '("inum,onum",2i4)',inum,onum
!
! -- Part I : Get the record sizes from the 'D'=Directory section
!
    n=0
    read_directory: do
!
! -- read in a line of the file as a text string
       read(inum,'(A80)',IOSTAT=ierr) line
       if(ierr.ne.0) stop 'iges: read err 1'
!
! -- check for directory

       if(line(73:73).eq.'D') then
!
! -- check for a readable entity (type 128,126,102,144,142)
          read(line(6:8),*) i
          if(debug) print '(a80)',line
          if(i.eq.128.or.i.eq.126.or.i.eq.102.or.i.eq.144.or.i.eq.142) then
             n = n+1
             if(n.gt.maxn) stop 'iges: too many entities'
!
! -- remember entity information
             ind(1,n) = i
             read(line(74:),*) ind(2,n)
             read(inum,*,IOSTAT=ierr) junk,junk,junk,ind(3,n)
             if(ierr.ne.0) stop 'iges: read err 2'
             if(ind(3,n).gt.maxlines) stop 'iges: record too long'
             if(debug) print '(" ent,type,position,length",4i4)',n,ind(:,n)
          end if
!
! -- check for end of directory
       else if(line(73:73).eq.'P') then
          exit read_directory
       end if
    end do read_directory
!
! -- write directory info
    ind(2,:) = ind(2,:)+onum
    if(debug) then
       print '("count of 128",i6)',count(ind(1,:).eq.128)
       print '("count of 126",i6)',count(ind(1,:).eq.126)
       print '("count of 102",i6)',count(ind(1,:).eq.102)
       print '("count of 142",i6)',count(ind(1,:).eq.142)
       print '("count of 144",i6)',count(ind(1,:).eq.144)
       pause
    end if
!
! -- Part II: Read/write data from the 'P'=Parameter section
!
    backspace(inum); n = 0
    read_parameters: do
       read(inum,'(A80)',IOSTAT=ierr) line
       if(ierr.ne.0) stop 'iges: read err 3'
       if(debug) print '(a80)',line
!                     !
! -- check for 128 -- !
!                     !
       if(line(73:73).eq.'P'.and.line(1:3).eq.'128') then
          read(line(68:72),*) nu
          if(nu.ne.ind(2,n+1)-onum) cycle
          n = n+1
!
! -- read data from line as numeric
          read(line,*) junk,nu,nv,ordu,ordv
          nu = nu+1; ordu = ordu+1
          nv = nv+1; ordv = ordv+1
          if(debug) print '(" surf",4i4)',nu,nv,ordu,ordv
!
! -- allocate arrays
          allocate(knotu(nu+ordu),knotv(nv+ordv),&
                   wght(nu,nv),xcp(3,nu,nv))
!
! -- read and store the complete record as one text string
          record(1:64) = line(1:64)
          do i=2,ind(3,n)
             read(inum,'(A80)',IOSTAT=ierr) line
             if(ierr.ne.0) stop 'iges: read err 4'
             record(i*64-63:i*64)= line(1:64)
          end do
!
! -- read data from record as numeric
          read(record,*) (junk,i=1,10),knotu,knotv,wght,xcp
!
! -- write surface data to an SRF scratch file
          num = ind(2,n)
          if(debug) then
             open(num)
          else
             open(num,status='scratch')
          end if
          write(num,'(i6,",",i6)') ordu,ordv
          write(num,'(i6,",",i6)') nu,nv
          write(num,'(99f7.3)') knotu
          write(num,'(99f7.3)') knotv
          do i=1,nu
             do j=1,nv
                write(num,'(4e16.8)') xcp(:,i,j),wght(i,j)
             end do
          end do
          rewind(num)
          deallocate(knotu,knotv,xcp,wght)
!                     !
! -- check for 126 -- !
!                     !
       else if(line(73:73).eq.'P'.and.line(1:3).eq.'126') then
          read(line(68:72),*) nu
          if(nu.ne.ind(2,n+1)-onum) cycle
          n = n+1
!
! -- read data from line as numeric
          read(line,*) junk,nu,ordu,i
          if(i.ne.1) cycle ! non planar (not interested)
          nu = nu+1; ordu = ordu+1
          if(debug) print '(" curve",4i4)',nu,ordu
!
! -- allocate arrays
          allocate(knotu(nu+ordu),wght(nu,1),xcp(3,nu,1))
!
! -- read and store the complete record as one text string
          record(1:64) = line(1:64)
          do i=2,ind(3,n)
             read(inum,'(A80)',IOSTAT=ierr) line
             if(ierr.ne.0) stop 'iges: read err 4'
             record(i*64-63:i*64)= line(1:64)
          end do
!
! -- read data from record as numeric
          read(record,*) (junk,i=1,7),knotu,wght,xcp
!
! -- write surface data to an CRV scratch file
          num = ind(2,n)
          if(debug) then
             open(num)
          else
             open(num,status='scratch')
          end if
          write(num,'(i6)') ordu
          write(num,'(i6)') nu
          write(num,'(99f6.3)') knotu
          do i=1,nu
             write(num,'(4e16.8)') xcp(1:2,i,1),wght(i,1)
          end do
          rewind(num)
          deallocate(knotu,xcp,wght)
!                     !
! -- check for 102 -- !
!                     !
       else if(line(73:73).eq.'P'.and.line(1:3).eq.'102') then
          read(line(68:72),*) nu
          if(nu.ne.ind(2,n+1)-onum) cycle
          n = n+1
!
! -- read data from line as numeric
          read(line,*) junk,nu
          if(debug) print '(" lines",4i4)',nu
!
! -- allocate arrays
          allocate(work(nu))
!
! -- read and store the complete record as one text string
          record(1:64) = line(1:64)
          do i=2,ind(3,n)
             read(inum,'(A80)',IOSTAT=ierr) line
             if(ierr.ne.0) stop 'iges: read err 4'
             record(i*64-63:i*64)= line(1:64)
          end do
          do i=1,ind(3,n)*64
             if(record(i:i).eq.';') record(i:i) = ' '
          end do
!
! -- read data from record as numeric
          read(record,*,IOSTAT=ierr) (junk,i=1,2),work
!
! -- write pointer data to a scratch file
          num = ind(2,n)
          if(debug) then
             open(num)
          else
             open(num,status='scratch')
          end if
          write(num,'(i6)') nu
          write(num,'(i6)') work+onum
          rewind(num)
          deallocate(work)
!                     !
! -- check for 142 -- !
!                     !
       else if(line(73:73).eq.'P'.and.line(1:3).eq.'142') then
          read(line(68:72),*) nu
          if(nu.ne.ind(2,n+1)-onum) cycle
          n = n+1
!
! -- read data from line as numeric
          read(line,*) junk,junk,nu,nv
          if(debug) print '(" surfcurve",2i4)',nu,nv
!
! -- write pointer data to a scratch file
          num = ind(2,n)
          if(debug) then
             open(num)
          else
             open(num,status='scratch')
          end if
          write(num,'(2i6)') nu+onum,nv+onum
          rewind(num)
!                     !
! -- check for 144 -- !
!                     !
       else if(line(73:73).eq.'P'.and.line(1:3).eq.'144') then
          read(line(68:72),*) nu
          if(nu.ne.ind(2,n+1)-onum) cycle
          n = n+1
!
! -- read data from line as numeric
          read(line,*) junk,ordu,nu,nv
          allocate(work(nu+nv))
          if(debug) print '(" trimsurf",3i4)',ordu,nu,nv
!
! -- read pointers
          do i=1,64
             if(line(i:i).eq.';') line(i:i) = ' '
          end do
          read(line,*) (junk,i=1,4),work
!
! -- write a scratch TRM file
          num = ind(2,n)
          if(debug) then
             open(num)
          else
             open(num,status='scratch')
          end if
          if(.not.any(ind(1,:).eq.128 &
                 .and.ind(2,:).eq.ordu+onum)) then
             stop '144 using non-128'
          end if
          write(num,'(i6)') ordu+onum ! surface
!
! -- look into the 142 files
          nu = 0
          do i=1,size(work)
             read(work(i)+onum,*) junk,ordv
             if(junk.ne.ordu+onum) stop '144 error 1'
!
! -- count curves in the 102 files and write total
             read(ordv,*) nv
             nu = nu+nv            
             rewind(ordv); rewind(work(i)+onum)
          end do
          write(num,'(i6)') nu     
!
! -- go through again and point directly to the curves
          do i=1,size(work)
             read(work(i)+onum,*) junk,ordv
             read(ordv,*) nv
             do j=1,nv
                read(ordv,*) junk
                write(num,'(i6)') junk
             end do
             rewind(ordv); rewind(work(i)+onum)
          end do
          rewind(num)
          if(debug.and.num.eq.171) pause
          deallocate(work)
!
! -- check for end of file
       else if(line(73:73).eq.'T') then
          exit read_parameters
       end if
    end do read_parameters
!
! -- Part III: Point to surface files
!
    i = count(ind(1,:).eq.144)
    if(debug) then
       open(onum)
    else
       open(onum,status='scratch')
    end if
!
! -- No trimmed surfaces, just print pointers to 128
    if(i.eq.0) then
       i = count(ind(1,:).eq.128)
       write(onum,'(i6)') i
       write_128: do i=1,n
          if(ind(1,i).eq.128) then
             write(onum,'(i6)') ind(2,i)
          end if
       end do write_128
       n = 128 ! flag !
!
! -- Print pointers to 144
    else
       write(onum,'(i6)') i
       write_144: do i=1,n
          if(ind(1,i).eq.144) then
             write(onum,'(i6)') ind(2,i)
          end if
       end do write_144
       n = 144 ! flag !
    end if
    rewind(onum)
    
  end subroutine iges_to_srf
end module iges
!!$!
!!$!-------------------------------------------------------!
!!$! -- Test routine for IGES module
!!$!
!!$program test_iges
!!$  use iges
!!$  implicit none
!!$  integer :: a,b
!!$  open(9,file='TMS_Surface_2.IGS')
!!$!  open(9,file='TMS_shape1.IGS')
!!$  call iges_to_srf(9,20,a)
!!$  write(*,*) a
!!$end program test_iges
