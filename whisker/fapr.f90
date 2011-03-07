module apriori
  use global
  use inout
  type(print_flags) :: pflags
contains
!
  subroutine init_apriori(u,p)
    implicit none
    real(8),dimension(ndims,ni,nj,nk),intent(in) :: u
    real(8),dimension(ni,nj,nk),intent(in)       :: p

    call io_read
    call io_read(f=pflags)
    call apriori_print(0,u,p)

  end subroutine init_apriori
!
  subroutine apriori_print(t,u,p)
    use mympi,   only: mympi_write_subarray,mympi_id
    use grid,    only: grid_position_array
    use utility, only: shift
    use freeint, only: freeint_get
    implicit none
    integer,intent(in)                           :: t
    real(8),dimension(ndims,ni,nj,nk),intent(in) :: u
    real(8),dimension(ni,nj,nk),intent(in)       :: p
    real(8),dimension(ni,nj,nk) :: v
    integer :: file,f,tmod,d
    real(8) :: time
!
! -- Turb file
    if(.not.pflags%prnt) return
    tmod = pflags%tmod
    if(mod(t,tmod).eq.0) then
       file = pflags%file
       time = time0+t*dt
       if(t.eq.0) then
          f = file
          do d=1,ndims
             v = grid_position_array(d,0); f = f+1
             write(file,1) f,time,merge('x',merge('y','z',d.eq.2),d.eq.1)
             call mympi_write_subarray(f,v)
          end do
       end if
       f = file+t/tmod+1000
       write(file,1) f,time,'f'
       call freeint_get(vof=v)
       call mympi_write_subarray(f,v)
       do d=1,ndims
          call shift(1,d,u(d,:,:,:),v)
          v = (v+u(d,:,:,:))*0.5
          f = f+1000
          write(file,1) f,time,merge('u',merge('v','w',d.eq.2),d.eq.1)
          call mympi_write_subarray(f,v)
       end do
       if(mympi_id().eq.0) call flush(file)
       f = f+1000
       write(file,1) f,time,'p'
       call mympi_write_subarray(f,p)
    end if

1   format("File = fort.",i5,", Time = ",e12.4,", Variable = ",a)
  end subroutine apriori_print
end module apriori
