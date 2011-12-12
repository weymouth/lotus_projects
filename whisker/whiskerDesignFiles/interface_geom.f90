!-------------------------------------------------------!
!------------------ Whisker motion ---------------------!
!-------------------------------------------------------!
program design_geom
  use analytic
  use path
  implicit none
  character(4)      :: action,query,string
  character(40)     :: name
  type(pathSet)     :: motion(3)
  logical           :: motion_on
  type(set),target  :: body,fint,velo
  logical,target    :: body_on=.false.,fint_on=.false.,velo_on=.false.
  type(set),pointer :: geom
  logical,pointer   :: geom_on
  integer           :: d,stat
!
! -- CL interface
  write(*,*) "Prompt gives allowed inputs. Other text will exit"
  CL_loop: do
     write(*,*) "action: load, display, make, save"
     read(*,*,IOSTAT=stat) action
     if(stat.ne.0) exit CL_loop
     action_case: select case(action)
!
! -- load file
     case('load')
        write(*,*) 'geom file to load'
        read(*,*,IOSTAT=stat) name
        open(7,file=name)
        call geomRead(7)
        write(*,*) 'loaded ',name
        close(7)
!
! -- display or make
     case('disp','make')
        write(*,*) "query: path, body, fint, velo"
        read(*,*,IOSTAT=stat) query
        write(*,*) query        
!
! -- path
        if(query.eq.'path') then
           if(action.eq.'disp') then
              if(motion_on) then
                 do d=1,3
                    call pathSet_write(6,motion(d))
                 end do
              else
                 write(*,*) 'not defined'
              end if
           else
              call pathMake
           end if
!
! -- geom
        else
           select case(query)
           case('body')
              geom => body
              geom_on => body_on
           case('fint')
              geom => fint
              geom_on => fint_on
           case('velo')
              geom => velo
              geom_on => velo_on
           case default
              cycle CL_loop
           end select
           if(action.eq.'disp') then
              if(geom_on) then
                 call set_write(6,geom)
              else
                 write(*,*) 'not defined'
              end if
           else
              call setMake(geom)
              geom_on = .true.
           end if
        end if
!
! -- save the current state
     case('save')
        write(*,*) 'geom file to load'
        read(*,*,IOSTAT=stat) name
        open(7,file=name)
        call geomWrite(7)
        write(*,*) 'saved ',name
        close(7)
     case default
        exit CL_loop
     end select action_case
  end do CL_loop
  rewind(7)
!
! -- write out
  

contains
  subroutine geomRead(fnum)
    integer fnum
    motion_on=.false.;body_on=.false.;fint_on=.false.;velo_on=.false.
    read_loop: do
       read(fnum,*,IOSTAT=stat) string
       if(stat.ne.0) return
       select case(string)
       case('path')
          do d=1,3
             call pathSet_read(fnum,motion(d))
          end do
          motion_on = .true.
       case('body')
          call set_read(fnum,body)
          body_on = .true.
       case('fint')
          call set_read(fnum,fint)
          fint_on = .true.
       case('velo')
          call set_read(fnum,velo)
          velo_on = .true.
       case default
          exit read_loop
       end select
    end do read_loop
  end subroutine geomRead

  subroutine pathMake
    real(8)   :: r(3)
    if(.not.motion_on) then 
       do d=1,3
          motion(d) = .set.line(0,0)
       end do
    end if
    make_loop: do
       write(*,*) "path axis: 1, 2, 3."
       read(*,*,IOSTAT=stat) d
       if(stat.ne.0) return
       write(*,*) " create: line, trig."
       read(*,*,IOSTAT=stat) string
       if(stat.ne.0) return
       select case(string)
       case('line')
          write(*,*) "slope,intercept."
          read(*,*,IOSTAT=stat) r(:2)
          if(stat.ne.0) return
          motion(d) = .set.line(r(1),r(2))
       case('trig')
          write(*,*) "amplitude,frequency(deg),phase shift(deg)."
          read(*,*,IOSTAT=stat) r(:3)
          if(stat.ne.0) return
          motion(d) = .set.trig(r(1),r(2)*2.*acos(-1.),r(3)*2.*acos(-1.))
       case default
          return   
       end select
       write(*,*) 'new path ='
       call pathSet_write(6,motion(d))
       motion_on = .true.
    end do make_loop
  end subroutine pathMake

  recursive subroutine setMake(geom)
    type(nurbsFlags) :: flags
    type(set) :: geom,a,b
    real(8)   :: r(6)
    integer   :: i(2)

    write(*,*) "operation: 1) or, 2) and, 3) difference, 4) compliment, 5) define geom"
    read(*,*,IOSTAT=stat) d
    if(stat.ne.0) return
    if(d.eq.5) then
       write(*,*) "geom: cylinder, sphere, wave, plane, nurbs"
       read(*,*,IOSTAT=stat) string
       if(stat.ne.0) return
       select case(string)
       case('cyli')
          write(*,*) "flag,spanwise axis, radius, center(3)."
          read(*,*,IOSTAT=stat) i(:2),r(:4)
          geom = .set.cylinder(i(1),i(2),r(1),r(2:4),0,0)
       case('sphe')
          write(*,*) "flag, radius, center(3)."
          read(*,*,IOSTAT=stat) i(1),r(:4)
          geom = .set.sphere(i(1),r(1),r(2:4),0,0)
       case('wave')
          write(*,*) "flag, normal axis, average height, amplitude, wave length(2), phase(2)(deg)."
          read(*,*,IOSTAT=stat) i(:2),r(:6)
          geom = .set.wave(i(1),i(2),r(1),r(2),r(3:4),r(5:6)*2.*acos(-1.))
       case('plan')
          write(*,*) "flag, normal vector(3), intercept(3)."
          read(*,*,IOSTAT=stat) i(1),r(:6)
          geom = .set.plane(i(1),r(1:3),r(4:6),0,0)
       case('nurb')
          write(*,*) "flag, IGS file name, scale(3), rotate(3), shift(3), trim?, blank surface index(10)."
          read(*,*,IOSTAT=stat) flags
          geom = .set.flags
       end select
    else
       write(*,*) 'define parent a'
       call setMake(a)
       if(d.le.3) then
       write(*,*) 'define parent a'
          call setMake(b)
          geom = set_operation(d,a,b)
       else
          geom = set_operation(d,a)
       end if
    end if
    call set_write(6,geom)
  end subroutine setMake

  subroutine geomWrite(fnum)
    integer fnum
    if(motion_on) then
       write(fnum,'(a4)') 'path'
       do d=1,3
          call pathSet_write(fnum,motion(d))
       end do
    end if
    if(body_on) then
       write(fnum,'(a4)') 'body'
       call set_write(fnum,body)
    end if
    if(fint_on) then
       write(fnum,'(a4)') 'fint'
       call set_write(fnum,fint)
    end if
    if(velo_on) then
       write(fnum,'(a4)') 'velo'
       call set_write(fnum,velo)
    end if
  end subroutine geomWrite

end program design_geom
