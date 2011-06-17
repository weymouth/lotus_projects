!-------------------------------------------------------!
!           Analytic Geometry Routines                  !
!                                                       ! 
! The main purpose of this module is to allow for the   !
!  creation of analytic distance functions. The method  ! 
!  is to create types for basic shapes (wave, sphere,   !
!  cylinder, and plane) routines to get distances for   !
!  those types and routines to combine those distance   ! 
!  functions, allowing the creation of compound geoms.  ! 
! The key bit of cleverness is to create a 'set' type   ! 
!  which points to a distance function operation and    ! 
!  1 or 2 parent sets. This allows an arbitrarily       ! 
!  complex compound geometry to be coded in a single    ! 
!  line and stored in a single set.                     ! 
! Prevalent use of operator overloading allows these    ! 
!  geoms to be created with a few commands.             !
!                                                       ! 
!   c = .set.a  -> creates set c pointing to shape a    ! 
!   c = a.or.b  -> returns the union of the two sets    ! 
!   c = a.and.b -> returns the intersection             ! 
!   c = a - b   -> returns the difference               ! 
!   c = .not.a  -> returns the complement of a          ! 
!                                                       ! 
! To query the properties of the shape, use:            ! 
!                                                       ! 
!   p = a.at.x  -> returns the props of a at point x    !
!                                                       ! 
! The property type contains the distance, normal,      ! 
!  velocity, curvature, and an integer 'flag' for the   ! 
!  shape which is arbitrary but used in BDIM to signal  ! 
!  the shape dimensionality and boundary condition.     ! 
!                                                       ! 
! The module is serial and has no dependancies. Other   !
!  examples and commands are given in program test.     ! 
!                                                       ! 
!-------------------------------------------------------!
!
module analytic
  implicit none
  real(8),private,parameter :: pi = 3.14159265359, zero = 0., one = 1.
!
! -- basic shape types
!
  type cylinder
     integer :: flag,axis
     real :: radius,center(3),drad,dcen(3)
  end type cylinder
  type sphere
     integer :: flag
     real :: radius,center(3),drad,dcen(3)
  end type sphere
  type wave
     integer :: flag,axis
     real :: height,amp,length(2),phase(2)
  end type wave
  type plane
     integer :: flag
     real :: norm(3),center(3),dint,dcen(3)
  end type plane
!
! -- type(set) holds instruction for making compound shapes
!
  type set
     private
     integer :: operation    ! operation pointer
                             ! = 1 -> .or.
                             ! = 2 -> .and.
                             ! = 3 ->  -
                             ! = 4 -> .not.
                             ! = 5 -> .set.
!
!   if(operation<5) then point to the parent sets
     type(set),pointer :: a  ! parent
     type(set),pointer :: b  ! parent
!
!   if(operation=5) then point to the shape and type
     integer                :: type ! type pointer
     type(cylinder),pointer :: cyl  ! type = 1
     type(sphere)  ,pointer :: sph  ! type = 2
     type(wave)    ,pointer :: wav  ! type = 3
     type(plane)   ,pointer :: pln  ! type = 4
  end type set
!
! -- type(property)
!
  type prop
     integer :: flag
     real(8) :: distance,normal(3),velocity(3),kappa(2)
  end type prop
!
! -- operator overloading
!
  interface operator(.at.)
     module procedure plane_prop,cylinder_prop,sphere_prop,wave_prop,set_prop
  end interface
  interface operator(.or.)
     module procedure union_prop,union_set
  end interface
  interface operator(.and.)
     module procedure intersection_prop,intersection_set
  end interface
  interface operator(.not.)
     module procedure compliment_prop,compliment_set
  end interface
  interface operator(-)
     module procedure difference_prop,difference_set
  end interface
  interface operator(.set.)
     module procedure set_cyl,set_wav,set_sph,set_pln
  end interface
contains
!
! -- property query functions for basic types, p=a.at.x
!
  pure function cylinder_prop(a,x) result(p)
    real(8),intent(in)        :: x(3)
    type(cylinder),intent(in) :: a
    type(prop) :: p
    real(8)    :: v(3),m
    v = x-a%center
    v(a%axis) = 0
    m = sqrt(sum(v**2))
    p%distance = m-a%radius
    p%normal   = v/merge(m,one,abs(m).gt.1e-9)
    p%velocity = a%drad*p%normal+a%dcen
    p%kappa    = (/0.,1./a%radius/)
    p%flag     = a%flag
  end function cylinder_prop
  pure function sphere_prop(a,x) result(p)
    real(8),intent(in)      :: x(3)
    type(sphere),intent(in) :: a
    type(prop) :: p
    real(8)    :: v(3),m
    v = x-a%center
    m = sqrt(sum(v**2))
    p%distance = m-a%radius
    p%normal   = v/merge(m,one,abs(m).gt.1e-9)
    p%velocity = a%drad*p%normal+a%dcen
    p%kappa    = (/1.,2./)/a%radius
    p%flag     = a%flag
  end function sphere_prop
  pure function plane_prop(a,x) result(p)
    real(8),intent(in)     :: x(3)
    type(plane),intent(in) :: a
    type(prop) :: p
    p%normal   = a%norm
    p%distance = sum(a%norm*(x-a%center))
    p%velocity = a%dcen+p%normal*a%dint
    p%kappa    = 0
    p%flag     = a%flag
  end function plane_prop
  pure function wave_prop(a,x) result(p)
    real(8),intent(in)    :: x(3)
    type(wave),intent(in) :: a
    type(prop) :: p
    real(8)    :: n(3),y(2),f
    integer    :: i
! note: this uses a linear approximation !
    i = a%axis-1
    y(:i) = x(:i)/a%length(:i)*2.*pi-a%phase(:i)
    f     =  a%amp*product(cos(y(:i)))
    n(:i) = 2.*pi/a%length(:i)*tan(y(:i))*f
    n(a%axis) = 1; n(a%axis+1:) = 0
    p%distance = x(a%axis)-a%height-f
    p%normal   = n/sqrt(sum(n**2))
    p%velocity = 0
    p%kappa    = 0
    p%flag     = a%flag
  end function wave_prop
!
! -- property operations, p = operation(a,b)
!
  pure function union_prop(a,b) result(p) 
    type(prop),intent(in) :: a,b
    type(prop) :: p
    p = a
    if(b%distance.lt.a%distance) p = b
  end function union_prop
  pure function intersection_prop(a,b) result(p)
    type(prop),intent(in) :: a,b
    type(prop) :: p
    p = a
    if(b%distance.gt.a%distance) p = b
  end function intersection_prop
  pure function compliment_prop(a) result(p)
    type(prop),intent(in) :: a
    type(prop) :: p
    p = a
    p%distance = -a%distance
    p%normal   = -a%normal
  end function compliment_prop
  pure function difference_prop(a,b) result(p)
    type(prop),intent(in) :: a,b
    type(prop) :: p
    p = a.and..not.b
  end function difference_prop
!
! -- set operations, c = operation(a,b) (operation < 5)
!
   function union_set(a,b) result(c)
     type(set),target,intent(in) :: a,b
     type(set),pointer :: c
     allocate(c)
     c%operation = 1 !-> .or.
     c%a         => a
     c%b         => b
   end function union_set
   function intersection_set(a,b) result(c)
     type(set),target,intent(in) :: a,b
     type(set),pointer :: c
     allocate(c)
     c%operation = 2 !-> .and.
     c%a         => a
     c%b         => b
   end function intersection_set
   function difference_set(a,b) result(c)
     type(set),target,intent(in) :: a,b
     type(set),pointer :: c
     allocate(c)
     c%operation = 3 !->  -
     c%a         => a
     c%b         => b
   end function difference_set
   function compliment_set(a) result(c)
     type(set),target,intent(in) :: a
     type(set),pointer :: c
     allocate(c)
     c%operation = 4 !-> .not.
     c%a         => a
   end function compliment_set
!
! -- identity sets, c = .set.a (operation = 5)
!
  function set_cyl(a) result(c)
     type(cylinder),target,intent(in) :: a
     type(set),pointer :: c
     allocate(c)
     c%operation = 5
     c%type      = 1 !-> cylinder
     c%cyl       => a 
   end function set_cyl
  function set_sph(a) result(c)
     type(sphere),target,intent(in) :: a
     type(set),pointer :: c
     allocate(c)
     c%operation = 5
     c%type      = 2 !-> sphere
     c%sph       => a 
   end function set_sph
  function set_wav(a) result(c)
     type(wave),target,intent(in) :: a
     type(set),pointer :: c
     allocate(c)
     c%operation = 5
     c%type      = 3 !-> wave
     c%wav       => a 
   end function set_wav
  function set_pln(a) result(c)
     type(plane),target,intent(in) :: a
     type(set),pointer :: c
     allocate(c)
     c%operation = 5
     c%type      = 4 !-> plane
     c%pln       => a 
   end function set_pln
!
! -- query function for sets, p = c.at.x = c%operation(c%a.at.x,c%b.at.x)
!
   pure recursive function set_prop(c,x) result(p)
    real(8),intent(in)   :: x(3)
    type(set),intent(in) :: c
    type(prop) :: a,b,p
!
! -- get a
    if(c%operation.lt.5) then
       a = c%a.at.x
    else if(c%type.eq.1) then
       a = c%cyl.at.x
    else if(c%type.eq.2) then
       a = c%sph.at.x
    else if(c%type.eq.3) then
       a = c%wav.at.x
    else if(c%type.eq.4) then
       a = c%pln.at.x
    end if
!
! -- get b
    if(c%operation.le.3) then
       b = c%b.at.x
    end if
!
! -- get p = operation(a,b)
    if     (c%operation.eq.1) then
       p = a.or.b
    else if(c%operation.eq.2) then
       p = a.and.b
    else if(c%operation.eq.3) then
       p = a-b
    else if(c%operation.eq.4) then
       p = .not.a
    else if(c%operation.eq.5) then
       p = a
    end if
    
  end function set_prop
!
! -- write/read a set
!
  recursive subroutine set_write(io_num,c)
    integer,intent(in)   :: io_num
    type(set),intent(in) :: c

    write(io_num,*) c%operation
    if(c%operation.le.4) then
       call set_write(io_num,c%a)
    else 
       write(io_num,*) c%type
       if(c%type.eq.1) then
          write(io_num,*) c%cyl
       else if(c%type.eq.2) then
          write(io_num,*) c%sph
       else if(c%type.eq.3) then
          write(io_num,*) c%wav
       else if(c%type.eq.4) then
          write(io_num,*) c%pln
       end if
    end if
    if(c%operation.le.3) then
       call set_write(io_num,c%b)
    end if
    
  end subroutine set_write
!
  recursive subroutine set_read(io_num,c)
    integer,intent(in)    :: io_num
    type(set),intent(out) :: c
    logical               :: debug=.false.

    read(io_num,*) c%operation
    if(c%operation.le.4) then
       if(debug) write(io_num+1,*) 'push a'
       allocate(c%a)
       call set_read(io_num,c%a)
       if(debug) call set_write(io_num+1,c%a)
       if(debug) write(io_num+1,*) 'pop a'
    else 
       read(io_num,*) c%type
       if(c%type.eq.1) then
          allocate(c%cyl)
          read(io_num,*) c%cyl
       else if(c%type.eq.2) then
          allocate(c%sph)
          read(io_num,*) c%sph
       else if(c%type.eq.3) then
          allocate(c%wav)
          read(io_num,*) c%wav
       else if(c%type.eq.4) then
          allocate(c%pln)
          read(io_num,*) c%pln
       end if
    end if
    if(c%operation.le.3) then
       if(debug) write(io_num+1,*) 'push b'
       allocate(c%b)
       call set_read(io_num,c%b)
       if(debug) call set_write(io_num+1,c%b)
       if(debug) write(io_num+1,*) 'pop b'
    end if
    
  end subroutine set_read
end module analytic
!!$!
!!$! -- test
!!$!
!!$module temp
!!$  use analytic
!!$  implicit none
!!$  type(plane)    :: a,b
!!$  type(cylinder) :: c
!!$  type(sphere)   :: d
!!$  type(wave)     :: e
!!$  type(set)      :: f,g,h
!!$contains
!!$  subroutine init_temp
!!$    type(set) :: z
!!$
!!$    a = plane(1,1,0,0,0)
!!$    b = plane(2,2,(/6,0,0/),0,0)
!!$    c = cylinder(3,2,1,(/7,0,0/),0,0)
!!$    d = sphere(4,5,0,0,0)
!!$    e = wave(5,2,0,2,6,0)
!!$    f = .set.a.or..not..set.b
!!$    g = .set.d-.set.e
!!$    h = .set.c.and.f.or.g
!!$
!!$    z = f.and.g-.set.c
!!$    call set_write(7,z)
!!$    close(7)
!!$
!!$  end subroutine init_temp
!!$end module temp
!!$program test
!!$  use temp
!!$  implicit none
!!$  type(prop) :: ax,bx,cx,dx,ex
!!$  real(8)    :: x(3)
!!$  type(set)  :: z
!!$
!!$  call init_temp()
!!$  x = (/5,5,0/)
!!$
!!$  ax = a.at.x
!!$  bx = b.at.x
!!$  cx = c.at.x
!!$  dx = d.at.x
!!$  ex = e.at.x
!!$  
!!$  print '("Test each basic shape (result should be 1,2,3,4,5)")' 
!!$  print 1, ax,bx,cx,dx,ex
!!$  print '("Test compound shapes against sets (should match)")' 
!!$  print 1, ax.or..not.bx,f.at.x
!!$  print 1, dx-ex,g.at.x
!!$  print 1, cx.and.(f.at.x).or.(g.at.x),h.at.x
!!$
!!$  print '("Test read/write (should match)")' 
!!$  call set_read(7,z)
!!$  print 1, (f.at.x).and.(g.at.x)-cx,z.at.x
!!$
!!$1 format(i3,9f6.2)
!!$end program test
