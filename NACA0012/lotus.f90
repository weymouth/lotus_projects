program flapping_NACA
 use fluidMod, only: fluid
  use bodyMod,  only: bodyUnion
  use mympiMod, only: init_mympi,mympi_end,mympi_rank
  use gridMod,  only: xg,composite
  use imageMod, only: display
  use geom_shape
  implicit none
!
! -- Define parameter, declare variables
  type(fluid)        :: flow
  type(bodyUnion)    :: geom
  type(set)          :: body1, body2
  real,parameter     :: c =64 , amp=c
  real,parameter     :: Re = 1.e4
  integer,parameter  :: n(3) = (/9.*c,6.*c,1./) ! number of points
  real,parameter     :: St1=0.4, St=0.4, pamp1 = atan(PI*St1)-PI/18. , pamp2 = atan(PI*St)-PI/18.
  logical            :: there = .false.
  logical :: root                 ! root processor
  integer:: b(3) = (/4,2,1/)   ! MPI blocks (product must equal n_procs)
!
! -- Initialize MPI (if MPI is OFF, b is set to 1)
  call init_mympi(ndims=2,set_blocks=b)
  root = mympi_rank()==0
!
! -- Initialize domain size and placement
  call xg(1)%stretch(n(1), -3.*c, -0.5*c, 6.*c, 10.*c, h_min=1., h_max=7.,prnt=root)
  call xg(2)%stretch(n(2), -5.*c, -2.*c, 2.*c, 5.*c, h_min=1., prnt=root)
!
! -- Initialize the sphere geometry (NACA0016)
body1=naca(c,thick=0.16 , pivot=0.25).map.init_rigid(6,p1).map.init_rigid(2,y1)
body2=naca(c,thick=0.16 , pivot=0.25).map.init_rigid(6,p2).map.init_rigid(2,y2).map.init_rigid(1,x)

call geom%add(body1)
call geom%add(body2)

!
! -- Initialize fluid
  call flow%init(n/b,geom,V=(/1.,0.,0./),nu=c/Re)
  call display(flow%velocity%vorticity_Z(),'flow',lim=0.5)
!
! -- Time update loop
  do while (flow%time<10*(2*amp/St).and..not.there)  ! run 10 cycles
     call geom%update(flow%time)
     call flow%update(geom)
! -- print force coefficients for each body in array
     write(9,'(f10.4,f8.4,80e16.8)') St*flow%time/(2*amp),flow%dt, &
     -2.*geom%bodies(2)%pforce(flow%pressure)/c
     flush(9)
     if(mod(flow%time,0.25*amp/St)<flow%dt) then
        if(root) print *,St*flow%time/(2*amp),flow%time/c,flow%dt
        call display(flow%velocity%vorticity_Z(),'flow',lim=0.5)
     end if
     inquire(file='.kill', exist=there)
  end do
  call flow%write()
  call mympi_end()
contains
!
!first foil motion
  real(8) pure function y1(t)
    real(8),intent(in) :: t
    y1 = amp*sin(pi*t*St1/amp)
  end function y1
  real(8) pure function p1(t)
    real(8),intent(in) :: t
    p1 = -pamp1*cos(pi*t*St1/amp)
  end function p1

real(8) pure function x(t)
    real(8),intent(in) :: t
    x=3*c
!+0.05*amp*cos(2*pi*t*St/amp)
end function x
! second foil motion
real(8) pure function y2(t)
    real(8),intent(in) :: t
    y2 = amp*sin(pi*t*St/amp + 1.75*pi)
  end function y2
real(8) pure function p2(t)
    real(8),intent(in) :: t
    p2 = -pamp2*cos(pi*t*St/amp + 1.75*pi)
  end function p2

  type(set) function naca(chord,thick,pivot,alpha)
    real,intent(in) :: chord
    real,intent(in),optional :: thick,pivot,alpha
    type(model_info) :: info
    ! the geometry is a NACA0012 defined from x=2.85-5.58, so
    real :: thick0=0.12, edge=2.8538, chord0=2.7303
    real :: t=0.12,piv=0.5,a=0

    if(present(thick)) t = thick
    if(present(pivot)) piv = pivot
    if(present(alpha)) a = alpha*pi/180. !convert to rad

    info%file = 'naca_square.IGS'
    info%x = (/-edge-piv*chord0,-10.271,-18.87/)
    info%r = (/a,0.,0./)
    info%s = chord/chord0*(/1.,t/thick0,-1./)
    info%xmax(1) = chord
    info%n = (/chord,chord,1./)
    ! surface_debug = .true.
    model_fill = .false.
    eps = 2.0
    naca = model_init(info)
  end function naca

end program flapping_NACA
