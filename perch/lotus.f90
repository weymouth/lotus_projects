!-------------------------------------------------------!
!---------------- 2D perching test case ----------------!
!-------------------------------------------------------!
program perch
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  use geom_shape, only: operator(.map.),init_rigid,pi
  use imageMod,   only: display
  implicit none
  real,parameter     :: L = 100            ! length scale
  real,parameter     :: Re = 2e3           ! Reynolds number
  real,parameter     :: Xi = 1./4.         ! Shape change number
!
  real,parameter     :: nu = L/Re          ! viscosity
  real,parameter     :: T = L/Xi           ! motion period
  integer            :: b(3) = (/2,1,1/)     ! blocks
  integer            :: n(3)
  real               :: tau,u
  logical            :: root
  type(fluid)        :: flow
  type(body)         :: foil
!
! -- Initialize MPI (if MPI is ON)
  call init_mympi(2,set_blocks=b)
  root = mympi_rank()==0
!
! -- Print run info
  if(root) print *, '-- Perch Test --'
  if(root) print '("   L=",f0.4,", nu=",f0.4)',L,nu
!
! -- Initialize domain size and placement
  n = composite(L*(/8.,8.,0./),prnt=root)
  xg%left = -n/2
!
! -- Initialize the foil geometry
  foil = naca(L,0.16).map.init_rigid(6,alpha,omega)
!
! -- Initialize fluid
  call flow%init(n/b,foil,V=(/1.,0.,0./),nu=nu)
  if(root) print *, '-- init complete --'
!
! -- Time update loop
  do while (flow%time < 6*L+T)
    tau = (flow%time-5.*L)/T      ! non-dimensional time
    u = min(1.,max(0.,1.-tau))    ! inflow velocity
!
!-- update and write foil and fluid
     if(tau+flow%dt/T > 0 .and. tau < 1+flow%dt/T ) &
       call foil%update(tau)
     call flow%update(foil,V=(/u,0.,0./))
     if(mod(flow%time,0.1*T)<flow%dt) then
        call flow%write(foil)
        call display(flow%velocity%vorticity_Z(),'out_vort',10./L)
     end if
!
! -- print force
     write(9,'(f10.4,f8.4,3e16.8)') &
        tau,flow%dt,2.*foil%pforce(flow%pressure)/L
     flush(9)
  end do
  if(root) write(6,*) '--- complete --- '
  call mympi_end()
contains
!
! -- geometry definitions
  type(set) function naca(c,o)
    use geom_shape, only: set,model_info,model_init,eps,assignment(=)
    real,intent(in) :: c,o
    type(model_info) :: info
    info%file = 'naca_square.IGS'               ! geometry file name
    info%x = (/-2.853*(1+o),-10.271,-18.876/)   ! geometry shift
    info%s = 0.36626*c*(/1,1,-1/)               ! geometry scale
    info%xmax(1) = c                            ! set size of geom box
    info%n = (/c,c,1./)                         ! set size of geom grid
    eps = 2.0                                   ! set BDIM width
    naca = model_init(info)                     ! initialize model
  end function naca
!
! -- motion definitions
  real(8) pure function alpha(tau)  ! rotation angle
    real(8),intent(in) :: tau
    if(tau<0) then
       alpha = 0.
    else if(tau>1) then
       alpha = pi/2.
    else
       alpha = pi/2.*(tau-sin(2.*pi*tau)/(2.*pi))
    end if
  end function alpha
  real(8) pure function omega(tau)  ! rotation velocity
    real(8),intent(in) :: tau
    omega = 0
    if(tau>0 .and. tau<1) omega = pi/2.*(1-cos(2.*pi*tau))/T
  end function omega
end program perch
