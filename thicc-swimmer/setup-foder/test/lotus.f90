program swimming_plate
  !
    use bodyMod,    only: body
    use fluidMod,   only: fluid
    use mympiMod,   only: init_mympi,mympi_end,mympi_rank
    use gridMod,    only: xg,composite
    use imageMod,   only: display
    use geom_shape
    implicit none
  !
  ! -- Physical parameters
    real,parameter     :: Re = 37947.33192202055
  !
    real,parameter     :: L=256.0, nu=L/Re
    real, parameter    :: finish=3
    integer            :: b(3) = [2,4,1]
  !
  ! -- Hyperparameters
    real, parameter    :: thicc=0.03*L
    real, parameter    :: A = 0.1*L, St_d = 0.3, k_x=16.5, k_z=16.0, h_roughness=0.01
    real, parameter    :: a_coeff = 0.28, &
                          b_coeff = -0.13, &
                          c_coeff = 0.05, &
                          k_coeff = 0.7101387960861528, &
                          f = St_d/(2.*A)
  !
  ! -- Dimensions
    integer            :: n(3), ndims=3
  !
  ! -- Setup solver
    logical            :: there = .false., root, p(3) = [.FALSE.,.FALSE.,.TRUE.]
    real               :: m(3), z
    type(fluid)        :: flow
    type(body)         :: geom
    !
    ! -- Outputs
      real            :: dt, t, pforce(3), vforce(3), ppower(3), enstrophy(3)
    !
    ! -- Initialize
      call init_mympi(ndims,set_blocks=b(1:ndims),set_periodic=p(1:ndims))
      root = mympi_rank()==0
      if(root) print *,'Setting up the grid, body and fluid'
      if(root) print *,'-----------------------------------'
    !
    ! -- Setup the grid
      if(ndims==2) then
        z = 0.0
      else
        z = 0.125
      end if
      m = [1.5, 1.5, z]
      n = composite(L*m,prnt=root)
      call xg(1)%stretch(n(1), -2.0*L, -.5*L, 2.5*L, 7.0*L,  h_min=4., h_max=10., prnt=root)
      call xg(2)%stretch(n(2), -2.0*L, -0.9*L, 0.9*L, 2.0*L, h_min=2., prnt=root)
      if(ndims==3) xg(3)%h = 4.
    !
    ! -- Call the geometry and kinematics
  
      geom = wavy_wall(L, thicc).and.upper(L,thicc).and.lower(L,thicc).map.init_warp(2,h,doth,dh)
    !
    ! -- Initialise fluid
      call flow%init(n/b,geom,V=[1.,0.,0.],nu=nu,exit=.true.)
      ! call avrg%init(flow)
      flow%time = 0
    !
      if(root) print *,'Starting time update loop'
      if(root) print *,'-----------------------------------'
    !
      ! call flow%write(geom, write_vtr=.false.)
      time_loop: do while (flow%time<finish/f.and..not.there)
        t = flow%time
        dt = flow%dt
        call geom%update(t+dt) ! update geom
        call flow%update(geom)
  
        write(9,'(f10.4,f8.4,4e16.8,4e16.8,4e16.8,4e16.8,4e16.8,4e16.8,&
              4e16.8,4e16.8,4e16.8,4e16.8,4e16.8,4e16.8,4e16.8,4e16.8)')&
              t*f,dt,pforce,ppower,vforce,enstrophy
  
        pforce = 2.*geom%pforce(flow%pressure)/(A*n(3)*xg(3)%h)
        ppower = 2.*geom%ppower(flow%pressure)/(A*n(3)*xg(3)%h)
        vforce = 2.*nu*geom%vforce_f(flow%velocity)/(A*n(3)*xg(3)%h)
        enstrophy = flow%velocity%enstrophy(mean=.true.)/(n(1)*n(2)*n(3))
  
        if((mod(t,1./f)<dt).and.(root)) print "('Time:',f15.3,'. Thrust:',f15.3,3f12.6)",&
        t*f, (pforce(1)+vforce(1))
  
        inquire(file='../.kill', exist=there)
        if (there) exit time_loop
        if((t>(finish-1)/f).and.(mod(t,0.1/f)<dt)) call flow%write(geom, write_vtr=.false.)
      end do time_loop

    if(root) print *,'Loop complete: writing restart files and exiting'
    if(root) print *,'-----------------------------------'
    call mympi_end
  contains
  !
  type(set) function wavy_wall(length, thickness) result(geom)
    real,intent(in) :: length, thickness
    geom = plane([1.,-0.,0.],[length,0.,0.]) & ! end cap
    .and.plane([-1.,0.,0.],[0.,0.,0.]) ! front cap
  end function
  !
  type(set) function upper(length, thickness) result(geom)
    real,intent(in) :: length, thickness
    geom = plane([0.,1.,0.],[0.,0.,0.])&
    .map.init_warp(2,egg_top,dotegg_top,degg_top)&
    .map.init_warp(2,naca_warp,dotnaca_warp,dnaca_warp)
  end function
  !
  type(set) function lower(length, thickness) result(geom)
    real,intent(in) :: length, thickness
    geom = plane([-0.,-1.,0.],[0.,0.,0.])&
    .map.init_warp(2,egg_bottom,dotegg_bottom,degg_bottom)&
    .map.init_warp(2,naca_warp_neg,dotnaca_warp_neg,dnaca_warp_neg)
  end function
  !
  ! -- General kinematics
  real pure function h(x)
    real,intent(in) :: x(3)
    h = amp(x(1))*sin(arg(x(1)))
  end function h
  real pure function doth(x)
    real,intent(in) :: x(3)
    doth = amp(x(1))*cos(arg(x(1)))*2*pi*f
  end function doth
  pure function dh(x)
    real,intent(in) :: x(3)
    real            :: dh(3)
    dh = 0
    dh(1) = damp(x(1))*sin(arg(x(1))) &
          - amp(x(1))*cos(arg(x(1)))*2*pi*k_coeff/L
  end function dh
  real pure function arg(x)
    real,intent(in) :: x
    real :: xp
    xp = min(max(x/L,0.),1.)
    arg = (2*pi*(f*flow%time - k_coeff*xp))
  end function arg
  real pure function amp(x)
    real,intent(in) :: x
    real :: xp
    xp = min(max(x/L,0.),1.)
    amp = A*(((a_coeff*(xp**2))+(b_coeff*(xp))+(c_coeff))/(a_coeff+b_coeff+c_coeff))
  end function amp
  real pure function damp(x)
    real,intent(in) :: x
    real :: xp
    xp = min(max(x/L,0.),1.)
    damp = A*(b_coeff+2.*(a_coeff*xp))/(L*(a_coeff+b_coeff+c_coeff))
  end function damp
  !
  ! -- Egg carton roughness distribution
  real pure function egg_top(x)
    real,intent(in) :: x(3)
    egg_top = (h_roughness*L)*sin((k_x)*(2*pi*x(1)/L)-pi/2)*cos((k_z)*(2*pi*x(3)/L))
  end function egg_top
  pure function degg_top(x)
    real,intent(in) :: x(3)
    real            :: degg_top(3)
    degg_top = 0
    degg_top(1) = (h_roughness*L)*cos((k_x)*(2*pi*x(1)/L)-pi/2)&
                                  *sin((k_z)*(2*pi*x(3)/L))*(k_x)*(2*pi/L)
    degg_top(3) = -(h_roughness*L)*cos((k_x)*(2*pi*x(1)/L-pi/2))&
                                  *sin((k_z)*(2*pi*x(3)/L))*(k_z)*(2*pi/L)
  end function degg_top
  real pure function dotegg_top(x)
    real,intent(in) :: x(3)
    dotegg_top = 0
  end function dotegg_top
  !
  real pure function egg_bottom(x)
    real,intent(in) :: x(3)
    egg_bottom = (h_roughness*L)*sin((k_x)*(2*pi*x(1)/L)-3*pi/2)&
                                 *cos((k_z)*(2*pi*x(3)/L))
  end function egg_bottom
  pure function degg_bottom(x)
    real,intent(in) :: x(3)
    real            :: degg_bottom(3)
    degg_bottom = 0
    degg_bottom(1) = (h_roughness*L)*cos((k_x)*(2*pi*x(1)/L)-3*pi/2)&
                                     *sin((k_z)*(2*pi*x(3)/L))*(k_x)*(2*pi/L)
    degg_bottom(3) = -(h_roughness*L)*cos((k_x)*(2*pi*x(1)/L)-3*pi/2)&
                                      *sin((k_z)*(2*pi*x(3)/L))*(k_z)*(2*pi/L)
  end function degg_bottom
  real pure function dotegg_bottom(x)
    real,intent(in) :: x(3)
    dotegg_bottom = 0
  end function dotegg_bottom
  !
  ! -- Create NACA shape by warping based on B-spline equation
  real pure function naca_warp(x)
    real,intent(in) :: x(3)
    real :: xp
    xp = min(max(x(1)/L,0.),1.)
    naca_warp = 5*0.12*L*(0.2969*sqrt(xp)-0.1260*(xp)&
                -0.3516*(xp)**2+0.2843*(xp)**3-0.1015*(xp)**4)
  end function naca_warp
  pure function dnaca_warp(x)
    real,intent(in) :: x(3)
    real            :: dnaca_warp(3)
    real :: xp
    xp = min(max(x(1)/L,0.),1.)
    dnaca_warp = 0
    ! dnaca_warp(1) = 0.6 * L * (-0.126 - 0.406 * xp**3 + 0.8529 * xp**2 - 0.7032 * xp + 0.14845/sqrt(xp))
  end function dnaca_warp
  real pure function dotnaca_warp(x)
    real,intent(in) :: x(3)
    dotnaca_warp = 0
  end function dotnaca_warp
  real pure function naca_warp_neg(x)
    real,intent(in) :: x(3)
    real :: xp
    xp = min(max(x(1)/L,0.),1.)
    naca_warp_neg = -5*0.12*L*(0.2969*sqrt(xp)-0.1260*(xp)&
                -0.3516*(xp)**2+0.2843*(xp)**3-0.1015*(xp)**4)
  end function naca_warp_neg
  pure function dnaca_warp_neg(x)
    real,intent(in) :: x(3)
    real            :: dnaca_warp_neg(3)
    real :: xp
    xp = min(max(x(1)/L,0.),1.)
    dnaca_warp_neg = 0
    ! dnaca_warp_neg(1) = -0.6 * L * (-0.126 - 0.406 * xp**3 + 0.8529 * xp**2 - 0.7032 * xp + 0.14845/sqrt(xp))
  end function dnaca_warp_neg
  real pure function dotnaca_warp_neg(x)
    real,intent(in) :: x(3)
    dotnaca_warp_neg = 0
  end function dotnaca_warp_neg
end program swimming_plate
  
