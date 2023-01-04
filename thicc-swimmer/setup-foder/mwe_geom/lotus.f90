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
    real,parameter     :: Re = 10000
  !
    real,parameter     :: L=256.0, nu=L/Re
    real, parameter    :: finish=3
    integer            :: b(3) = [2,4,1]
  !
  ! -- Hyperparameters
    real, parameter    :: thicc=0.03*L
    real, parameter    :: A = 0.1*L, St_d = 0.3, k_x=4.5, k_z=4.0, h_roughness=0.01
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
      z = 0.125
      m = [1.5, 1.5, z]
      n = composite(L*m,prnt=root)
      call xg(1)%stretch(n(1), -2.0*L, -.5*L, 2.5*L, 7.0*L,  h_min=4., h_max=10., prnt=root)
      call xg(2)%stretch(n(2), -2.0*L, -0.9*L, 0.9*L, 2.0*L, h_min=2., prnt=root)
      if(ndims==3) xg(3)%h = 4.
      geom = upper(L,thicc).and.lower(L,thicc).and.wavy_wall(L,thicc)
      !
    ! -- Initialise fluid
      call flow%init(n/b,geom,V=[1.,0.,0.],nu=nu,exit=.true.)
      flow%time = 0
      call flow%write(geom, write_vtr=.false.)
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
    .map.init_warp(2,naca_warp,dotnaca_warp,dnaca_warp)
  end function
  !
  type(set) function lower(length, thickness) result(geom)
    real,intent(in) :: length, thickness
    geom = plane([-0.,-1.,0.],[0.,-0.05*L,0.])&
    .map.init_warp(2,naca_warp_neg,dotnaca_warp_neg,dnaca_warp_neg)
  end function
  !
  ! -- Create NACA shape by warping
  real pure function naca_warp(x)
    real,intent(in) :: x(3)
    real :: xp, a, b
    a = 0.21273426557106268
    b = -0.20456259679086652
    xp = min(max(x(1)/L,0.),1.)
    naca_warp = L*(a * xp    &    
                  + b * xp**2&
                  )
    
  end function naca_warp
  pure function dnaca_warp(x)
    real,intent(in) :: x(3)
    real            :: dnaca_warp(3)
    real :: xp, a, b
    a = 0.21273426557106268
    b = -0.20456259679086652
    xp = min(max(x(1)/L,0.),1.)
    dnaca_warp = 0
    dnaca_warp(1) = (a             &    
                      + 2 * b * xp   &
                      )
  end function dnaca_warp
  real pure function dotnaca_warp(x)
    real,intent(in) :: x(3)
    dotnaca_warp = 0
  end function dotnaca_warp
  real pure function naca_warp_neg(x)
    real,intent(in) :: x(3)
    real :: xp, a, b
    a = 0.21273426557106268
    b = -0.20456259679086652
    xp = min(max(x(1)/L,0.),1.)
    naca_warp_neg = -L*(a * xp   &    
                  + b * xp**2   &
                  )
                
  end function naca_warp_neg
  pure function dnaca_warp_neg(x)
    real,intent(in) :: x(3)
    real            :: dnaca_warp_neg(3)
    real :: xp, a, b
    a = 0.21273426557106268
    b = -0.20456259679086652
    xp = min(max(x(1)/L,0.),1.)
    dnaca_warp_neg = 0
    dnaca_warp_neg(1) = -(a         &    
                      + 2 * b * xp   &
                      )
  
  end function dnaca_warp_neg
  real pure function dotnaca_warp_neg(x)
    real,intent(in) :: x(3)
    dotnaca_warp_neg = 0
  end function dotnaca_warp_neg
end program swimming_plate
  
