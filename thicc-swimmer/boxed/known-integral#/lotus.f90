program enstrophy_test
  use fluidMod,   only: fluid
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  implicit none

  type(fluid)       :: flow
  real              :: Re=1.e5,kn,tke,omega,nu,omegabox
  integer,parameter :: L=4096,num=9
  integer           ::n(3),b(3)=[2,2,1]
  real              :: lbox(3)=[0.1,0.1,1./L],hbox(3)=[1.,1.,1./L]
!
! -- Set parameters
    nu = L/(Re)
    call init_mympi(2,set_blocks=b)
    n = composite(L*[1.,1.,0.],prnt=mympi_rank()==0)
!
! -- Initialize the velocity field
    call flow%init(n/b,nu=nu)
    call flow%velocity%eval(tgv)
    call flow%reset_u0()
    omega = flow%velocity%tke(mean=.true.)
    omegabox = flow%velocity%tke(mean=.true.,lowerb=lbox*L, upperb=hbox*L)

    if(mympi_rank()==0) write(*,*) '--- Enstrophy Test --- '
    if(mympi_rank()==0) write(*,*) 'Enstrophy',omega,omegabox
    call mympi_end()
  contains
    pure function tgv(x) result(v)
      real,intent(in) :: x(3)
      real :: v(3)
      v(1) = 0
      v(2) = 1.
      v(3) = 0
    end function tgv

end program enstrophy_test