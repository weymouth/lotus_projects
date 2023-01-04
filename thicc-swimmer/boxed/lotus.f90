program tke_test
  use fluidMod,   only: fluid
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  implicit none

  type(fluid)       :: flow
  real              :: Re=1.e5,tke,nu,tkebox
  integer,parameter :: L=2048,num=9,ndims=2
  integer           ::n(3),b(3)=[1,1,1]
  real              :: lbox(3)=[0.,0.,0.],hbox(3)=[1.,1.,0.]
!
! -- Set parameters
    nu = L/(Re)
    call init_mympi(ndims,set_blocks=b)
    n = composite(L*[1.,1.,0.],prnt=mympi_rank()==0)
    ! call xg(1)%stretch(n(1), -1.*L, -.25*L, 0.25*L, 1.*L, prnt=mympi_rank()==0)
    ! call xg(2)%stretch(n(2), -1.*L, -.25*L, 0.25*L, 1.*L, prnt=mympi_rank()==0)
    ! call xg(1)%stretch(n(1), 0., 0., 0.5*L, 2.0*L, prnt=mympi_rank()==0)
    ! call xg(2)%stretch(n(2), 0., 0., 0.5*L, 2.0*L, prnt=mympi_rank()==0)
    ! if(ndims==3) xg(3)%h = 4.
!
! -- Initialize the velocity field
    call flow%init(n/b,nu=nu)
    call flow%velocity%eval(tgv)
    ! call flow%reset_u0()
    tke = flow%velocity%tke()
    tkebox = flow%velocity%tke(mean=.true.,lowerb=lbox*L, upperb=hbox*L)

    if(mympi_rank()==0) write(*,*) '--- TKE Test --- '
    if(mympi_rank()==0) write(*,*) 'tke',tke,tkebox!,flow%velocity%tke(lowerb=lbox*L, upperb=hbox*L),flow%velocity%tke()
    call mympi_end()
  contains
    pure function tgv(x) result(v)
      real,intent(in) :: x(3)
      real :: v(3)
      v(1) = 1.
      v(2) = 0.
      v(3) = 0
    end function tgv

end program tke_test