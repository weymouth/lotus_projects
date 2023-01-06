program profile_test
    use fluidMod, only: Cyl_test, TGV_test
    implicit none
  
    real,parameter  :: D = 64, tmx= 10., tprnt =.25
    real, parameter :: Re = 1.e5
    integer         :: b(3) = (/2,2,2/), n_tgv(3)
  
    n_tgv = (/128,128,128/)
    call TGV_test(Re,int(D),n_tgv,b,tmx,tprnt,ldisplay=.false.)
    ! call Cyl_test(D,b,tmx,tprnt)
  
  end program profile_test