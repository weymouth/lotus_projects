!-------------------------------------------------------!
!------------------ Whisker motion ---------------------!
!-------------------------------------------------------!
program design_geom
  use analytic
  implicit none
  open(7,file='inp.geom')
!
! -- write pointer to the IGS,ISG2 files and make srf,grd files
  write(7,'(a4)') 'body'
  call set_write(7,.set.(/'inp.IGS ','inp.IGS2'/))
!
! -- free interface set
  write(7,'(a4)') 'fint'
  call set_write(7,.set.plane(1,(/0,0,1/),0,0,0))
  close(7)

end program design_geom
