program test
  use nurbs
  use analytic
  implicit none
  character(20) :: string
  type(nurbsSet) :: nrb
  logical :: there

  call getarg(1,string)
  string = trim(string)
  inquire(file=string,exist=there)
  if(.not.there) stop 'missing 7'
  print *,'7=',string
  open(7,file=string)

  call getarg(2,string)
  string = trim(string)
  inquire(file=string,exist=there)
  if(.not.there) stop 'missing 8'
  print *,'8=',string
  open(8,file=string)

  open(9,file="srf.dat")
  open(10,file="grd.dat")
  call init_nurbs(nrb,7,8,.true.)
  call nurbs_plot(nrb,9,fgrd=10)
  stop 'clean exit'
end program test
