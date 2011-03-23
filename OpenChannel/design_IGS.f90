program test
  use nurbs
  use analytic
  implicit none
  character(20) :: string
  type(nurbsSet) :: nrb
  type(prop) :: pr
  integer :: i,n=41
  real(8) :: x(3) = 0

  call getarg(1,string)
  string = trim(string)
  print *,'7=',string
  open(7,file=string)
  call getarg(2,string)
  string = trim(string)
  print *,'8=',string
  open(8,file=string)
  open(9,file="iges.dat")
  call init_nurbs(nrb,7,8)
  call nurbs_plot(nrb,9)

  open(10,file="iges2.dat")
  write(10,*)'VARIABLES=x,y,z,dis,u,v,w'
  write(10,1) n
  x(3) = 2.0
  do i=1,n
     x(1) = (i-1.D0)/(n-1.D0)*3.D0-1.5D0
     pr = nrb.at.x
     write(10,2) x,pr%distance,pr%normal
  end do

  stop 'clean exit'
911 stop 'missing input IGES file'
1 format('ZONE, I =',i5,', F=POINT')
2 format(7e14.6)
end program test
