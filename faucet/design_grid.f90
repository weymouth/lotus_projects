!------------------------------------------------------------!
!                   Cartesian Grid Designer                  !
!                                                            !
! This program outputs the metrics for a multi-section       !
!  Cartesian grid. The output metrics in each direction are: !
!                                                            !
!   the starting point     - x(0)                            !
!   the number of sections - m                               !
!   for every section:                                       !
!     the number of points   - n                             !
!     the initial grid size  - a                             !
!     the expansion ratio    - r                             !
!                                                            !
!  such that x(i) = x(i-1)+a*r**(i-1)  for i=1,n             !
!  and each new section starts where the last left off.      !
!                                                            !
! This program actually outputs a fundamentally 6 section    !
!  grid. In sequence we have 1 uniform, 1 contraction, 2     !
!  uniform, 1 expansion, and 1 more uniform. The point in    !
!  between the two uniform sections defines the origin x=0.  !
!  However, sections can be empty (n=0) or combined.         !
! This 6 section grid is designed for external flows with    !
!  the object of interest centered at the origin where the   !
!  grid is finest and a coarser grid used in the far-field.  !
!  The two exterior uniform sections ensure the grid cells   !
!  do not exceed a given aspect ratio.                       !
!                                                            !
! The inputs are:                                            !
!                                                            !
!   f     - the scaling factor                               !
!   h     - the f=1 fine section grid size                   !
!   r     - the f=1 expansion ratio (should be around 1.01)  !
!   c     - the aspect ratio cutoff (should be around 4)     !
!   n     - the total number of points                       !
!   l1,l2 - the length of the fine sections                  !
!   alpha - grid length left of origin / right of origin     !
!                                                            !
!  where the scaling factor is used to make sets of grids    !
!  for grid size dependance studies. ie, a set with f=1,2,4  !
!  or f=1,sqrt{2},2. Given the value of f, the other metrics !
!  are scaled appropriately to achieve identical topology.   !
! The least intuitive scaled metric is b, the starting grid  !
!  size for the expansion sections. If, for f=1, we set      !
!  b = h*r' then the length of the section is                !
!                                                            !
!    l = h*r'*sum{r**(i-1),i=1,n} = h*r'*(1-r**n)/(1-r)      !
!                                                            !
!  To maintain exactly this length for all f we must use     !
!                                                            !
!    b = h*r'*(1-R)/(1-r)  => L = b*(1-R**N)/(1-R)           !
!                                                            !
!  where R=r**f is the new expansion ratio and N=n/f is the  !
!  new number of points in the section. Then R**N=r**n and   !
!  L=l as required. Setting r'= sqrt(r) ensures              !
!                                                            !
!    h*f < b < h*f*R                                         !
!                                                            !
!  for all f and r<1.                                        !
! Assuming that for a given problem that: n is a function of !
!  h and l, l1 l2 and alpha are relatively fixed by the      !
!  problem geometry, and r and c remain essentially fixed,   !
!  the only free parameters is h (and f which is used to     !
!  determine the dependance on h). As it should be.          !
!                                                            !
!------------------------------------------------------------!
program design_grid
  implicit none
  character(20) :: string
  real(8) :: f,h,l1,l2,a,b,r1,alpha,r,c
  integer :: n
!
! -- set up file
  open(7,file="inp.grid")
!
! -- read command line input
  call getarg(1,string)
  read(string,*,IOSTAT=n) f
  if(n.ne.0) then
     f = 1
  else
     call getarg(2,string)
     read(string,*,IOSTAT=n) h
     if(n==0) f = f/h
  end if
  print '("Grid size scaling, f =",f8.4)',f
!
! -- set global inputs
  r = 1.02   ! about 1.01
  c = 4.00   ! about 4  NOTE: c.eq.'cut midvalue' on coarsest grid
  r1 = r**f  ! -> R
  h = 1./128.
!
! -- set x-inputs
  n = 128*2
  l1 = 0
  l2 = 1.
  alpha = 0.  ! NOTE: alpha.le.'real alpha' on coarsest grid
!
! -- get metrics
  n = n/f !N
  a = h*f
  b = h*sqrt(r)*(1.-r1)/(1.-r)
  call make_print(n,l1,l2,r1,r1,b,a,b,alpha,c)
!
! -- same for y
  call make_print(n,l1,l2,r1,r1,b,a,b,alpha,c)
!
! -- set z-inputs
  n = 128*7
  l1 = 3.0
  l2 = 0.75
  alpha = 1.5
!
! -- get metrics
  n = n/f !N
  a = h*f
  b = h*sqrt(r)*(1.-r1)/(1.-r)
  call make_print(n,l1,l2,r1,r1,b,a,b,alpha,c)

end program design_grid
!----------------------------------------------------------------!
!
! Given fixed n,ds,and r this routine distributes the points
!  between the sections to get as close to the requested lengths
!  as possible and prints the output grid metrics to fort.7
! Section notation here does not match above. Section 2a,2b are 
!  the fine uniform sections. 1 and 3 are the expansion sections
!  and 0 and 4 are the coarse uniform sections.
! The maximum number of expansion points is determined by solving
!
!        ds0 = cut*ds2 > ds1*r1**(n1-1)
!
!  for n1, and the equivalent equation for n3.
!
subroutine make_print(n,l2a,l2b,r1,r3,ds1,ds2,ds3,alpha_in,cut)
  implicit none
  real(8) :: l2a,l2b,r1,r3,ds1,ds2,ds3,alpha_in,cut
  real(8) :: l0,l1,l2,l3,l4,ds0,ds4,alpha,left,right,ds1m,ds3m
  integer :: n,n1m,n3m,n0,n1,n2,n3,n4,m,mm
  logical :: debug=.true.
!
! -- get section 2
  l2 = l2a+l2b     ! total length of section 2
  n2 = l2/ds2      ! the section must have an integer number of points
  l2 = ds2*n2/l2   ! which adjusts the requested length by this factor
  l2a = l2a*l2     ! pass this adjustment onto l2a
  l2b = l2b*l2     ! pass this adjustment onto l2b
  l2 = l2a+l2b     ! new total length of section 2
!
! -- maximum number of expansion points
  ds0 = ds2*cut
  ds4 = ds2*cut
  n1m = log(ds0/ds1)/log(r1)+1
  n3m = log(ds4/ds3)/log(r3)+1
!
! -- get n1,n3
  n0=0; n1=0; n3=0; n4=0
  l0=0; l1=0; l3=0; l4=0; 
  left = l2a; right = l2b
  do m=0,n-n2-1
     n1 = min(n1m,m)
     n0 = m-n1
     n3 = min(n-n2-n1-n0,n3m)
     n4 = n-n2-n1-n0-n3
!
! -- get lengths
     l0 = ds0*n0
     l1 = ds1*(1.-r1**n1)/(1.-r1)
     l3 = ds3*(1.-r3**n3)/(1.-r3)
     l4 = ds4*n4
     left = l0+l1+l2a
     right = l2b+l3+l4
!
! -- balance
     alpha = left/right
     if(alpha.ge.alpha_in) exit
  end do
!
! -- final expansion sizes
  ds1m = ds1*r1**(n1-1)
  ds3m = ds3*r3**(n3-1)
!
! -- print segment coeffs to screen
  print 1,0, n0,ds0,l0
  print 1,1, n1,ds1,l1
  print 1,2, n2,ds2,l2
  print 1,3, n3,ds3,l3
  print 1,4, n4,ds4,l4
!
! -- factors
  mm = log(dble(n))/log(2.D0)
  do m=1,mm
     if(mod(n,2**m).ne.0) exit
  end do
  print '("  total points",i4)', n
  print '("  factors, 2^",i2,"*",i2)', m-1,n/2**(m-1)
1 format('seg',i2,', n=',i4,', ds=',f8.4,', l=',f8.4)
!
! -- test stuff
  if(debug) then
     print '("  maximums, n1= ",i4,", n3= ",i4)', n1m,n3m
     print '("  cut midvalue",2f8.4)', 0.5*(ds1*r1**(n1m-1)/ds2+ds1*r1**(n1m)/ds2)
     print '("  total lengths",2f8.4)', left,right
     print '("  final alpha",f8.4)', alpha
  end if
!
! -- write to file
  write(7,'(i4,19X," ! total points")') n
  write(7,'(f10.6,13X," ! grid starting position")') -left
  m = count((/n0.ne.0,n1.ne.0,n2.ne.0,n3.ne.0,n4.ne.0/))
  write(7,'(i3,20X,   " ! number of segments")') m
  if(n0.ne.0) &
  write(7,'(i3,2f10.6," ! segment coeffs (n,a,r)")') n0,ds0,1.
  if(n1.ne.0) &
  write(7,'(i3,2f10.6," ! segment coeffs (n,a,r)")') n1,ds1m,1./r1
  write(7,'(i3,2f10.6," ! segment coeffs (n,a,r)")') n2,ds2,1.
  if(n3.ne.0) &
  write(7,'(i3,2f10.6," ! segment coeffs (n,a,r)")') n3,ds3,r3
  if(n4.ne.0) &
  write(7,'(i3,2f10.6," ! segment coeffs (n,a,r)")') n4,ds4,1.

end subroutine make_print
