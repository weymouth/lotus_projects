!-------------------------------------------------------!
!------------------- Cylinder Array  -------------------!
!-------------------------------------------------------!
program array
  use fluidMod, only: fluid
  use bodyMod,  only: bodyUnion
  use mympiMod, only: init_mympi,mympi_end,mympi_rank
  use gridMod,  only: xg,composite
  use imageMod, only: display
  use geom_shape
  implicit none

  integer,parameter :: rows = 2, ndims = 2, finish = 300
  real,parameter    :: D = 32, DG = 7*D, Re_D = 100, nu = D/Re_D
  logical,parameter :: p(3) = (/.false.,.false.,.true./)

  type(fluid)      :: flow
  type(bodyUnion)  :: geom
  integer     :: n(3),i,b(3) = (/4,4,1/)
  real        :: A, pforce(40,3), vforce(40,3)
  logical     :: root, there=.false.
!
! -- Initialize simulation
  call init_mympi(ndims,set_blocks=b(:ndims),set_periodic=p(:ndims))
  root = mympi_rank()==0

  n = composite((/6.0*DG,2.4*DG,0./), prnt=root)
  call xg(1)%stretch(n(1), -10*DG, -0.6*DG, 0.6*DG, 20*DG, h_max=DG/D, prnt=root)
  call xg(2)%stretch(n(2), -10*DG, -0.6*DG, 0.6*DG, 10*DG, prnt=root)

  geom = make_array()
  if(root) print *, '-- made array --'
  A = D*n(3)*xg(3)%h

  call flow%init(n/b, geom, V=(/1.,0./), nu=nu)
  if(root) print *, '-- init complete --'
!
! -- Time update loop
  do while (flow%time/DG<finish .and. .not.there)
     call flow%update()
!
! -- print force coefficients for each body in array
     do i=1,geom%n
       pforce(i,:) = -2.*geom%bodies(i)%pforce(flow%pressure)/A
       vforce(i,:) = 2.*nu*geom%bodies(i)%vforce(flow%velocity)/A
     end do
     write(9,'(f10.4,f8.4,80e16.8)') flow%time/DG,flow%dt, &
                pforce(:geom%n,:2),vforce(:geom%n,:2)
     flush(9)
     if(mod(flow%time,D)<flow%dt) then
        call display(flow%velocity%vorticity_Z(),'vort',lim=0.25)
        if(root) print '(f10.4,f8.4)',flow%time/DG,flow%dt
      end if
     inquire(file='.kill', exist=there)
  end do
  call flow%write()
  call mympi_end
contains
!
! -- make an array of cylinders with group diameter DG
  type(BodyUnion) function make_array()
!    integer,parameter :: n(4) = (/6,13,19,25/)
    integer,parameter :: n(4) = (/6,12,18,26/) ! edit for symmetry
    real,parameter    :: R = 0.5*DG-0.5*D
    real    :: theta,xc,yc
    integer :: i,j
    call make_array%add(place_cyl(0.,0.))
    do j=1,rows
       do i=1,n(j)
          theta = 2.*pi*(i-1.)/real(n(j))
          xc = R*sin(theta)*real(j)/rows
          yc = R*cos(theta)*real(j)/rows
          call make_array%add(place_cyl(xc,yc))
       end do
    end do
  end function make_array
  type(cylinder) function place_cyl(xc,yc)
    real,intent(in) :: xc,yc
    place_cyl = cylinder(axis=3,radius=0.5*D,center=(/xc,yc,0./))
  end function place_cyl
end program array
