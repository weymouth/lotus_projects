!-------------------------------------------------------!
!---------------------- Ring thruster ------------------!
!-------------------------------------------------------!
program thruster
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  use geom_shape  ! to define geom (set,eps,plane, etc)
  implicit none
  real,parameter     :: f = 2              ! scaling factor
  real,parameter     :: D = 320/f          ! length scale
  real,parameter     :: dd = D/10          ! length scale
  real,parameter     :: xi = -10           ! rotation rate
  real(8),parameter  :: x = 0.45*D         ! x location
  real,parameter     :: Re = 200           ! Reynolds number
!
  integer,parameter  :: ndims = 2                       ! dimensions
  logical,parameter  :: p(2) = .false.                  ! periodic BCs
  real,parameter     :: nu = D/Re                       ! viscosity
  real,parameter     :: Ufric = sqrt(0.013/Re**(1./7.)) ! friction est.
  real,parameter     :: omega = 2.*xi/dd                ! spin freq
  integer            :: b(2) = (/4,4/)                  ! blocks
  integer            :: n(3)
  real               :: t1,dt,dtPrint=0.05,force(3)
!
  type(fluid)        :: flow
  type(set)          :: ring
  type(body)         :: bodies
!
! -- Initialize MPI (if MPI is ON)
#if MPION
  call init_mympi(ndims,set_blocks=b(:ndims),set_periodic=p(:ndims))
#else
  b=1
#endif
!
! -- Print run info
  if(mympi_rank()==0) print *, '-- Ring thurster test --'
  if(mympi_rank()==0) print '("   D=",f0.4,", nu=",f0.4,", y+=",f0.4)',D,nu,Ufric/nu
!
! -- Initialize array size
  n(:2) = composite(D*(/6.4,3.2/)/b); n(3) = 1
!
! -- Initialize and print grid
  call xg(1)%init(n(1)*b(1),0.6*D,4.1*D,1.0,f=f,r=1.02,d=4.)
  call xg(2)%init(n(2)*b(2),1.0*D,1.0*D,1.0,f=f,r=1.02)
  if(mympi_rank()==0) then
     call xg(1)%write(D)
     call xg(2)%write(D)
     print '("   total points=",i0)', product(n(1:2)*b)
  end if
!!$  call mympi_end()
!!$  stop
!
! -- Initialize the tandem geometry
  ring =  ((cylinder(1,1,3,dd/2.,0.,0.,0.).map.init_rigid(6,zip, spin).map.(init_affn()+(/x, x,0.D0/))) &
       .or.(cylinder(1,1,3,dd/2.,0.,0.,0.).map.init_rigid(6,zip,nspin).map.(init_affn()+(/x,-x,0.D0/))))
  bodies = ring.or.cylinder(1,1,3,D/2.,0.,0.,0.)
!
! -- Initialize fluid
  call flow%init(n,bodies,V=(/1.,0.,0./),nu=nu)
  call flow%write(bodies)
  call bodies%update(0.)
  if(mympi_rank()==0) print *, '-- init complete --'
!
! -- Time update loop
  do while (flow%time/D<3)
!
! -- update body and fluid
     dt = flow%dt
     t1 = flow%time+dt
     call flow%update(bodies)
!
! -- print force
     force = bodies%pforce(flow%pressure); force(3) = 0
     write(9,'(f10.4,f8.4,3e16.8)') t1/D,dt,2.*force/D
     flush(9)
!
! -- full output
     if(mod(t1,dtPrint*D)<dt) then
        call flow%write
        if(mympi_rank()==0) print 1,t1/D,2.*force/D
1       format("   t=",f0.4," force=",3f8.4)
     end if
  end do
  if(mympi_rank()==0) write(6,*) '--- complete --- '
  !
! -- Finalize MPI
#if MPION
  call mympi_end
#endif
contains
!
! -- motion definitions
  real(8) pure function spin(ts)
    implicit none
    real(8),intent(in) :: ts
    spin = omega
  end function spin
  real(8) pure function nspin(ts)
    implicit none
    real(8),intent(in) :: ts
    nspin = -omega
  end function nspin
  real(8) pure function zip(ts)
    implicit none
    real(8),intent(in) :: ts
    zip = 0.
  end function zip
end program thruster
