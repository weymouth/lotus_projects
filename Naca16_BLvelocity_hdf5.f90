program flapping_NACA
  use fluidMod, only: fluid
  use bodyMod,  only: body
  use imageMod, only: display
  use fieldMod, only: field
  use vectorMod
  use geom_shape
  use gridMod,  only: xg,composite
  use mympiMod
  use ioMod, only:write_vtk
  use iohdfMod
  use MGsolverMod, only: write_stats
  implicit none
!
! -- Define parameter, declare variables
  type(fluid)        :: flow
  type(body)         :: geom
  real,parameter     :: c = 700.,  D=0.16, St=0.30, k=0.6,cycles=30 !y+~5
  real,parameter     :: Re = 1.e4, nu=c/Re, amp=(St/k)/2.*c
  integer,parameter  :: ndims = 3
  real,parameter     :: Freq=k/c, amp2=0.25*c, pamp = asin(amp2/(0.75*c))
  real               :: f(3) = c*(/3.,1.2,0.1/), force(3), moment(3),power,vForce(3),Area
  integer            :: b(3) = (/32,13,1/), n(3)
  !
  integer,parameter  :: phases=32, cols = 12
  integer,parameter  :: cols13=phases+2, rows13=832, Prows = 500 ! Prows: pressure points along x/c
  integer,parameter  :: cols14=phases+2, rows14=24128, nnorm=80, ntan=100,Vrows = product((/3,nnorm,ntan/))
                        ! nnorm:velocity points along the height of BL, ntan:along x/c
                        ! Requirement: ntan << Prows , rows** > data's Prows and/or Vrows,
                        !       & mod(rows**,product(b)==0)
  integer            :: rows,row13,row14,n_ph
  real               :: flowtime0
  real               :: phase,n_count(phases),x_c(Prows),x_c_v(Vrows),BLdis(Vrows)
  real               :: p_probe(Prows),p_avrg(Prows,phases),p_in(phases+2,rows13)
  real               :: v_probe(3,nnorm,ntan),v_avrg(Vrows,phases),v_in(phases+2,rows14)
  logical            :: there = .false., root, p(3) = (/.FALSE.,.FALSE.,.TRUE./),continued=.false. !continue iteration, only for f_ave & v_ave
  logical            :: startwr = .true.
  type(field)        :: ave_p0,distF
  type(vfield)       :: ave_v0,distN
  integer(kind=8)    :: data9_id,data13_id,data14_id
  real,allocatable   :: data9(:,:),data13(:,:),data14(:,:)
  character(1)       :: n_i
  character(LEN=*),parameter :: folder_name='..//' !folder name of prev simulation between '//' only for f_ave & v_ave
  character(LEN=80)  :: prefolder
!
! -- Define MPI
#if MPION
  call init_mympi(ndims,set_blocks=b(1:ndims),set_periodic=p(1:ndims))
#else
  b=1
#endif
  root = mympi_rank()==0
!
! -- Initialize domain size and placement
  if(ndims==2) f(3)=1
  n = composite(f,prnt=root)
  call xg(1)%stretch(n(1),-3.*c,-.4*c,3.2*c,11.*c,h_min=2.,h_max=9.,prnt=root)
  call xg(2)%stretch(n(2),-3.*c,-.41999*c,.41999*c,3.*c,h_min=1.,prnt=root)
  if(ndims==3) xg(3)%h = 2.
  Area = c*xg(3)%h*n(3)

  if(root) print *, 'St=',St, ' k=', Freq*c, ' A/c(heave)=', amp/c,' A/c(pitch)=', amp2/c,  ' pamp(deg)=', pamp/pi*180
!
! -- Initialize the geometry (NACA0016)
  geom = naca(c,thick=D , pivot=0.25).map.init_rigid(6,pitch).map.init_rigid(2,heave)
!
! -- Initialize fluid
  call flow%init(n/b,geom,V=(/1.,0.,0./),nu=c/Re)
  flowtime0 = flow%time+flow%dt
  flow%time = 0.
!
!-- if continued simulation, read data from prev hf file & dump them into p_in array
 call getarg(1,prefolder)
! array initialization
 n_count=0.
 x_c=0.
 p_in=0.
 v_in=0.
 v_probe=0.
 p_probe=0.
 v_avrg = 0.
 p_avrg = 0.
 !! to uncomment for continuing simulation & reading from previous data
 ! if(len(trim(prefolder)) .ne. 0) then
 !   call read_dataset("fort13","Result",p_in)
 !   call read_dataset("fort14","Result",v_in)
 !   do n_ph=1,phases
 !     call dyn_average(p_in(n_ph+2,2:Prows+1),p_avrg,n_ph,n_init=int(p_in(n_ph+2,1)))
 !   !   !if(root) print *,size(v_in,1), size(v_in,2)
 !     call dyn_average_v(reshape(v_in(n_ph+2,2:Vrows+1),(/3,nnorm,ntan/)),v_avrg,n_ph,n_init=int(v_in(n_ph+2,1)))
 !     n_count(n_ph)=p_in(n_ph+2,1)
 !   end do
 ! end if
!
! -- other initialization
  call distN%init(n/b)
  call distN%applyBC()
  call distF%init(n/b)
  call distF%applyBC()
!
! allocate save data, rowlimit is from hdf5mod
  allocate(data9(rowlimit,cols),data13(rows13,cols13),data14(rows14,cols14))
  call init_hdf5array(rows,data9,data9_id,cols,rowlimit,'fort9')
  call init_hdf5array(row13,data13,data13_id,cols13,rows13,'fort13')
  call init_hdf5array(row14,data14,data14_id,cols14,rows14,'fort14')
!
! --continued parameter
  if(continued.eqv..true.) then
	call continue_sims(folder_name)
	! --re_init flow
  	call flow%init(n/b,geom,V=(/1.,0.,0./),nu=c/Re,exit=.true.)
  	flow%time = 0
  end if
!
  if(root) print *, '-- init complete --'
!
! -- Time update loop
  flow%time=flowtime0+flow%time
  do while (flow%time*Freq<cycles .and..not.there)
    call geom%update(flow%time+flow%dt)
    call flow%update(geom,write8=.false.)
!
! -- write to fort.9
    force = -geom%pforce(flow%pressure)/(0.5*Area)
    moment = -geom%pmoment(flow%pressure)/(0.5*Area*c)
    power = -geom%ppower(flow%pressure)/(0.5*Area)
    vForce = -nu*geom%vforce(flow%velocity)/(0.5*Area)
    !write(9,'(f10.4,f8.4,3e16.8,3e16.8,3e16.8,3e16.8)') flow%time*freq,flow%dt,force,moment,power,vForce
    !flush(9)
! -----------write data to allocated processes in hdf5 style -----!!
    data9(rows,1)    = flow%time*freq
    data9(rows,2)    = flow%dt
    data9(rows,3:5)  = force
    data9(rows,6:8)  = moment
    data9(rows,9)    = power
    data9(rows,10:12) = vforce
    if ((rows==rowlimit) .or. (flow%time*Freq.ge.cycles)) then
      if (startwr .eqv. .true.) then
          call write_data_hdf5(data9_id,data9,cols,rwlimit=rowlimit,in_parallel=.false.)
      else
          call write_data_hdf5(data9_id,data9,cols,rwlimit=rowlimit,extent = .true.,in_parallel=.false.)
      end if
      ! then erase data & re-init column number
      rows = 0
      data9 = 0.
      startwr = .false.
    end if
    rows=rows+1
    !
    if(mod(flow%time,0.03125/freq)<flow%dt) then
      if(mod(flow%time,0.25/freq)<flow%dt) then
        if(root) print *,flow%time*Freq,flow%time/c,flow%dt
        ! write last cycles every 0.03125
        if(flow%time*freq/cycles>(cycles-1.0)/cycles) then
          call flow%write(lambda=.true.,geom=geom,average=.true.)
        end if
      end if
      if(flow%time*freq/cycles>0.03125/cycles) then !start averaging only after cycle 1/32
         phase=floor((flow%time*freq-floor(flow%time*freq))*real(phases))+1.
         n_count(int(phase))=n_count(int(phase))+1
         distN = geom%fnormals()
         call geom%distance_field(distF)
         call distF%applyBC()
         if (ndims==2) then
           call pcontour(geom,distN,flow%velocity,flow%ub,distF,flow%pressure,int(phase), &
              p_probe,v_probe,x_of_c=x_c,x_of_c_v=x_c_v,thickness=BLdis)
         else !in 3D average all fields into 2D
           call pcontour(geom,distN%average(),flow%velocity%average(),flow%ub%average(),&
              distF%average(),flow%pressure%average(),int(phase),p_probe,v_probe,&
              x_of_c=x_c,x_of_c_v=x_c_v,thickness=BLdis)
         end if
         call dyn_average(p_probe,p_avrg,int(phase))
         call dyn_average_v(v_probe,v_avrg,int(phase)) ! v_avrg=[flattened(dim,thick,x/c),ph]
      end if
    end if
    inquire(file='.kill', exist=there)
  end do
  !
  call write_array(p_avrg,n_count,data13_id,data13,cols13,rows13,col1=x_c)
  call write_array(v_avrg,n_count,data14_id,data14,cols14,rows14,col1=x_c_v,col2=BLdis) ! the v_probe array is flattened
  !call flow%write(lambda=.true.,geom=geom,average=.true.)

  !
  ! ------------ closing hdf5---------!!
  call write_stats(ended=.true.)
  deallocate(data9,data13,data14)
  call close_data(data9_id)
  call close_data(data13_id)
  call close_data(data14_id)
  call close_file(hdfile_id)
  !
  call mympi_end()
contains
  subroutine pcontour(geom,Normals,flowV,bodyV,Dists,pnds,ph,pout,vout,x_of_c,x_of_c_v,thickness)
    use geom_global, only: eps
    use mympiMod, only: mympi_sum
    implicit none
    type(body), intent(in)  :: geom
    type(vfield), intent(in) :: Normals, flowV, bodyV
    type(field), intent(in) :: Dists,pnds !pnds now just flow%pressure instead of dot(p,kernel)
    integer,intent(in) :: ph
    real,intent(inout) :: pout(:),vout(:,:,:)
    real,intent(inout),optional :: x_of_c(:),x_of_c_v(:),thickness(:)
    real :: cosT,sinT,Hv,xe,xs,gap,x01,xfoil(3),xfoil2(3),xdyn(3),xdyn2(3)
    real :: xeps(3),xeps2(3),press,x_temp(3),xcheck(3)
    real :: dn(3),dist0,norml0(3),xeps_t(3),velo(3),vtan(3)
    integer :: d,i,j,k,n_iter,vthick,v_iter,l
    if(.not.geom%upToDate) return
    n_iter=size(pout)
    v_iter=size(vout,3) ! size in x/c
    vthick=size(vout,2) ! size in thickness
    !
    sinT = sin(pitch(real(flow%time+flow%dt,8)))
    cosT = cos(pitch(real(flow%time+flow%dt,8)))
    Hv = heave(real(flow%time,8))
    !
    xs = -0.25*c-2.*eps
    xe = 0.75*c+2.*eps
    gap = abs(xe-xs)/(n_iter-1.)
    !
    xfoil = 0.
    xdyn = 0.
    xfoil2 = 0.
    xdyn2 = 0.
    xeps = 0.
    xeps2 = 0.
    vtan = 0.
    vout = 0.
    k = 0
    l = 1
    !
    do i=1,n_iter
      x01 = (gap/(c+4.*eps))* float(i-1) !uniform sequence array [0,1]
      !
      !xfoil(i,1) is for guessing foil's points xfoil(i,2), xfoil2=lower foil
      xfoil(1) = gap* float(i-1) + xs ! this is x*c
      xfoil(2) = yfoil(x01)*(c+4.*eps)
      xfoil2(1) = xfoil(1)
      xfoil2(2) = -1.* yfoil(x01)*(c+4.*eps)
      !
      ! x_dynamic: transfer xfoil to pitch then heave
      xdyn(1) = xfoil(1)*cosT-xfoil(2)*sinT
      xdyn(2) = xfoil(1)*sinT+xfoil(2)*cosT+Hv
      xdyn2(1) = xfoil2(1)*cosT-xfoil2(2)*sinT
      xdyn2(2) = xfoil2(1)*sinT+xfoil2(2)*cosT+Hv
      !
      dist0 = Dists%at(xdyn) ! distance
      norml0(1) = Normals%e(1)%at(xdyn) ! Normals
      norml0(2) = Normals%e(2)%at(xdyn) ! Normals
      norml0(3) = 0.
      !
      ! coordinates at 1.5eps (3pts) away of body & dynamic
      xeps(1) = xdyn(1) + (1.5*eps - dist0) * norml0(1)
      xeps(2) = xdyn(2) + (1.5*eps - dist0) * norml0(2)
      xeps(3) = xdyn(3) + (1.5*eps - dist0) * norml0(3)

      press = pnds%at(xeps) !pressure @ eps at upper foil only
      pout(i)= press
      if(present(x_of_c)) x_of_c(i) = x01 ! normalised to [0,1]
      !
      ! velocity loop
      if(mod(i,n_iter/v_iter)==0) then ! points along the chord only, less than pressure points
        ! tangential direction of Normals pointing to trailing edge
        k = k+1
        dn(1)=norml0(2)*norml0(2)
        dn(2)=-norml0(1)*norml0(2)
        dn(3)=0.
        dn = dn/norm2(dn)
        if(norm2(dn)==0.) dn=0.
        !
        do j=1,vthick
          ! points in tangential height, measured above/= 1.5 eps, iterated every 0.5eps
          xeps_t(1) = xeps(1) + (real(j-1)*eps/2.) * norml0(1)
          xeps_t(2) = xeps(2) + (real(j-1)*eps/2.) * norml0(2)
          xeps_t(3) = xeps(3) + (real(j-1)*eps/2.) * norml0(3)
          !
          ! subtract(fluid velocity , body velocity)
          velo(1) = flowV%e(1)%at(xeps_t)-bodyV%e(1)%at(xeps_t)
          velo(2) = flowV%e(2)%at(xeps_t)-bodyV%e(2)%at(xeps_t)
          velo(3) = 0.
          !
          vtan=(/DOT_PRODUCT(velo,dn),DOT_PRODUCT(velo,norml0),0./)
          do d = 1,ndims
            !-- vout(d,j,k) vector velocity dims (d) at height (j) at position x/c (k-iter)
            vout(d,j,k) = vtan(d)
            if(present(x_of_c_v)) x_of_c_v(l) = x01
            if(present(thickness)) thickness(l) = 1.5*eps+real(j-1.)*eps/2. ! less accurate than dis%at(xeps_t)
            l=l+1
          end do
          if (ndims == 2) then
            vout(3,j,k) = 0.
            if(present(x_of_c_v)) x_of_c_v(l) = x01
            if(present(thickness)) thickness(l) = 1.5*eps+real(j-1.)*eps/2. ! less accurate than dis%at(xeps_t)
            l=l+1
          end if
        end do
      end if
    end do
  end subroutine pcontour
  !
  real function yfoil(x) result(y)
    real, intent(in) :: x
    y= 5.*D*(0.2969*(x**0.5)-0.1260*x-0.3516*(x**2)+0.2843*(x**3)-0.1015*(x**4))
  end function yfoil
  !
  ! moving average of pressure on different phases. This routine can only handle 1 variable averaging.
  subroutine dyn_average(pin,pinout,ph,n_init)
    implicit none
    real,intent(in)    :: pin(:)
    real,intent(inout) :: pinout(:,:)
    integer,intent(in) :: ph
    integer,intent(in),optional :: n_init
    real,allocatable,save :: n0(:),p_av(:,:)
    logical,save :: init_n0 = .true.
    real,allocatable :: temp_p(:,:)
    if(init_n0 .eqv. .true.) then
       allocate(n0(size(pinout,2)),p_av(size(pinout,1),size(pinout,2)))
       n0 = 0.
       p_av = 0.
       init_n0 = .false.
    end if
    allocate(temp_p(size(pinout,1),size(pinout,2)))
    if(n0(ph)==0.) then
      p_av(:,ph) = pin
    else
      p_av(:,ph)=p_av(:,ph)*(n0(ph)/(n0(ph)+1))
      temp_p(:,ph) = pin
      p_av(:,ph)=p_av(:,ph)+temp_p(:,ph)/(n0(ph)+1)
    end if
    pinout=p_av
    if(present(n_init)) then
      n0(ph) = n0(ph)+n_init
    else
      n0(ph) = n0(ph)+1
    end if
    deallocate(temp_p)
  end subroutine dyn_average
  !
  ! moving average for velocity fields on different phases
  subroutine dyn_average_v(pin,pinout,ph,n_init)
    implicit none
    real,intent(in)    :: pin(:,:,:)
    real,intent(inout) :: pinout(:,:) ! array of [flattened pin, ph]
    integer,intent(in) :: ph
    integer,intent(in),optional :: n_init
    real,allocatable,save :: nv0(:), v_av(:,:)
    integer :: i,j,k,ijk(3)
    logical,save :: init_nv0 = .true.
    real,allocatable :: temp_p(:,:),flat_pin(:)
    !
    !flattening input array
    allocate(flat_pin(size(pin)))
    flat_pin= pack(pin,.true.)
    !
    if(init_nv0 .eqv. .true.) then
       allocate(nv0(size(pinout,2)))
       allocate(v_av(size(pinout,1),size(pinout,2)))
       nv0 = 0.
       v_av = 0.
       init_nv0 = .false.
    end if
    allocate(temp_p(size(pinout,1),size(pinout,2)))
    if(nv0(ph)==0) then
      v_av(:,ph) = flat_pin
    else
      v_av(:,ph)=v_av(:,ph)*(nv0(ph)/(nv0(ph)+1))
      temp_p(:,ph) = flat_pin
      v_av(:,ph)=v_av(:,ph)+temp_p(:,ph)/(nv0(ph)+1)
    end if
    pinout=v_av
    if(present(n_init)) then
      nv0(ph) = nv0(ph)+n_init
    else
      nv0(ph) = nv0(ph)+1
    end if
    deallocate(temp_p)
  end subroutine dyn_average_v
  !
  ! write pressure to hdf5
  subroutine write_array(a_in,icount,data_in_id,data_in,col_in,row_in,col1,col2)
    real,intent(inout) :: a_in(:,:), icount(:)
    real,intent(inout),optional:: col1(:),col2(:)
    real,intent(inout) :: data_in(:,:)
    integer(kind=8),intent(inout) :: data_in_id
    integer,intent(in) :: col_in,row_in
    integer :: i,j,k,l,ij(2)
    ! data_in[row,column] as appear in hdfview
    ij(:)=(/size(a_in,1),size(a_in,2)/)
    data_in(1,3:) = (/ (icount(k), k=1,size(icount)) /) ! first row is number of iteration
    do i=1,ij(1)
      if (present(col1)) data_in(i+1,1) = col1(i) ! first column is x/c[0,1]
      if (present(col2)) data_in(i+1,2) = col2(i) ! first column is x/c[0,1]
      do j=1,ij(2)
          data_in(i+1,j+2) = a_in(i,j)
      end do
    end do
    call write_data_hdf5(data_in_id,data_in,col_in,rwlimit=row_in,in_parallel=.false.)
  end subroutine write_array
  !
  subroutine init_hdf5array(row_n,data_n,data_n_id,col_n,rowlim_n,arrayname)
    integer,intent(inout) :: row_n
    real,intent(inout)    :: data_n(:,:)
    integer(kind=8),intent(inout) :: data_n_id
    integer, intent(in) :: col_n, rowlim_n
    character(len=*),intent(in) :: arrayname
    row_n=1
    data_n=0.
    call create_data_hdf5(hdfile_id,data_n_id,col_n,group_id=sim_group_id, &
        dataset_name=arrayname,rwlimit=rowlim_n)
  end subroutine init_hdf5array
  !
  ! Motion
  real(8) pure function heave(t)
    real(8),intent(in) :: t
    heave = amp*sin(2*pi*Freq*t)
  end function heave
  real(8) pure function pitch(t)
    real(8),intent(in) :: t
    pitch = -pamp*cos(2*pi*Freq*t)
  end function pitch
  !
  !NACA geometry subroutine
  type(set) function naca(chord,thick,pivot,alpha)
    real,intent(in) :: chord
    real,intent(in),optional :: thick,pivot,alpha
    type(model_info) :: info
    ! the geometry is a NACA0012 defined from x=2.85-5.58, so
    real :: thick0=0.12, edge=2.8538, chord0=2.7303
    real :: t=0.12,piv=0.5,a=0

    if(present(thick)) t = thick
    if(present(pivot)) piv = pivot
    if(present(alpha)) a = alpha*pi/180. !convert to rad

    info%file = 'naca_square.IGS'
    info%x = (/-edge-piv*chord0,-10.271,-18.87/)
    info%r = (/a,0.,0./)
    info%s = chord/chord0*(/1.,t/thick0,-1./)
    !info%xmax(1) = chord
    !info%n = (/chord,chord,1./)
    ! surface_debug = .true.
    model_fill = .false.
    eps = 2.0
    naca = model_init(info)
  end function naca
  !
  ! -- field averaging, if needed
    subroutine f_ave(u_in,u_out,check,n_in,ph,nmin1)
      type(field),intent(in)             :: u_in
      type(field),save                   :: u_ave(4)
      type(field),intent(out)            :: u_out
      type(field)                        :: temp
      integer, save                      :: m(4) = 0
      logical,intent(in)                 :: check
      integer,intent(in)                 :: n_in(3),ph
      integer,intent(in),optional        :: nmin1
      if (check.eqv..false.) then
        ! initialize m and u_ave
        call u_ave(ph)%init(n_in/b)
        call u_ave(ph)%applyBC()
        m(ph)=0
      else
        ! m and u_ave are continued from last calling
        if (m(ph) == 0) then
    	    u_ave(ph) = u_in
        else
          u_ave(ph)%p = u_ave(ph)%p*(float(m(ph))/float(m(ph)+1))
  	      temp = u_in
  	      u_ave(ph)%p = u_ave(ph)%p+temp%p/float(m(ph)+1)
        end if
        u_out=u_ave(ph)
        if(present(nmin1)) then
          m(ph)=m(ph)+nmin1
        else
          m(ph)=m(ph)+1
        end if
      end if
      return
    end subroutine f_ave
  !
  ! -- vector averaging, if needed
    subroutine v_ave(u_in,u_out,check,n_in,ph,nmin1)
      type(vfield),intent(in)             :: u_in
      type(vfield),save                   :: u_ave(4)
      type(vfield),intent(out)            :: u_out
      type(vfield)                        :: temp
      integer, save                       :: m(4) = 0
      integer                             :: l
      logical,intent(in)                  :: check
      integer,intent(in)                  :: n_in(3),ph
      integer,intent(in),optional         :: nmin1
      !
      if (check.eqv..false.) then
        !initialize m and u_ave
        call u_ave(ph)%init(n_in/b)
        call u_ave(ph)%applyBC()
        m(ph)=0
      else
        ! m and u_ave are continued from last calling
        if (m(ph) == 0) then
          u_ave(ph) = u_in
        else
          call temp%init(n_in/b)
          temp=u_in
          do l=1,u_in%ndims
              u_ave(ph)%e(l)%p = u_ave(ph)%e(l)%p*(float(m(ph))/float(m(ph)+1))
              u_ave(ph)%e(l)%p = u_ave(ph)%e(l)%p+temp%e(l)%p/float(m(ph)+1)
          end do
        end if
        u_out=u_ave(ph)
        if(present(nmin1)) then
          m(ph)=m(ph)+nmin1
        else
          m(ph)=m(ph)+1
        end if
      end if
      return
    end subroutine v_ave
  !
   subroutine continue_sims(fldr) ! An example for calling f_ave & v_ave when cuntinuing simulation
      use ioMod, only:read_vtk,write_vtk
      type(field)        :: ave_p, ave_p0
      type(vfield)       :: ave_v, ave_v0
      character(LEN=*),intent(in):: fldr
      integer            :: i
      character(1)       :: ni
      do i=1,4 ! only for 4 phases/cycle
        	call ave_v%init(n/b)
        	call ave_p%init(n/b)
        	write(ni,"(i1)") i
        	call read_vtk(ave_v,ave_p,flow%time,name_in='Ave_'//ni//'.',folder_in=fldr)
        	call ave_v%applyBC()
        	call ave_p%applyBC()
        	call f_ave(ave_p,u_out=ave_p0,check=.true.,n_in=n,ph=i,nmin1=it(i))
        	call v_ave(ave_v,u_out=ave_v0,check=.true.,n_in=n,ph=i,nmin1=it(i))
      end do
    end subroutine continue_sims
end program flapping_NACA
