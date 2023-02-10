module neutraldata3Dobj_fclaw

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: stderr=>error_unit
use phys_consts, only: wp,debug,pi,Re
use inputdataobj, only: inputdata
use neutraldataobj, only: neutraldata
use neutraldata3Dobj, only: neutraldata3D
use meshobj, only: curvmesh
use gemini3d_config, only: gemini_cfg
use reader, only: get_simsize3,get_simsize2,get_grid2,get_precip
use timeutils, only: dateinc,date_filename
use grid, only: gridflag

implicit none (type, external)
private
public :: neutraldata3D_fclaw


!> type definition for 3D neutral data that will be provided from a parallel model (i.e. one that runs with GEMINI)
type, extends(neutraldata3D) :: neutraldata3D_fclaw
  ! these are for storing information about locations that are being communicated to the neutral model
  real(wp), dimension(:), pointer :: zlocsi,xlocsi,ylocsi
  integer, dimension(:,:), pointer :: ilocsi
  real(wp), dimension(:,:), pointer :: dataxyzinow    ! will need to be rotated prior to placing in final arrays

  contains
    ! procedures specific to this class
    procedure :: get_locationsi       ! get a list of interpolation sites that are in bounds with regards to the neutral modoel
    procedure :: set_datainow          ! place a set of interpolated data into the data array at indices corresponding to locations
                                      !   provided 

    ! overriding procedures
    procedure :: update
    procedure :: init_storage

    ! bindings for deferred procedures
    procedure :: init=>init_neu3D_fclaw
    procedure :: set_coordsi=>set_coordsi_neu3D_fclaw
    procedure :: load_data=>load_data_neu3D_fclaw
    procedure :: load_grid=>load_grid_neu3D_fclaw
    procedure :: load_size=>load_size_neu3D_fclaw

    ! destructor
    final :: destructor
end type neutraldata3D_fclaw


contains
  !> initialize arrays for storing object data once the sizes are set
    subroutine init_storage(self)
    class(neutraldata3D_fclaw), intent(inout) :: self
    integer :: lc1,lc2,lc3
    integer :: lc1i,lc2i,lc3i
    integer :: l0D
    integer :: l1Dax1,l1Dax2,l1Dax3
    integer :: l2Dax23,l2Dax12,l2Dax13
    integer :: l3D

    ! check sizes are set
    if (.not. self%flagsizes) error stop 'inpudata:init_storage(); must set sizes before allocations...'

    ! local size variables for convenience
    lc1=self%lc1; lc2=self%lc2; lc3=self%lc3;
    lc1i=self%lc1i; lc2i=self%lc2i; lc3i=self%lc3i;
    l0D=self%l0D
    l1Dax1=self%l1Dax1; l1Dax2=self%l1Dax2; l1Dax3=self%l1Dax3;
    l2Dax23=self%l2Dax23; l2Dax12=self%l2Dax12; l2Dax13=self%l2Dax13;
    l3D=self%l3D

    ! NOTE: type extensions are reponsible for zeroing out any arrays they will use...

!    ! input data coordinate arrays (presume plaid)
!    allocate(self%coord1(lc1),self%coord2(lc2),self%coord3(lc3))

    ! interpolation site arrays (note these are flat, i.e. rank 1), if one needed to save space by not allocating unused block
    !   could override this procedure...
    allocate(self%coord1i(lc1i*lc2i*lc3i),self%coord2i(lc1i*lc2i*lc3i),self%coord3i(lc1i*lc2i*lc3i))

    ! No singleton array objects being allocated by this extension; but this doesn't hurt anything so leave in place
!    ! coordinate sites for singleton axes depend on mangling of data
!    if (self%flagdipmesh) then    ! mangle 2,3 sizes
!      allocate(self%coord1iax1(lc1i),self%coord2iax2(lc3i),self%coord3iax3(lc2i))
!    else
!      allocate(self%coord1iax1(lc1i),self%coord2iax2(lc2i),self%coord3iax3(lc3i))
!    end if
!
!    allocate(self%coord2iax23(lc2i*lc3i),self%coord3iax23(lc2i*lc3i))
!    allocate(self%coord1iax13(lc1i*lc3i),self%coord3iax13(lc1i*lc3i))
!    allocate(self%coord1iax12(lc1i*lc2i),self%coord2iax12(lc1i*lc2i))

!    ! allocate object arrays for input data at a reference time.  FIXME: do we even need to store this perm. or can be local to
!    ! load_data?
!    allocate(self%data0D(l0D))
!    allocate(self%data1Dax1(lc1,l1Dax1), self%data1Dax2(lc2,l1Dax2), self%data1Dax3(lc3,l1Dax3))
!    allocate(self%data2Dax23(lc2,lc3,l2Dax23), self%data2Dax12(lc1,lc2,l2Dax12), self%data2Dax13(lc1,lc3,l2Dax13))
!    allocate(self%data3D(lc1,lc2,lc3,l3D))
!
!    ! allocate object arrays for interpolation sites at reference times
!    allocate(self%data0Di(l0D,2))
!    allocate(self%data1Dax1i(lc1i,l1Dax1,2), self%data1Dax2i(lc2i,l1Dax2,2), self%data1Dax3i(lc3i,l1Dax3,2))
!    allocate(self%data2Dax23i(lc2i,lc3i,l2Dax23,2), self%data2Dax12i(lc1i,lc2i,l2Dax12,2), self%data2Dax13i(lc1i,lc3i,l2Dax13,2))
!    allocate(self%data3Di(lc1i,lc2i,lc3i,l3D,2))

    ! allocate object arrays at interpolation sites for current time.  FIXME: do we even need to store permanently?
    allocate(self%data0Dinow(l0D))
    allocate(self%data1Dax1inow(lc1i,l1Dax1), self%data1Dax2inow(lc2i,l1Dax2), self%data1Dax3inow(lc3i,l1Dax3))
    allocate(self%data2Dax23inow(lc2i,lc3i,l2Dax23), self%data2Dax12inow(lc1i,lc2i,l2Dax12), self%data2Dax13inow(lc1i,lc3i,l2Dax13))
    allocate(self%data3Dinow(lc1i,lc2i,lc3i,l3D))

    allocate(self%coord1i(lc1i*lc2i*lc3i),self%coord2i(lc1i*lc2i*lc3i),self%coord3i(lc1i*lc2i*lc3i))   
    allocate(self%ximat(lc1i,lc2i,lc3i),self%yimat(lc1i,lc2i,lc3i),self%zimat(lc1i,lc2i,lc3i))
    allocate(self%proj_ezp_e1(lc1i,lc2i,lc3i),self%proj_ezp_e2(lc1i,lc2i,lc3i),self%proj_ezp_e3(lc1i,lc2i,lc3i))
    allocate(self%proj_eyp_e1(lc1i,lc2i,lc3i),self%proj_eyp_e2(lc1i,lc2i,lc3i),self%proj_eyp_e3(lc1i,lc2i,lc3i))
    allocate(self%proj_exp_e1(lc1i,lc2i,lc3i),self%proj_exp_e2(lc1i,lc2i,lc3i),self%proj_exp_e3(lc1i,lc2i,lc3i))

    self%flagalloc=.true.
  end subroutine init_storage


    !> set pointer variables to locations for storage of interpolated data (3D always for neutral input)
  subroutine setptrs_grid(self)
    class(neutraldata3D_fclaw), intent(inout) :: self

!    ! set aliases for prev data
!    self%dnOiprev=>self%data3Di(:,:,:,1,1)
!    self%dnN2iprev=>self%data3Di(:,:,:,2,1)
!    self%dnO2iprev=>self%data3Di(:,:,:,3,1)
!    self%dvn1iprev=>self%data3Di(:,:,:,4,1)
!    self%dvn2iprev=>self%data3Di(:,:,:,5,1)
!    self%dvn3iprev=>self%data3Di(:,:,:,6,1)
!    self%dTniprev=>self%data3Di(:,:,:,7,1)
!
!    ! set pointers for next data
!    self%dnOinext=>self%data3Di(:,:,:,1,2)
!    self%dnN2inext=>self%data3Di(:,:,:,2,2)
!    self%dnO2inext=>self%data3Di(:,:,:,3,2)
!    self%dvn1inext=>self%data3Di(:,:,:,4,2)
!    self%dvn2inext=>self%data3Di(:,:,:,5,2)
!    self%dvn3inext=>self%data3Di(:,:,:,6,2)
!    self%dTninext=>self%data3Di(:,:,:,7,2)

    ! set aliases for interpolated data that is "outward facing"
    self%dnOinow=>self%data3Dinow(:,:,:,1)
    self%dnN2inow=>self%data3Dinow(:,:,:,2)
    self%dnO2inow=>self%data3Dinow(:,:,:,3)
    self%dvn1inow=>self%data3Dinow(:,:,:,4)
    self%dvn2inow=>self%data3Dinow(:,:,:,5)
    self%dvn3inow=>self%data3Dinow(:,:,:,6)
    self%dTninow=>self%data3Dinow(:,:,:,7)
  end subroutine setptrs_grid


  !> initialize object for this type of neutral input data.  In this case the source arrays are not needed at all since
  !    we expect the top-level app to populate the neutral data for us.
  subroutine init_neu3D_fclaw(self,cfg,sourcedir,x,dtmodel,dtdata,ymd,UTsec)
    class(neutraldata3D_fclaw), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    character(*), intent(in) :: sourcedir               ! will not be used but part of call signature for uniformity, could replace
                                                        !   with "generic" procedure
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dtmodel,dtdata
    integer, dimension(3), intent(in) :: ymd            ! target date of initiation
    real(wp), intent(in) :: UTsec                       ! target time of initiation
    integer :: lc1,lc2,lc3
    character(:), allocatable :: strname    ! allow auto-allocate for strings

    ! force 3D interpolation regardless of working subarray size
    self%flagforcenative=.true.

    ! tell our object where its data are and give the dataset a name
    !   Note that we don't have a source location per se
    !call self%set_source(sourcedir)
    call self%set_source('')
    strname='neutral perturbations fclaw (3D)'
    call self%set_name(strname)
    call self%set_cadence(dtdata)
    ! I believe that this won't even be used but set to false anyway
    !self%flagdoinput=cfg%flagdneu/=0
    self%flagdoinput=.false.

    ! set sizes, we have 7 arrays all 3D (irrespective of 2D vs. 3D neutral input).  for 3D neutral input
    !    the situation is more complicated that for other datasets because you cannot compute the number of
    !    source grid points for each worker until you have root compute the entire grid and slice everything up
    allocate(self%lc1,self%lc2,self%lc3)                                     ! these are pointers, even though scalar
    self%lzn=>self%lc1; self%lxn=>self%lc2; self%lyn=>self%lc3;              ! these referenced while reading size and grid data
    call self%set_sizes( &
             0, &          ! number scalar parts to dataset
             0, 0, 0, &    ! number 1D data along each axis
             0, 0, 0, &    ! number 2D data
             7, &          ! number 3D datasets
             x)            ! The main purpose of this is to set the number of 3D datasets (other params already set)

    ! allocate space for arrays, for this use case we only need the "now" arrays
    call self%init_storage()

    ! define interpolation site coordinates
    call self%set_coordsi(cfg,x)

    ! we no longer need to load these data; we just populate sites on target grid
    !call self%load_sizeandgrid_neu3D(cfg)          ! cfg needed to form source neutral grid

!    ! set aliases to point to correct source data arrays
!    self%dnO=>self%data3D(:,:,:,1)
!    self%dnN2=>self%data3D(:,:,:,2)
!    self%dnO2=>self%data3D(:,:,:,3)
!    self%dvnz=>self%data3D(:,:,:,4)
!    self%dvnx=>self%data3D(:,:,:,5)
!    self%dvny=>self%data3D(:,:,:,6)
!    self%dTn=>self%data3D(:,:,:,7)

    ! call to base class procedure to set pointers for prev,now,next
    call self%setptrs_grid()

!    ! initialize previous data so we get a correct starting value
!    self%dnOiprev=0
!    self%dnN2iprev=0
!    self%dnO2iprev=0
!    self%dvn1iprev=0
!    self%dvn2iprev=0
!    self%dvn3iprev=0
!    self%dTniprev=0

    ! initialize previous data so we get a correct starting value
    self%dnOinow=0
    self%dnN2inow=0
    self%dnO2inow=0
    self%dvn1inow=0
    self%dvn2inow=0
    self%dvn3inow=0
    self%dTninow=0

    ! We don't really even have to keep up with datasets times since controlling app is going to give us whatever data
    !   we need and we dont' have to make decisions about when to load new data.  
    ! set to start time of simulation - not needed since assigned by update on first call.  FIXME: a bit messy
    !self%ymdref(:,1)=cfg%ymd0; self%ymdref(:,2)=cfg%ymd0;
    !self%UTsecref(1)=cfg%UTsec0; self%UTsecref(2)=cfg%UTsec0;

    ! No priming required
    !call self%prime_data(cfg,x,dtmodel,ymd,UTsec)
  end subroutine init_neu3D_fclaw



  !> set coordinates for target interpolation points; for neutral inputs we are forced to do some of the property array allocations here
  subroutine set_coordsi_neu3D_fclaw(self,cfg,x)
    class(neutraldata3D_fclaw), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp) :: theta1,phi1,theta2,phi2,gammarads,theta3,phi3,gamma1,gamma2,phip
    real(wp) :: xp,yp
    real(wp), dimension(3) :: ezp,eyp,tmpvec,exprm
    real(wp) :: tmpsca
    integer :: ix1,ix2,ix3,iyn,izn,ixn,iid

    ! Space for coordinate sites and projections in neutraldata3D object
    self%zi=>self%coord1i; self%xi=>self%coord2i; self%yi=>self%coord3i;     ! alias coordinates of interpolation sites

    !Neutral source locations specified in input file, here referenced by spherical magnetic coordinates.
    phi1=cfg%sourcemlon*pi/180
    theta1=pi/2 - cfg%sourcemlat*pi/180

    !Convert plasma simulation grid locations to z,rho values to be used in interoplation.  altitude ~ zi; lat/lon --> rhoi.  Also compute unit vectors and projections
!    if (mpi_cfg%myid==0) then
!      print *, 'Computing alt,radial distance values for plasma grid and completing rotations'
!    end if

    self%zimat=x%alt     !vertical coordinate is just altitude array already stored in grid object
    do ix3=1,x%lx3
      do ix2=1,x%lx2
        do ix1=1,x%lx1
          ! interpolation based on geomag
          theta2=x%theta(ix1,ix2,ix3)                    !field point zenith angle

          !print*, ' center NS set',shape(self%zi),shape(self%zimat),theta2,x%theta(ix1,ix2,ix3)

          if (x%lx2/=1) then
            phi2=x%phi(ix1,ix2,ix3)                      !field point azimuth, full 3D calculation
          else
            phi2=phi1                                    !assume the longitude is the samem as the source in 2D, i.e. assume the source epicenter is in the meridian of the grid
          end if

          !we need a phi locationi (not spherical phi, but azimuth angle from epicenter), as well, but not for interpolation - just for doing vector rotations
          theta3=theta2
          phi3=phi1
          gamma1=cos(theta2)*cos(theta3)+sin(theta2)*sin(theta3)*cos(phi2-phi3)
          if (gamma1 > 1) then     !handles weird precision issues in 2D
            gamma1 = 1
          else if (gamma1 < -1) then
            gamma1 = -1
          end if
          gamma1=acos(gamma1)

          gamma2=cos(theta1)*cos(theta3)+sin(theta1)*sin(theta3)*cos(phi1-phi3)
          if (gamma2 > 1) then     !handles weird precision issues in 2D
            gamma2= 1
          else if (gamma2 < -1) then
            gamma2= -1
          end if
          gamma2=acos(gamma2)
          xp=Re*gamma1
          yp=Re*gamma2     !this will likely always be positive, since we are using center of earth as our origin, so this should be interpreted as distance as opposed to displacement

          ! coordinates from distances
          if (theta3>theta1) then       !place distances in correct quadrant, here field point (theta3=theta2) is is SOUTHward of source point (theta1), whreas yp is distance northward so throw in a negative sign
            yp= -yp            !do we want an abs here to be safe
          end if
          if (phi2<phi3) then     !assume we aren't doing a global grid otherwise need to check for wrapping, here field point (phi2) less than source point (phi3=phi1)
            xp= -xp
          end if

          self%ximat(ix1,ix2,ix3)=xp     !eastward distance
          self%yimat(ix1,ix2,ix3)=yp     !northward distance

          !PROJECTIONS FROM NEUTURAL GRID VECTORS TO PLASMA GRID VECTORS
          !projection factors for mapping from axisymmetric to dipole (go ahead and compute projections so we don't have to do it repeatedly as sim runs
          ezp=x%er(ix1,ix2,ix3,:)

          tmpvec=ezp*x%e2(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_ezp_e2(ix1,ix2,ix3)=tmpsca

          tmpvec=ezp*x%e1(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_ezp_e1(ix1,ix2,ix3)=tmpsca

          tmpvec=ezp*x%e3(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)    !should be zero, but leave it general for now
          self%proj_ezp_e3(ix1,ix2,ix3)=tmpsca

          eyp= -x%etheta(ix1,ix2,ix3,:)

          tmpvec=eyp*x%e1(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_eyp_e1(ix1,ix2,ix3)=tmpsca

          tmpvec=eyp*x%e2(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_eyp_e2(ix1,ix2,ix3)=tmpsca

          tmpvec=eyp*x%e3(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_eyp_e3(ix1,ix2,ix3)=tmpsca

          exprm=x%ephi(ix1,ix2,ix3,:)   !for 3D interpolation need to have a unit vector/projection onto x-direction (longitude)

          tmpvec=exprm*x%e1(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_exp_e1(ix1,ix2,ix3)=tmpsca

          tmpvec=exprm*x%e2(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_exp_e2(ix1,ix2,ix3)=tmpsca

          tmpvec=exprm*x%e3(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_exp_e3(ix1,ix2,ix3)=tmpsca
        end do
      end do
    end do

    !Assign values for flat lists of grid points
!    if (mpi_cfg%myid==0) then
!      print*, '...Packing interpolation target points...'
!    end if
    self%zi=pack(self%zimat,.true.)     !create a flat list of grid points to be used by interpolation functions
    self%yi=pack(self%yimat,.true.)
    self%xi=pack(self%ximat,.true.)

    ! FIXME: do we need to have the new grid code clear its unit vectors?  Or maybe this isn't a huge waste of memory???
!    if (mpi_cfg%myid==0) then
!      print*, '...Clearing out unit vectors (after projections)...'
!    end if
    !call clear_unitvecs(x)

!    if(mpi_cfg%myid==0) then
!      print*, 'Interpolation coords:  ',minval(self%zi),maxval(self%zi), &
!                                        minval(self%xi),maxval(self%xi), &
!                                        minval(self%yi),maxval(self%yi)
!      print*, 'Projection checking:  ',minval(self%proj_exp_e1),maxval(self%proj_exp_e1), &
!                                       minval(self%proj_exp_e2),maxval(self%proj_exp_e2), &
!                                       minval(self%proj_exp_e3),maxval(self%proj_exp_e3)
!    end if

    self%flagcoordsi=.true.
  end subroutine set_coordsi_neu3D_fclaw


  !> do nothing stub - type extensions must override this to perform whatever load steps are needed for their data types
  subroutine load_size_neu3D_fclaw(self)
    class(neutraldata3D_fclaw), intent(inout) :: self

    return
  end subroutine load_size_neu3D_fclaw


  !> do nothing stub
  subroutine load_grid_neu3D_fclaw(self)
    class(neutraldata3D_fclaw), intent(inout) :: self

    return
  end subroutine load_grid_neu3D_fclaw


  subroutine load_data_neu3D_fclaw(self,t,dtmodel,ymdtmp,UTsectmp)
    class(neutraldata3D_fclaw), intent(inout) :: self
    real(wp), intent(in) :: t,dtmodel
    integer, dimension(3), intent(inout) :: ymdtmp
    real(wp), intent(inout) :: UTsectmp
    integer :: iid
    integer :: lhorzn                        !number of horizontal grid points
    real(wp), dimension(:,:,:), allocatable :: paramall
    character(:), allocatable :: fn

    ! this should not be a no-op, i.e. no automatic updating because we expect external controlling app to do this  

!    lhorzn=self%lyn
!    ymdtmp = self%ymdref(:,2)
!    UTsectmp = self%UTsecref(2)
!    call dateinc(self%dt,ymdtmp,UTsectmp)                !get the date for "next" params
!
!    !read in the data from file
!    if (mpi_cfg%myid==0) then    !root
!      !in the 3D case we cannot afford to send full grid data and need to instead use neutral subgrid splits defined earlier
!      allocate(paramall(self%lzn,self%lxnall,self%lynall))     ! space to store a single neutral input parameter
!
!      !print*, '  date and time (neutral3D):  ',ymdtmp,UTsectmp
!
!      fn = date_filename(self%sourcedir,ymdtmp,UTsectmp) // ".h5"
!
!      if (debug) print *, 'READ neutral 3D data from file: ',fn
!      call hf%open(fn, action='r')
!
!      call hf%read('/dn0all', paramall)
!      if (.not. all(ieee_is_finite(paramall))) error stop 'dnOall: non-finite value(s)'
!      if (debug) print*, 'Min/max values for dnOall:  ',minval(paramall),maxval(paramall)
!      call dneu_root2workers(paramall,tag%dnO,self%slabsizes,self%indx,self%dnO)
!      call hf%read('/dnN2all', paramall)
!      if (.not. all(ieee_is_finite(paramall))) error stop 'dnN2all: non-finite value(s)'
!      if (debug) print*, 'Min/max values for dnN2all:  ',minval(paramall),maxval(paramall)
!      call dneu_root2workers(paramall,tag%dnN2,self%slabsizes,self%indx,self%dnN2)
!      call hf%read('/dnO2all', paramall)
!      if (.not. all(ieee_is_finite(paramall))) error stop 'dnO2all: non-finite value(s)'
!      if (debug) print*, 'Min/max values for dnO2all:  ',minval(paramall),maxval(paramall)
!      call dneu_root2workers(paramall,tag%dnO2,self%slabsizes,self%indx,self%dnO2)
!      call hf%read('/dTnall', paramall)
!      if (.not. all(ieee_is_finite(paramall))) error stop 'dTnall: non-finite value(s)'
!      if (debug) print*, 'Min/max values for dTnall:  ',minval(paramall),maxval(paramall)
!      call dneu_root2workers(paramall,tag%dTn,self%slabsizes,self%indx,self%dTn)
!      call hf%read('/dvnrhoall', paramall)
!      if (.not. all(ieee_is_finite(paramall))) error stop 'dvnrhoall: non-finite value(s)'
!      if (debug) print*, 'Min/max values for dvnrhoall:  ',minval(paramall),maxval(paramall)
!      call dneu_root2workers(paramall,tag%dvnrho,self%slabsizes,self%indx,self%dvny)
!      call hf%read('/dvnzall', paramall)
!      if (.not. all(ieee_is_finite(paramall))) error stop 'dvnzall: non-finite value(s)'
!      if (debug) print*, 'Min/max values for dvnzall:  ',minval(paramall),maxval(paramall)
!      call dneu_root2workers(paramall,tag%dvnz,self%slabsizes,self%indx,self%dvnz)
!      call hf%read('/dvnxall', paramall)
!      if (.not. all(ieee_is_finite(paramall))) error stop 'dvnxall: non-finite value(s)'
!      if (debug) print*, 'Min/max values for dvnxall:  ',minval(paramall),maxval(paramall)
!      call dneu_root2workers(paramall,tag%dvnx,self%slabsizes,self%indx,self%dvnx)
!
!      call hf%close()
!      deallocate(paramall)
!    else     !workers
!      !receive a subgrid copy of the data from root
!      call dneu_workers_from_root(tag%dnO,self%dnO)
!      call dneu_workers_from_root(tag%dnN2,self%dnN2)
!      call dneu_workers_from_root(tag%dnO2,self%dnO2)
!      call dneu_workers_from_root(tag%dTn,self%dTn)
!      call dneu_workers_from_root(tag%dvnrho,self%dvny)
!      call dneu_workers_from_root(tag%dvnz,self%dvnz)
!      call dneu_workers_from_root(tag%dvnx,self%dvnx)
!    end if
!
!
!    if (mpi_cfg%myid==mpi_cfg%lid/2 .and. debug) then
!      print *, 'Min/max values for dnO:  ',mpi_cfg%myid,minval(self%dnO),maxval(self%dnO)
!      print *, 'Min/max values for dnN:  ',mpi_cfg%myid,minval(self%dnN2),maxval(self%dnN2)
!      print *, 'Min/max values for dnO2:  ',mpi_cfg%myid,minval(self%dnO2),maxval(self%dnO2)
!      print *, 'Min/max values for dvnx:  ',mpi_cfg%myid,minval(self%dvnx),maxval(self%dvnx)
!      print *, 'Min/max values for dvnrho:  ',mpi_cfg%myid,minval(self%dvny),maxval(self%dvny)
!      print *, 'Min/max values for dvnz:  ',mpi_cfg%myid,minval(self%dvnz),maxval(self%dvnz)
!      print *, 'Min/max values for dTn:  ',mpi_cfg%myid,minval(self%dTn),maxval(self%dTn)
!    !  print*, 'coordinate ranges:  ',minval(zn),maxval(zn),minval(rhon),maxval(rhon),minval(zi),maxval(zi),minval(rhoi),maxval(rhoi)
!    end if
  end subroutine load_data_neu3D_fclaw


  !> overriding procedure for updating neutral atmos (need additional rotation steps)
  subroutine update(self,cfg,dtmodel,t,x,ymd,UTsec)
    class(neutraldata3D_fclaw), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    real(wp), intent(in) :: dtmodel             ! need both model and input data time stepping
    real(wp), intent(in) :: t                   ! simulation absoluate time for which perturabation is to be computed
    class(curvmesh), intent(in) :: x            ! mesh object
    integer, dimension(3), intent(in) :: ymd    ! date for which we wish to calculate perturbations
    real(wp), intent(in) :: UTsec               ! UT seconds for which we with to compute perturbations

    ! this should not be a no-op, i.e. no automatic updating because we expect external controlling app to do this

!    ! execute a basic update
!    call self%update_simple(cfg,dtmodel,t,x,ymd,UTsec)
!
!    ! FIXME: more efficient to rotate the winds only when interpolations are done...
!    ! now we need to rotate velocity fields following interpolation (they are magnetic ENU prior to this step)
!    call self%rotate_winds()
!
!    ! FIXME: check if we need to further rotate these winds into geographic coordinates
!
!    if (mpi_cfg%myid==mpi_cfg%lid/2 .and. debug) then
!      print*, ''
!      print*, 'neutral data size:  ',mpi_cfg%myid,self%lzn,self%lxn,self%lyn
!      print*, 'neutral data time:  ',ymd,UTsec
!      print*, ''
!      print *, 'Min/max values for dnOinext:  ',mpi_cfg%myid,minval(self%dnOinext),maxval(self%dnOinext)
!      print *, 'Min/max values for dnNinext:  ',mpi_cfg%myid,minval(self%dnN2inext),maxval(self%dnN2inext)
!      print *, 'Min/max values for dnO2inext:  ',mpi_cfg%myid,minval(self%dnO2inext),maxval(self%dnO2inext)
!      print *, 'Min/max values for dvn1inext:  ',mpi_cfg%myid,minval(self%dvn1inext),maxval(self%dvn1inext)
!      print *, 'Min/max values for dvn2inext:  ',mpi_cfg%myid,minval(self%dvn2inext),maxval(self%dvn2inext)
!      print *, 'Min/max values for dvn3inext:  ',mpi_cfg%myid,minval(self%dvn3inext),maxval(self%dvn3inext)
!      print *, 'Min/max values for dTninext:  ',mpi_cfg%myid,minval(self%dTninext),maxval(self%dTninext)
!      print*, ''
!      print *, 'Min/max values for dnOinow:  ',mpi_cfg%myid,minval(self%dnOinow),maxval(self%dnOinow)
!      print *, 'Min/max values for dnNinow:  ',mpi_cfg%myid,minval(self%dnN2inow),maxval(self%dnN2inow)
!      print *, 'Min/max values for dnO2inow:  ',mpi_cfg%myid,minval(self%dnO2inow),maxval(self%dnO2inow)
!      print *, 'Min/max values for dvn1inow:  ',mpi_cfg%myid,minval(self%dvn1inow),maxval(self%dvn1inow)
!      print *, 'Min/max values for dvn2inow:  ',mpi_cfg%myid,minval(self%dvn2inow),maxval(self%dvn2inow)
!      print *, 'Min/max values for dvn3inow:  ',mpi_cfg%myid,minval(self%dvn3inow),maxval(self%dvn3inow)
!      print *, 'Min/max values for dTninow:  ',mpi_cfg%myid,minval(self%dTninow),maxval(self%dTninow)
!    end if
  end subroutine update


  !> Find and return a list of interpolation sites that are in bounds for the external neutral model grid.  This
  !    procedure needs to allocate space to store the (unknown upon entry) number of locations needed.  
  subroutine get_locationsi(self,zlims,xlims,ylims)
    class(neutraldata3D_fclaw), intent(inout) :: self
    real(wp), dimension(2), intent(in) :: zlims,xlims,ylims    ! global boundary of neutral grid we are accepting data from
    integer :: ix1,ix2,ix3,lx1,lx2,lx3
    integer :: ipts,lpts,itarg

    lx1=self%lc1i; lx2=self%lc2i; lx3=self%lc3i;

    ! get rid of any data that might be sitting in our output array
    if (associated(self%zlocsi)) deallocate(self%zlocsi)
    if (associated(self%xlocsi)) deallocate(self%xlocsi)
    if (associated(self%ylocsi)) deallocate(self%ylocsi)
    if (associated(self%ilocsi)) deallocate(self%ilocsi)
    if (associated(self%dataxyzinow)) deallocate(self%dataxyzinow)

    ! count the number of in bounds points so we can do allocation
    lpts=0
    do ipts=1,lx1*lx2*lx3
      if (self%xi(ipts) > xlims(1) .and. self%xi(ipts) < xlims(2) .and. &
            self%yi(ipts) > ylims(1) .and. self%yi(ipts) < ylims(2) .and. &
            self%zi(ipts) > zlims(1) .and. self%zi(ipts) < zlims(2)) then
        lpts=lpts+1
      end if
    end do
    allocate(self%zlocsi(lpts))
    allocate(self%xlocsi,self%ylocsi, mold=self%zlocsi)
    allocate(self%ilocsi(lpts,3))
    allocate(self%dataxyzinow(lpts,7))

    ! now make another pass through the data to copy out the locations
    itarg=1
    do ipts=1,lx1*lx2*lx3
      if (self%xi(ipts) > xlims(1) .and. self%xi(ipts) < xlims(2) .and. &
            self%yi(ipts) > ylims(1) .and. self%yi(ipts) < ylims(2) .and. &
            self%zi(ipts) > zlims(1) .and. self%zi(ipts) < zlims(2)) then
        self%zlocsi(itarg)=self%zi(ipts)
        self%xlocsi(itarg)=self%xi(ipts)
        self%ylocsi(itarg)=self%yi(ipts)
        self%ilocsi(itarg,:)=[ix1,ix2,ix3]
        itarg=itarg+1
        if (itarg > lpts) return    ! we are done and can stop iterating through the list of points
      end if
    end do
  end subroutine get_locationsi


  ! Populate object arrays with information from an external model corresponding to locations defined by a call to
  !   get_locationsi
  subroutine set_datainow(self)
    class(neutraldata3D_fclaw), intent(inout) :: self
    integer lpts,ipts,ix1,ix2,ix3

    if (.not. associated(self%dataxyzinow)) then
      error stop 'neutraldata3D_fclaw:  attempting to assign unallocated space to inputdata'
    end if

    ! number of points being placed in object
    lpts=size(self%zlocsi)

    ! points not being specified by input model need to be zeroed out
    self%dnOinow=0
    self%dnN2inow=0
    self%dnO2inow=0
    self%dvn1inow=0
    self%dvn2inow=0
    self%dvn3inow=0
    self%dTninow=0

    ! place data into object
    do ipts=1,lpts
      ix1=self%ilocsi(ipts,1)
      ix2=self%ilocsi(ipts,2)
      ix3=self%ilocsi(ipts,3)
      self%dnOinow(ix1,ix2,ix3)=self%dataxyzinow(ipts,1)
      self%dnN2inow(ix1,ix2,ix3)=self%dataxyzinow(ipts,2)
      self%dnO2inow(ix1,ix2,ix3)=self%dataxyzinow(ipts,3)
      self%dvn1inow(ix1,ix2,ix3)=self%dataxyzinow(ipts,4)
      self%dvn2inow(ix1,ix2,ix3)=self%dataxyzinow(ipts,5)
      self%dvn3inow(ix1,ix2,ix3)=self%dataxyzinow(ipts,6)
      self%dTninow(ix1,ix2,ix3)=self%dataxyzinow(ipts,7)
    end do

    ! insure winds are correctly rotated before returning
    call self%rotate_winds()

    ! deallocate the temp space for the data exchange now that we are done populating
    deallocate(self%zlocsi,self%xlocsi,self%ylocsi,self%ilocsi,self%dataxyzinow)
  end subroutine set_datainow


  !> destructor for when object goes out of scope
  subroutine destructor(self)
    type(neutraldata3D_fclaw) :: self

    ! deallocate arrays from base inputdata class
    call self%dissociate_pointers()

    ! null pointers specific to parent neutraldata class
    call self%dissociate_neutral_pointers()

    ! now deallocate arrays specific to this extension
    deallocate(self%proj_ezp_e1,self%proj_ezp_e2,self%proj_ezp_e3)
    deallocate(self%proj_eyp_e1,self%proj_eyp_e2,self%proj_eyp_e3)
    deallocate(self%proj_exp_e1,self%proj_exp_e2,self%proj_exp_e3)
    deallocate(self%ximat,self%yimat,self%zimat)

    ! root has some extra data
!    if (mpi_cfg%myid==0) then
!      deallocate(self%extents,self%indx,self%slabsizes)
!      deallocate(self%xnall,self%ynall)
!    end if

    ! set pointers to null
    nullify(self%xi,self%yi,self%zi);
    nullify(self%xn,self%yn,self%zn);
    nullify(self%dnO,self%dnN2,self%dnO2,self%dvnz,self%dvnx,self%dvny,self%dTn)
  end subroutine destructor
end module neutraldata3Dobj_fclaw
