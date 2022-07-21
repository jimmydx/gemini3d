!> non-mpi related grid quantities are stored here for use by other modules.  This could arguably be an
!    object...
module grid

use, intrinsic :: iso_c_binding, only: C_PTR,C_INT,c_loc,c_f_pointer
use meshobj, only: curvmesh
use meshobj_dipole, only: dipolemesh
use meshobj_cart, only: cartmesh
use phys_consts, only: wp
use reader, only: get_simsize3

integer, protected :: lx1,lx2,lx3,lx2all,lx3all
!! this is a useful shorthand for most program units using this module,
!! occasionally a program unit needs to define its own size in which case an only statement
!! is required when using this module.
integer, protected :: gridflag
!! for cataloguing the type of grid that we are using, open, closed, inverted, etc.  0 - closed dipole, 1 - inverted open, 2 - standard open.
real(wp), protected :: glonctr=-720._wp,glatctr=-180._wp

private
public :: lx1,lx2,lx3,lx2all,lx3all,gridflag, &
             get_grid3_coords_hdf5, &
             set_total_grid_sizes,set_subgrid_sizes,set_gridflag,grid_size, &
             grid_from_extents, grid_internaldata_alloc, grid_internaldata_generate, &
             grid_internaldata_ungenerate, &
             get_grid3_coords,read_size_gridcenter,detect_gridtype,set_size_gridcenter, &
             meshobj_alloc!, generate_worker_grid

interface ! readgrid_*.f90
  module subroutine get_grid3_coords_hdf5(path,x1,x2all,x3all,glonctr,glatctr)
    character(*), intent(in) :: path
    real(wp), dimension(:), intent(inout) :: x1,x2all,x3all
    real(wp), intent(out) :: glonctr,glatctr
  end subroutine
end interface

!! some overloading for situations needing vs. not needing an allocation step
!interface grid_from_extents
!  module procedure grid_from_extents_noalloc,grid_from_extents_alloc
!end interface

contains
  !> detect the type of grid that we are dealing with based solely on native coordinate values
  function detect_gridtype(x1,x2,x3) result(xtype)
    real(wp), dimension(-1:), intent(in) :: x1,x2,x3
    integer :: xtype

    if (maxval(abs(x2))<100) then
      print '(a)', 'Detected dipole grid...'
      xtype=2
    else
      print '(a)', 'Detected Cartesian grid...'
      xtype=1
    end if
  end function detect_gridtype


  !> Force a size and grid center location into module variables, if desired.  In general some other method
  !    should be used like read_size_gridcenter().  
  subroutine set_size_gridcenter(lx1in,lx2allin,lx3allin,glonctrin,glatctrin)
    integer, intent(in) :: lx1in,lx2allin,lx3allin
    real(wp), intent(in) :: glonctrin,glatctrin

    lx1=lx1in; lx2all=lx2allin; lx3all=lx3allin;
    glonctr=glonctrin; glatctr=glatctrin;
  end subroutine set_size_gridcenter


  !> Query the coordinates file and pull out the center geographic location for the entire grid (used for
  !    generation of Cartesian meshes and put in in a module-scope variable.
  subroutine read_size_gridcenter(indatsize,outdir)
    character(*), intent(in) :: indatsize,outdir
    real(wp), dimension(:), allocatable :: x1,x2all,x3all

    call get_simsize3(indatsize,lx1,lx2all,lx3all)
    allocate(x1(-1:lx2all+2),x2all(-1:lx2all+2),x3all(-1:lx3all+2))
    call get_grid3_coords(outdir,x1,x2all,x3all,glonctr,glatctr)
    deallocate(x1,x2all,x3all)
  end subroutine read_size_gridcenter


  !> Generate grid from a set of extents and sizes - e.g. similar to what is used in forestcalw.  input
  !    sizes should include ghost cells.  WARNING: this function will always just assume you are using a
  !    local grid, i.e. one that doesn't need knowledge of the full grid extents!  This requires that
  !    the grid type/class already be defined
  subroutine grid_from_extents(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,x)
    real(wp), dimension(2), intent(in) :: x1lims,x2lims,x3lims
    integer, intent(in) :: lx1wg,lx2wg,lx3wg
    class(curvmesh), intent(inout) :: x
    integer :: ix1,ix2,ix3
    real(wp), dimension(:), allocatable :: x1,x2,x3

    ! error checking
    if (glatctr<-90._wp .or. glatctr>90._wp) then
      error stop 'ERROR:grid_from_extents:  prior to calling must use read_size_gridcenter or set_size_gridcenter' // &
        'to assign module variables glonctr,glatctr'
    end if

    ! create temp space
    allocate(x1(lx1wg),x2(lx2wg),x3(lx3wg))

    ! make uniformly spaced coordinate arrays
    x1=[(x1lims(1) + (x1lims(2)-x1lims(1))/(lx1wg-1)*(ix1-1),ix1=1,lx1wg)]
    x2=[(x2lims(1) + (x2lims(2)-x2lims(1))/(lx2wg-1)*(ix2-1),ix2=1,lx2wg)]
    x3=[(x3lims(1) + (x3lims(2)-x3lims(1))/(lx3wg-1)*(ix3-1),ix3=1,lx3wg)]

    !call generate_worker_grid(x1,x2,x3,x2,x3,glonctr,glatctr,x)
    call grid_internaldata_alloc(x1,x2,x3,x2,x3,glonctr,glatctr,x)
    call grid_internaldata_generate(x)

    ! get rid of temp. arrays
    deallocate(x1,x2,x3)
  end subroutine grid_from_extents

!  ! FIXME: split into grid_alloc and generate_worker_grid; add interfaces for both to libgemini and C
!  !> this version additionally allocates the input argument, which is now a pointer
!  subroutine grid_from_extents_alloc(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,x,xtype,xC)
!    real(wp), dimension(2), intent(in) :: x1lims,x2lims,x3lims
!    integer, intent(in) :: lx1wg,lx2wg,lx3wg
!    class(curvmesh), intent(inout), pointer :: x
!    integer, intent(inout) :: xtype
!    type(c_ptr), intent(inout) :: xC
!    integer :: ix1,ix2,ix3
!    real(wp), dimension(:), allocatable :: x1,x2,x3
!
!    ! error checking
!    if (glatctr<-90._wp .or. glatctr>90._wp) then
!      error stop ' grid_from_extents:  prior to calling must use read_size_gridcenter or set_size_gridcenter to assign &
!                   module variables glonctr,glatctr'
!    end if
!
!    ! create temp space
!    allocate(x1(lx1wg),x2(lx2wg),x3(lx3wg))
!
!    ! make uniformly spaced coordinate arrays
!    x1=[(x1lims(1) + (x1lims(2)-x1lims(1))/(lx1wg-1)*(ix1-1),ix1=1,lx1wg)]
!    x2=[(x2lims(1) + (x2lims(2)-x2lims(1))/(lx2wg-1)*(ix2-1),ix2=1,lx2wg)]
!    x3=[(x3lims(1) + (x3lims(2)-x3lims(1))/(lx3wg-1)*(ix3-1),ix3=1,lx3wg)]
!
!    ! generate a subgrid from these
!    !if (present(xtype) .and. present(xC)) then   ! optional arguments confuses overloading for some types of calls :/
!      call meshobj_alloc(x1,x2,x3,x,xtype,xC)
!    !else
!    !  call meshobj_alloc(x1,x2,x3,x2,x3,x)
!    !end if
!    call generate_worker_grid(x1,x2,x3,x2,x3,glonctr,glatctr,x)
!
!    ! get rid of temp. arrays
!    deallocate(x1,x2,x3)
!  end subroutine grid_from_extents_alloc


  !> Find the type of the grid and allocate the correct type/class, return a C pointer if requested
  !    via optional arguments
  subroutine meshobj_alloc(x1,x2,x3,x,xtype,xC)
    real(wp), dimension(:), intent(in) :: x1,x2,x3
    class(curvmesh), pointer, intent(inout) :: x
    integer(C_INT), intent(inout), optional :: xtype
    type(C_PTR), intent(inout), optional :: xC
    integer :: gridtype
    type(cartmesh), pointer :: xcart
    type(dipolemesh), pointer :: xdipole

    !! allocate and read correct grid type
    gridtype=detect_gridtype(x1,x2,x3)
    select case (gridtype)
      case (2)
        !allocate(dipolemesh::x)
        allocate(xdipole)
        x=>xdipole
        if (present(xC) .and. present(xtype)) then
          xC = c_loc(xdipole)
          xtype = gridtype
        end if
      case (1)
        !allocate(cartmesh::x)
        allocate(xcart)
        x=>xcart
        if (present(xC) .and. present(xtype)) then
          xC = c_loc(xcart)
          xtype = gridtype
        end if
      case default
        error stop 'grid:meshobj_alloc - Unable to identify grid type'
    end select
  end subroutine 


!  !> Generate a "worker" grid based on coordinate arrays and grid center, polymorphic grid object must already
!  !    exist, i.e. already be allocated with some dynamic type.  Note that you can set x2all=x2 and
!  !    (or) x3all=x3 if you are only doing "local" grid operations in your GEMINI application, e.g. as with
!  !    trees-GEMINI.  The dynamic type of x must be set prior to calling this function; this can be 
!  !    accomplished e.g. through a wrapper
!  subroutine generate_worker_grid(x1,x2,x3,x2all,x3all,glonctr,glatctr,x)
!    real(wp), dimension(:), intent(in) :: x1,x2,x3,x2all,x3all
!    real(wp), intent(in) :: glonctr,glatctr
!    class(curvmesh), intent(inout) :: x
!
!    ! Create the grid object
!    call x%set_center(glonctr,glatctr)
!    call x%set_coords(x1,x2,x3,x2all,x3all)    ! store coordinate arrays
!    call x%init()                              ! allocate space for subgrid variables
!    call x%make()                              ! fill auxiliary arrays
!
!    call set_gridflag(x%gridflag)
!  end subroutine generate_worker_grid


  !> Trigger allocation of grid class internal data once the class itself has been allocated and typed
  subroutine grid_internaldata_alloc(x1,x2,x3,x2all,x3all,glonctr,glatctr,x)
    real(wp), dimension(:), intent(in) :: x1,x2,x3,x2all,x3all
    real(wp), intent(in) :: glonctr,glatctr
    class(curvmesh), intent(inout) :: x

    call x%set_center(glonctr,glatctr)         ! set center location for grid (in case used, e.g. for Cartesian)
    call x%set_coords(x1,x2,x3,x2all,x3all)    ! store coordinate arrays
    call x%init()                              ! allocate space for subgrid variables
  end subroutine grid_internaldata_alloc


  !> Trigger a generation of all grid internal data
  subroutine grid_internaldata_generate(x)
    class(curvmesh), intent(inout) :: x

    call x%make()                              ! trigger generation of all internal data arrays
    call set_gridflag(x%gridflag)              ! set module variable to match the type stored in the grid class
  end subroutine 


  !> Force deallocation of grid data at least to the point where it can be "remade", e.g. for AMR-like operations
  subroutine grid_internaldata_ungenerate(x)
    class(curvmesh), intent(inout) :: x

    call x%dissociate_pointers()
  end subroutine grid_internaldata_ungenerate


  !> Read in native coordinates from a grid file
  subroutine get_grid3_coords(path,x1,x2all,x3all,glonctr,glatctr)
    character(*), intent(in) :: path
    real(wp), dimension(:), intent(inout) :: x1,x2all,x3all
    real(wp) :: glonctr,glatctr

    character(:), allocatable :: fmt

    fmt = path(index(path, '.', back=.true.) : len(path))
    select case (fmt)
      case ('.h5')
        call get_grid3_coords_hdf5(path,x1,x2all,x3all,glonctr,glatctr)
      case default
        error stop 'grid:read:get_grid3: unknown grid format: ' // fmt
    end select

    if(size(x1) < 1) error stop 'grid:get_grid3_coords: size(x1) must be strictly positive'
    if(size(x2all) < 1) error stop 'grid:get_grid3_coords: size(x2all) must be strictly positive'
    if(size(x3all) < 1) error stop 'grid:get_grid3_coords: size(x3all) must be strictly positive'
  end subroutine get_grid3_coords


  subroutine set_total_grid_sizes(lx1in,lx2allin,lx3allin)
    integer, intent(in) :: lx1in,lx2allin,lx3allin

    lx1=lx1in; lx2all=lx2allin; lx3all=lx3allin;
  end subroutine set_total_grid_sizes


  subroutine set_subgrid_sizes(lx2in,lx3in)
    integer, intent(in) :: lx2in,lx3in

    lx2=lx2in; lx3=lx3in;
  end subroutine set_subgrid_sizes


  subroutine set_gridflag(gridflagin)
    integer, intent(in) :: gridflagin

    gridflag=gridflagin
  end subroutine set_gridflag


!  subroutine bind_grav_ptrs(g1in,g2in,g3in)
!    real(wp), dimension(:,:,:), pointer, intent(in) :: g1in,g2in,g3in
!
!    g1=>g1in; g2=>g2in; g3=>g3in
!  end subroutine bind_grav_ptrs


  subroutine grid_size(indatsize)
  !! CHECK THE SIZE OF THE GRID TO BE LOADED AND SET SIZES IN THIS MODULE (NOT IN STRUCTURE THOUGH)
    character(*), intent(in) :: indatsize

    call get_simsize3(indatsize, lx1, lx2all, lx3all)
    print *, 'grid_size_root: full grid size:  ',lx1,lx2all,lx3all
    call set_total_grid_sizes(lx1,lx2all,lx3all)    !! set module global sizes for use on other contexts
  end subroutine grid_size
end module grid
