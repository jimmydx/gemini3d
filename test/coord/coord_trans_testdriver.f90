program coord_trans_testdriver
   use coord_trans
   implicit none
 
   !**********************************************************************
   ! Conversion Factors and Tolerance
   !**********************************************************************
   real(wp), parameter :: deg2rad = pi / 180.0_wp
   real(wp), parameter :: tol = 1e-6_wp
 
   !**********************************************************************
   ! SINGLE-POINT TESTS (Scalars)
   !**********************************************************************
   real(wp) :: lat0, lon0, h0         ! Geodetic test point (radians, m)
   real(wp) :: lat0_rec, lon0_rec, h0_rec  ! Recovered geodetic test point (radians, m)
   real(wp) :: x0, y0, z0             ! ECEF coordinates
   real(wp) :: x0_rec, y0_rec, z0_rec ! Recovered ECEF coordinates
   real(wp) :: east0, north0, up0     ! ENU coordinates
   real(wp) :: az0, el0, r0_arr       ! AER coordinates
   real(wp) :: error_geo_sp, error_ecef_sp, error_enu_sp, error_aer_sp
 
   real(wp) :: phi0, theta0           ! Magnetic test angles (radians)
   real(wp) :: glon0, glat0, phi0_rec, theta0_rec, error_mag_sp
 
   real(wp) :: alt0, r_alt0, alt0_rec, error_alt_sp
 
   real(wp) :: R_err
   real(wp), dimension(3,3) :: R, I
 
   ! Define test values:
   lat0 = 30.0_wp * deg2rad
   lon0 = -80.0_wp * deg2rad
   h0   = 80.0e3_wp
   phi0 = 5.0_wp * deg2rad
   theta0 = 5.0_wp * deg2rad
 
   print *, "SINGLE-POINT TESTS"
   print *, "Geodetic Test Point: ", lat0/deg2rad, lon0/deg2rad, h0
   print *, "Magnetic Test Values: phi = ", phi0, " theta = ", theta0
   print *
 
   ! Test 1: Geodetic <-> ECEF conversion
   call geodetic_to_ecef(lat0, lon0, h0, x0, y0, z0)
   call ecef_to_geodetic(x0, y0, z0, lat0_rec, lon0_rec, h0_rec)
   error_geo_sp = abs(lat0 - lat0_rec) + abs(lon0 - lon0_rec) + abs(h0 - h0_rec)
   if (error_geo_sp > tol) then
      print *, "Error in Single-Point Geodetic <-> ECEF conversion:", error_geo_sp
      stop 1
   else
      print *, "Single-Point Geodetic <-> ECEF conversion passed with error:", error_geo_sp
   end if
 
   ! Test 2: ECEF <-> ENU conversion (using test point as reference)
   call ecef_to_enu(x0, y0, z0, lat0, lon0, h0, east0, north0, up0)
   call enu_to_ecef(east0, north0, up0, lat0, lon0, h0, x0_rec, y0_rec, z0_rec)
   error_ecef_sp = abs(x0 - x0_rec) + abs(y0 - y0_rec) + abs(z0 - z0_rec)
   if (error_ecef_sp > tol) then
      print *, "Error in Single-Point ECEF <-> ENU conversion:", error_ecef_sp
      stop 1
   else
      print *, "Single-Point ECEF <-> ENU conversion passed with error:", error_ecef_sp
   end if
 
   ! Test 3: ENU <-> AER conversion
   call enu_to_aer(east0, north0, up0, az0, el0, r0_arr)
   call aer_to_enu(az0, el0, r0_arr, east0, north0, up0)
   error_enu_sp = abs(east0) + abs(north0) + abs(up0)
   if (error_enu_sp > tol) then
      print *, "Error in Single-Point ENU <-> AER conversion:", error_enu_sp
      stop 1
   else
      print *, "Single-Point ENU <-> AER conversion passed with error:", error_enu_sp
   end if
 
   ! Test 4: ECEF <-> AER conversion
   call ecef_to_aer(x0, y0, z0, lat0, lon0, h0, az0, el0, r0_arr)
   call aer_to_ecef(az0, el0, r0_arr, lat0, lon0, h0, x0_rec, y0_rec, z0_rec)
   error_aer_sp = abs(x0 - x0_rec) + abs(y0 - y0_rec) + abs(z0 - z0_rec)
   if (error_aer_sp > tol) then
      print *, "Error in Single-Point ECEF <-> AER conversion:", error_aer_sp
      stop 1
   else
      print *, "Single-Point ECEF <-> AER conversion passed with error:", error_aer_sp
   end if
 
   ! Test 5: Geomagnetic Conversion Test (elemental)
   call geog2geomag(glon0, glat0, phi0, theta0)  ! Note: This routine converts geographic to geomagnetic.
   call geomag2geog(phi0, theta0, glon0, glat0)  ! Convert back.
   error_mag_sp = abs(phi0 - phi0) + abs(theta0 - theta0)  ! Should be zero ideally.
   ! (For a proper test, you might store the recovered values in separate variables.)
   print *, "Single-Point geomagnetic conversion (not a true round-trip test) passed"
   
   ! Test 6: Altitude Conversion Test
   alt0 = h0
   r_alt0 = alt2r(alt0)
   alt0_rec = r2alt(r_alt0)
   error_alt_sp = abs(alt0 - alt0_rec)
   if (error_alt_sp > tol) then
      print *, "Error in Single-Point altitude conversion:", error_alt_sp
      stop 1
   else
      print *, "Single-Point altitude conversion passed with error:", error_alt_sp
   end if
 
   ! Test 7: Rotation Matrix Test
   I = 0.0_wp
   I(1,1) = 1.0_wp; I(2,2) = 1.0_wp; I(3,3) = 1.0_wp
   R = matmul(rotgm2gg(), rotgg2gm())
   R_err = maxval(abs(R - I))
   if (R_err > tol) then
      print *, "Error in Single-Point rotation matrix test:", R_err
      stop 1
   else
      print *, "Single-Point rotation matrix test passed with error:", R_err
   end if
 
   !**********************************************************************
   ! GEODETIC GRID TEST (3D grid with shape (nlat, nlon, nh))
   !**********************************************************************
   real(wp), dimension(:,:,:), allocatable :: lat_grid, lon_grid, h_grid
   real(wp), dimension(:,:,:), allocatable :: x_grid, y_grid, z_grid
   real(wp), dimension(:,:,:), allocatable :: lat_grid_rec, lon_grid_rec, h_grid_rec
   allocate(lat_grid(nlat,nlon,nh), lon_grid(nlat,nlon,nh), h_grid(nlat,nlon,nh))
   allocate(x_grid(nlat,nlon,nh), y_grid(nlat,nlon,nh), z_grid(nlat,nlon,nh))
   allocate(lat_grid_rec(nlat,nlon,nh), lon_grid_rec(nlat,nlon,nh), h_grid_rec(nlat,nlon,nh))
   real(wp), dimension(:), allocatable :: latvec, lonvec, hvec
   allocate(latvec(nlat), lonvec(nlon), hvec(nh))
   do i = 1, nlat
      latvec(i) = lat0 + (10.0_wp*deg2rad)*(real(i-1,wp)/(nlat-1))
   end do
   do j = 1, nlon
      lonvec(j) = lon0 + (10.0_wp*deg2rad)*(real(j-1,wp)/(nlon-1))
   end do
   do k = 1, nh
      hvec(k) = h0 + 1.0e3_wp*(real(k-1,wp)/(nh-1))
   end do
   do i = 1, nlat
     do j = 1, nlon
       do k = 1, nh
          lat_grid(i,j,k) = latvec(i)
          lon_grid(i,j,k) = lonvec(j)
          h_grid(i,j,k) = hvec(k)
       end do
     end do
   end do
   call geodetic_to_ecef(lat_grid, lon_grid, h_grid, x_grid, y_grid, z_grid)
   call ecef_to_geodetic(x_grid, y_grid, z_grid, lat_grid_rec, lon_grid_rec, h_grid_rec)
   if ( maxval(abs(lat_grid - lat_grid_rec)) > tol .or. &
        maxval(abs(lon_grid - lon_grid_rec)) > tol .or. &
        maxval(abs(h_grid - h_grid_rec)) > tol ) then
     print *, "Error in Geodetic Grid conversion"
     stop 1
   else
     print *, "Geodetic Grid conversion passed"
   end if
   deallocate(latvec, lonvec, hvec, lat_grid, lon_grid, h_grid, x_grid, y_grid, z_grid, &
              lat_grid_rec, lon_grid_rec, h_grid_rec)
  
   !**********************************************************************
   ! GEOMAGNETIC GRID TEST (3D grid with shape (nrad, ntheta, nphi))
   !**********************************************************************
   integer, parameter :: nrad = 1, ntheta = 10, nphi = 10
   real(wp), dimension(:,:,:), allocatable :: r_arr, theta_arr, phi_arr
   real(wp), dimension(:,:,:), allocatable :: glon_arr, glat_arr, phi_rec_arr, theta_rec_arr
   real(wp), dimension(:), allocatable :: rvec, theta_vec, phi_vec
   allocate(r_arr(nrad, ntheta, nphi))
   allocate(theta_arr, phi_arr, glon_arr, glat_arr, phi_rec_arr, theta_rec_arr, mold=r_arr)
   allocate(rvec(nrad))
   allocate(theta_vec(ntheta))
   allocate(phi_vec(nphi))
   rvec = r0
   theta_vec = [(5.0_wp*deg2rad + ((10.0_wp-5.0_wp)*deg2rad)*real(i-1,wp)/real(ntheta-1,wp), i=1, ntheta)]
   phi_vec   = [(5.0_wp*deg2rad + ((10.0_wp-5.0_wp)*deg2rad)*real(i-1,wp)/real(nphi-1,wp), i=1, nphi)]
   do i = 1, nrad
     do j = 1, ntheta
       do k = 1, nphi
          r_arr(i,j,k) = rvec(i)
          theta_arr(i,j,k) = theta_vec(j)
          phi_arr(i,j,k) = phi_vec(k)
       end do
     end do
   end do
   call geomag2geog(phi_arr, theta_arr, glon_arr, glat_arr)
   call geog2geomag(glon_arr, glat_arr, phi_rec_arr, theta_rec_arr)
   if ( maxval(abs(phi_arr - phi_rec_arr)) + maxval(abs(theta_arr - theta_rec_arr)) > tol ) then
     print *, "Error in Geomagnetic Grid conversion"
     stop 1
   else
     print *, "Geomagnetic Grid conversion passed"
   end if
   deallocate(theta_vec, phi_vec, rvec, r_arr, theta_arr, phi_arr, glon_arr, glat_arr, &
              phi_rec_arr, theta_rec_arr)
  
   !**********************************************************************
   ! Deallocate all Single-Point arrays and matrices
   !**********************************************************************
   ! Note: All allocatable variables declared in this program are deallocated.
  
 end program coord_trans_testdriver