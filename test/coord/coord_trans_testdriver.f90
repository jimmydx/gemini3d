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
   real(wp) :: lat0, lon0, h0                 ! Geodetic test point [radians, m]
   real(wp) :: lat0_rec, lon0_rec, h0_rec     ! Geodetic test point [radians, m] (recovered)
   real(wp) :: x0, y0, z0                     ! ECEF coordinates
   real(wp) :: x0_rec, y0_rec, z0_rec         ! ECEF coordinates (recovered)
   real(wp) :: east0, north0, up0             ! ENU coordinates
   real(wp) :: east0_rec, north0_rec, up0_rec ! ENU coordinates (recovered)
   real(wp) :: az0, el0, r0_arr               ! AER coordinates
   real(wp), dimension(3) :: err_geo_sp, err_ecef_sp, err_enu_sp, err_aer_sp
 
   real(wp) :: alt0, r0, phi0, theta0           ! Magnetic test angles (radians)
   real(wp) :: glon0, glat0, phi0_rec, theta0_rec
   real(wp), dimension(3) :: err_mag_sp 
   real(wp) :: r_alt0, alt0_rec, r0_rec ! altitude test values
 
   real(wp) :: R_err
   real(wp), dimension(3,3) :: R, I

   integer, parameter :: nlat = 10, nlon = 10, nh = 10
   real(wp), dimension(2) :: latlims, lonlims, hlims
   integer :: ii, jj, kk, N, idx ! Loop indices and counters
   real(wp), dimension(:,:,:), allocatable :: lat_grid, lon_grid, h_grid
   real(wp), dimension(:,:,:), allocatable :: x_grid, y_grid, z_grid
   real(wp), dimension(:,:,:), allocatable :: lat_grid_rec, lon_grid_rec, h_grid_rec
   real(wp), dimension(:), allocatable :: latvec, lonvec, hvec

   integer, parameter :: nrad = 10, ntheta = 10, nphi = 10
   real(wp), dimension(2) :: rlims, thetalims, philims
   real(wp), dimension(:,:,:), allocatable :: r_arr, alt_arr, theta_arr, phi_arr
   real(wp), dimension(:,:,:), allocatable :: glon_arr, glat_arr, r_rec_arr, phi_rec_arr, theta_rec_arr
   real(wp), dimension(:), allocatable :: rvec, theta_vec, phi_vec
   
   ! Define test values:
   lat0 = 30.0_wp * deg2rad
   lon0 = -80.0_wp * deg2rad
   h0   = 80.0e3_wp
   r0 = Re + 80.0e3_wp
   phi0 = 5.0_wp * deg2rad
   theta0 = 5.0_wp * deg2rad

   latlims=[lat0,lat0+10.0*deg2rad]
   lonlims=[lon0,lon0+10.0*deg2rad]
   hlims=[h0,h0+10.0e3]
    
   rlims = [r0, r0 + 900.0e3_wp]
   thetalims = [theta0, pi - theta0]
   philims = [phi0, 2.0_wp * pi - phi0]
 
   print *, "SINGLE-POINT TESTS"
   print *, "Geodetic Test Point: ", lat0/deg2rad, lon0/deg2rad, h0
   print *, "Magnetic Test Values: phi = ", phi0, " theta = ", theta0
   print *

   call set_earth_model(SPHERICAL)
 
   ! Test 1: Geodetic <-> ECEF conversion
  call geodetic_to_ecef(lat0, lon0, h0, x0, y0, z0)
  call ecef_to_geodetic(x0, y0, z0, lat0_rec, lon0_rec, h0_rec)
  err_geo_sp = [ abs(lat0 - lat0_rec), abs(lon0 - lon0_rec), abs(h0 - h0_rec) ]
  if (any(err_geo_sp > tol)) then
     print *, "Error in Single-Point Geodetic <-> ECEF conversion: ", err_geo_sp
     stop 1
  else
     write(*,'(A,ES12.5)') "Single-Point Geodetic <-> ECEF conversion passed with max error: ", maxval(err_geo_sp)
  end if

  ! Test 2: ECEF <-> ENU conversion (reference = test point)
  call ecef_to_enu(x0, y0, z0, lat0, lon0, h0, east0, north0, up0)
  call enu_to_ecef(east0, north0, up0, lat0, lon0, h0, x0_rec, y0_rec, z0_rec)
  err_ecef_sp = [ abs(x0 - x0_rec), abs(y0 - y0_rec), abs(z0 - z0_rec) ]
  if (any(err_ecef_sp > tol)) then
     print *, "Error in Single-Point ECEF <-> ENU conversion: ", err_ecef_sp
     stop 1
  else
     write(*,'(A,ES12.5)') "Single-Point ECEF <-> ENU conversion passed with max error: ", maxval(err_ecef_sp)
  end if

  ! Test 3: ENU <-> AER conversion
  call enu_to_aer(east0, north0, up0, az0, el0, r0_arr)
  call aer_to_enu(az0, el0, r0_arr, east0_rec, north0_rec, up0_rec)
  err_enu_sp = [ abs(east0 - east0_rec), abs(north0 - north0_rec), abs(up0 - up0_rec) ]
  if (any(err_enu_sp > tol)) then
     print *, "Error in Single-Point ENU <-> AER conversion: ", err_enu_sp
     stop 1
  else
     write(*,'(A,ES12.5)') "Single-Point ENU <-> AER conversion passed with max error: ", maxval(err_enu_sp)
  end if

  ! Test 4: ECEF <-> AER conversion
  call ecef_to_aer(x0, y0, z0, lat0, lon0, h0, az0, el0, r0_arr)
  call aer_to_ecef(az0, el0, r0_arr, lat0, lon0, h0, x0_rec, y0_rec, z0_rec)
  err_aer_sp = [ abs(x0 - x0_rec), abs(y0 - y0_rec), abs(z0 - z0_rec) ]
  if (any(err_aer_sp > tol)) then
     print *, "Error in Single-Point ECEF <-> AER conversion: ", err_aer_sp
     stop 1
  else
     write(*,'(A,ES12.5)') "Single-Point ECEF <-> AER conversion passed with max error: ", maxval(err_aer_sp)
  end if

  ! Test 5: Geomagnetic Conversion Test
  ! Convert from geomagnetic to geographic and then back.
  ! Here we assume the starting geomagnetic coordinates are given by phi0 and theta0.
  call geomag2geog(phi0, theta0, glon0, glat0)
  alt0 = r2alt(r0)
  call geog2geomag(glon0, glat0, phi0_rec, theta0_rec)
  r0_rec = alt2r(alt0)
  err_mag_sp = [ abs(r0 - r0_rec), abs(phi0 - phi0_rec), abs(theta0 - theta0_rec) ]
  if ( any(err_mag_sp > tol) ) then
     print *, "Error in Single-Point geomagnetic conversion: ", err_mag_sp
     stop 1
  else
     write(*,'(A,ES12.5)') "Single-Point geomagnetic <-> geographic conversion passed with max error: ", maxval(err_mag_sp)
  end if

  ! Test 6: Altitude Conversion Test
  alt0 = h0
  r_alt0 = alt2r(alt0)
  alt0_rec = r2alt(r_alt0)
  if (abs(alt0 - alt0_rec) > tol) then
     print *, "Error in Single-Point altitude conversion:", abs(alt0 - alt0_rec)
     stop 1
  else
     write(*,'(A,ES12.5)') "Single-Point altitude conversion passed with error: ", abs(alt0 - alt0_rec)
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
     write(*,'(A,ES12.5)') "Single-Point rotation matrix test passed with error: ", R_err
  end if


  !**********************************************************************
   ! GEODETIC GRID TEST (3D grid with shape (nlat, nlon, nh))
   !**********************************************************************

allocate(lat_grid(nlat,nlon,nh), lon_grid(nlat,nlon,nh), h_grid(nlat,nlon,nh))
allocate(x_grid(nlat,nlon,nh), y_grid(nlat,nlon,nh), z_grid(nlat,nlon,nh))
allocate(lat_grid_rec(nlat,nlon,nh), lon_grid_rec(nlat,nlon,nh), h_grid_rec(nlat,nlon,nh))
allocate(latvec(nlat), lonvec(nlon), hvec(nh))

latvec = [(lat0 + (latlims(2)-latlims(1))/nlat*(ii-1),ii=1,nlat)]
lonvec = [(lon0 + (lonlims(2)-lonlims(1))/nlon*(jj-1),jj=1,nlon)]
hvec = [(h0 + (hlims(2)-hlims(1))/nh*(kk-1),kk=1,nh)]

do ii = 1, nlat
   do jj = 1, nlon
      do kk = 1, nh
          lat_grid(ii,jj,kk) = latvec(ii)
          lon_grid(ii,jj,kk) = lonvec(jj)
          h_grid(ii,jj,kk) = hvec(kk)
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
     write(*,'(A,3(ES12.5))') "Geodetic Grid conversion passed with max errors: ", &
              maxval(abs(lat_grid - lat_grid_rec)), maxval(abs(lon_grid - lon_grid_rec)), &
              maxval(abs(h_grid - h_grid_rec))
  end if
  deallocate(latvec, lonvec, hvec, lat_grid, lon_grid, h_grid, x_grid, y_grid, z_grid, &
             lat_grid_rec, lon_grid_rec, h_grid_rec)

  !**********************************************************************
  ! GEOMAGNETIC GRID TEST (3D grid with shape (nrad, ntheta, nphi))
  !**********************************************************************
  allocate(r_arr(nrad, ntheta, nphi))
  allocate(theta_arr, phi_arr, glon_arr, glat_arr, phi_rec_arr, theta_rec_arr, mold=r_arr)
  allocate(rvec(nrad))
  allocate(theta_vec(ntheta))
  allocate(phi_vec(nphi))

  rvec = [ (r0 + (rlims(2) - rlims(1))/nrad*(ii-1), ii=1,nrad) ]
  theta_vec = [ (theta0 + (thetalims(2) - thetalims(1))/ntheta*(jj-1), jj=1,ntheta) ]
  phi_vec = [ (phi0 + (philims(2) - philims(1))/nphi*(kk-1), kk=1,nphi) ]

  do ii = 1, nrad
   do jj = 1, ntheta
      do kk = 1, nphi
          r_arr(ii,jj,kk) = rvec(ii)
          theta_arr(ii,jj,kk) = theta_vec(jj)
          phi_arr(ii,jj,kk) = phi_vec(kk)
      end do
   end do
end do
  
  call geomag2geog(phi_arr, theta_arr, glon_arr, glat_arr)
  alt_arr = r2alt(r_arr)
  call geog2geomag(glon_arr, glat_arr, phi_rec_arr, theta_rec_arr)
  r_rec_arr = alt2r(alt_arr)
  if ( any([maxval(abs(r_arr - r_rec_arr)),maxval(abs(phi_arr - phi_rec_arr)),maxval(abs(theta_arr - theta_rec_arr))] > tol) ) then
     print *, "Error in Geomagnetic Grid conversion"
     stop 1
  else
     write(*,'(A,3(ES12.5))') "Geomagnetic Grid conversion passed with max errors: ", &
             maxval(abs(r_arr - r_rec_arr)), maxval(abs(phi_arr - phi_rec_arr)), maxval(abs(theta_arr - theta_rec_arr))
  end if
  deallocate(theta_vec, phi_vec, rvec, r_arr, theta_arr, phi_arr, glon_arr, glat_arr, &
             phi_rec_arr, theta_rec_arr)

  !**********************************************************************
  ! Deallocate all Single-Point Arrays and Matrices
  !**********************************************************************
  ! Note: All allocatable variables in this program are deallocated automatically at program end.
  
end program coord_trans_testdriver