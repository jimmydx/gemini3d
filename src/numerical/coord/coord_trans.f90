!======================================================================
! Module: coord_trans
! Purpose: Provides coordinate transformation routines for converting
!          geodetic, ECEF, ENU, AER, and geomagnetic coordinates.
!======================================================================
module coord_trans  ! Module for coordinate transformations (geodetic, ECEF, ENU, AER, geomagnetic)
  implicit none  ! Ensure all variables are explicitly declared

  !----------------------------------------------------------------------
  ! Public interface
  !----------------------------------------------------------------------
  public :: set_earth_model, geodetic_to_ecef, ecef_to_geodetic, ecef_to_enu, enu_to_ecef, &
            enu_to_aer, aer_to_enu, ecef_to_aer, aer_to_ecef, & ! ECEFspher2ENU, &
            geog2geomag, geomag2geog, r2alt, alt2r, rotz, roty, rotgm2gg, rotgg2gm

  integer, parameter :: wp = selected_real_kind(15, 307)

  !======================================================================
  ! Earth Model and Ellipsoid Parameters
  !======================================================================
  !---------------------------
  ! WGS84 Ellipsoid Parameters
  !---------------------------
  real(wp), parameter :: a  = 6378137.0_wp  ! Semi-major axis (WGS84)
  real(wp), parameter :: f  = 1.0_wp/298.257223563_wp  ! Flattening factor (WGS84)
  real(wp), parameter :: b  = a*(1.0_wp - f)  ! Semi-minor axis computed from a and f
  real(wp), parameter :: e2 = (a**2 - b**2)/a**2  ! Square of eccentricity
  real(wp), parameter :: ep2= (a**2 - b**2)/b**2  ! Second eccentricity squared

  !---------------------------
  ! Spherical Earth Parameter
  !---------------------------
  real(wp), parameter :: Re = 6371000.0_wp  ! Mean Earth radius for spherical model

  !---------------------------
  ! Earth Model Types
  !---------------------------
  integer, parameter :: SPHERICAL = 0, WGS84 = 1  ! Earth model types: 0 for SPHERICAL, 1 for WGS84
  integer :: earth_model = WGS84  ! Default earth model is WGS84

  !---------------------------
  ! Fundamental Constants
  !---------------------------
  real(wp), parameter :: pi = 3.14159265358979323846_wp  ! Pi constant
  real(wp), parameter :: two_pi = 2.0_wp*pi  ! 2*pi constant

  !---------------------------
  ! Geomagnetic Parameters
  !---------------------------
  real(wp), parameter :: thetan = 11.0_wp*pi/180.0_wp  ! Geomagnetic colatitude in radians
  real(wp), parameter :: phin   = 289.0_wp*pi/180.0_wp  ! Geomagnetic longitude in radians

  !======================================================================
  ! Module Procedures: Subroutines and Functions
  !======================================================================
contains

  !----------------------------------------------------------------------
  ! Subroutine: set_earth_model
  ! Purpose: Sets the earth model type based on the input value.
  !----------------------------------------------------------------------
  subroutine set_earth_model(model_choice)
    integer, intent(in) :: model_choice
    if (model_choice == SPHERICAL .or. model_choice == WGS84) then  ! Check if the provided earth model is valid
      earth_model = model_choice  ! Set earth model to the chosen valid option
    else
      print*, "Invalid earth model choice. Defaulting to SPHERICAL."  ! Warn invalid choice
      earth_model = SPHERICAL  ! Default to SPHERICAL model
    end if
  end subroutine set_earth_model

  !----------------------------------------------------------------------
  ! Subroutine: geodetic_to_ecef
  ! Purpose: Converts geodetic coordinates (lat, lon, h) to ECEF coordinates (x, y, z)
  !----------------------------------------------------------------------
  elemental subroutine geodetic_to_ecef(lat, lon, h, x, y, z)
    real(wp), intent(in) :: lat, lon, h
    real(wp), intent(out) :: x, y, z
    if (earth_model == WGS84) then
      x = ( a / sqrt(1.0_wp - e2*sin(lat)**2) + h ) * cos(lat) * cos(lon)  ! Compute x coordinate using WGS84 parameters
      y = ( a / sqrt(1.0_wp - e2*sin(lat)**2) + h ) * cos(lat) * sin(lon)  ! Compute y coordinate using WGS84 parameters
      z = ( a / sqrt(1.0_wp - e2*sin(lat)**2)*(1.0_wp - e2) + h ) * sin(lat)  ! Compute z coordinate using WGS84 parameters
    else
      x = ( Re + h ) * cos(lat) * cos(lon)  ! Compute x coordinate using spherical Earth model
      y = ( Re + h ) * cos(lat) * sin(lon)  ! Compute y coordinate using spherical Earth model
      z = ( Re + h ) * sin(lat)  ! Compute z coordinate using spherical Earth model
    end if
  end subroutine geodetic_to_ecef

  !----------------------------------------------------------------------
  ! Subroutine: ecef_to_geodetic
  ! Purpose: Converts ECEF coordinates (x, y, z) to geodetic coordinates (lat, lon, h)
  !----------------------------------------------------------------------
  elemental subroutine ecef_to_geodetic(x, y, z, lat, lon, h)
    real(wp), intent(in) :: x, y, z
    real(wp), intent(out) :: lat, lon, h
    real(wp) :: p, theta, N_local
    if (earth_model == WGS84) then
      p = sqrt(x**2 + y**2)  ! Compute horizontal distance from Earth's axis
      lon = atan2(y, x)  ! Determine longitude using atan2
      if (p == 0.0_wp) then
        if (z >= 0.0_wp) then
          lat = pi/2.0_wp  ! North pole case
        else
          lat = -pi/2.0_wp  ! South pole case
        end if
        h = abs(z) - b  ! Altitude computed from z and semi-minor axis
      else
        theta = atan2(z*a, p*b)  ! Estimate parametric latitude (theta)
        lat = atan2(z + ep2*b*sin(theta)**3, p - e2*a*cos(theta)**3)  ! Compute geodetic latitude
        N_local = a / sqrt(1.0_wp - e2*sin(lat)**2)  ! Local radius of curvature
        h = p/cos(lat) - N_local  ! Calculate altitude above ellipsoid
      end if
    else
      p = sqrt(x**2 + y**2)  ! Compute horizontal distance
      lon = atan2(y, x)  ! Determine longitude
      if (p == 0.0_wp) then
        if (z >= 0.0_wp) then
          lat = pi/2.0_wp  ! Handle pole case
        else
          lat = -pi/2.0_wp
        end if
        h = abs(z) - Re  ! Altitude from spherical model
      else
        lat = atan2(z, p)  ! Compute latitude for spherical Earth
        h = sqrt(x**2 + y**2 + z**2) - Re  ! Altitude above spherical Earth
      end if
    end if
  end subroutine ecef_to_geodetic

  !----------------------------------------------------------------------
  ! Subroutine: ecef_to_enu
  ! Purpose: Converts ECEF coordinates (x, y, z) to local ENU coordinates (east, north, up)
  !          relative to a reference point (ref_lat, ref_lon, ref_h)
  !----------------------------------------------------------------------
  elemental subroutine ecef_to_enu(x, y, z, ref_lat, ref_lon, ref_h, east, north, up)
    real(wp), intent(in) :: x, y, z
    real(wp), intent(in) :: ref_lat, ref_lon, ref_h
    real(wp), intent(out) :: east, north, up
    real(wp) :: dx, dy, dz
    dx = x - ((a / sqrt(1.0_wp - e2*sin(ref_lat)**2) + ref_h)*cos(ref_lat)*cos(ref_lon))  ! x difference from reference point
    dy = y - ((a / sqrt(1.0_wp - e2*sin(ref_lat)**2) + ref_h)*cos(ref_lat)*sin(ref_lon))  ! y difference from reference point
    dz = z - (((a*(1.0_wp - e2)) / sqrt(1.0_wp - e2*sin(ref_lat)**2) + ref_h)*sin(ref_lat))  ! z difference from reference point
    east = -sin(ref_lon)*dx + cos(ref_lon)*dy  ! Compute east component
    north = -sin(ref_lat)*cos(ref_lon)*dx - sin(ref_lat)*sin(ref_lon)*dy + cos(ref_lat)*dz  ! Compute north component
    up = cos(ref_lat)*cos(ref_lon)*dx + cos(ref_lat)*sin(ref_lon)*dy + sin(ref_lat)*dz  ! Compute up component
  end subroutine ecef_to_enu

  !----------------------------------------------------------------------
  ! Subroutine: enu_to_ecef
  ! Purpose: Converts local ENU coordinates (east, north, up) back to ECEF coordinates (x, y, z)
  !          relative to a reference point (ref_lat, ref_lon, ref_h)
  !----------------------------------------------------------------------
  elemental subroutine enu_to_ecef(east, north, up, ref_lat, ref_lon, ref_h, x, y, z)
    real(wp), intent(in) :: east, north, up
    real(wp), intent(in) :: ref_lat, ref_lon, ref_h
    real(wp), intent(out) :: x, y, z
    real(wp) :: dx, dy, dz
    dx = -sin(ref_lon)*east - sin(ref_lat)*cos(ref_lon)*north + cos(ref_lat)*cos(ref_lon)*up  ! Inverse rotation: x offset
    dy = cos(ref_lon)*east - sin(ref_lat)*sin(ref_lon)*north + cos(ref_lat)*sin(ref_lon)*up   ! Inverse rotation: y offset
    dz = cos(ref_lat)*north + sin(ref_lat)*up  ! Inverse rotation: z offset
    x = ((a / sqrt(1.0_wp - e2*sin(ref_lat)**2) + ref_h)*cos(ref_lat)*cos(ref_lon)) + dx  ! Add reference x-coordinate
    y = ((a / sqrt(1.0_wp - e2*sin(ref_lat)**2) + ref_h)*cos(ref_lat)*sin(ref_lon)) + dy  ! Add reference y-coordinate
    z = (((a*(1.0_wp - e2)) / sqrt(1.0_wp - e2*sin(ref_lat)**2) + ref_h)*sin(ref_lat)) + dz  ! Add reference z-coordinate
  end subroutine enu_to_ecef

  !----------------------------------------------------------------------
  ! Subroutine: enu_to_aer
  ! Purpose: Converts ENU coordinates (east, north, up) to AER coordinates (azimuth, elevation, range)
  !----------------------------------------------------------------------
  elemental subroutine enu_to_aer(east, north, up, az, el, r)
    real(wp), intent(in) :: east, north, up
    real(wp), intent(out) :: az, el, r
    real(wp) :: dist
    dist = sqrt(east**2 + north**2 + up**2)  ! Compute Euclidean distance in ENU coordinates
    if (dist == 0.0_wp) then
      az = 0.0_wp  ! Default azimuth if zero distance
      el = 0.0_wp  ! Default elevation if zero distance
      r = 0.0_wp   ! Default range if zero distance
    else
      az = atan2(east, north)  ! Compute azimuth from east and north
      if (az < 0.0_wp) az = az + two_pi  ! Adjust azimuth to be positive
      el = asin(up/dist)  ! Compute elevation angle
      r = dist  ! Set range equal to computed distance
    end if
  end subroutine enu_to_aer

  !----------------------------------------------------------------------
  ! Subroutine: aer_to_enu
  ! Purpose: Converts AER coordinates (azimuth, elevation, range) to ENU coordinates (east, north, up)
  !----------------------------------------------------------------------
  elemental subroutine aer_to_enu(az, el, r, east, north, up)
    real(wp), intent(in) :: az, el, r
    real(wp), intent(out) :: east, north, up
    east = r * cos(el) * sin(az)  ! Convert AER to east component
    north = r * cos(el) * cos(az)  ! Convert to north component
    up = r * sin(el)  ! Convert to up component
  end subroutine aer_to_enu

  !----------------------------------------------------------------------
  ! Subroutine: ecef_to_aer
  ! Purpose: Converts ECEF coordinates (x, y, z) to AER coordinates (azimuth, elevation, range)
  !          via an intermediate conversion to ENU.
  !----------------------------------------------------------------------
  elemental subroutine ecef_to_aer(x, y, z, ref_lat, ref_lon, ref_h, az, el, r)
    real(wp), intent(in) :: x, y, z
    real(wp), intent(in) :: ref_lat, ref_lon, ref_h
    real(wp), intent(out) :: az, el, r
    real(wp) :: east, north, up
    call ecef_to_enu(x, y, z, ref_lat, ref_lon, ref_h, east, north, up)  ! Convert from ECEF to ENU
    call enu_to_aer(east, north, up, az, el, r)  ! Convert from ENU to AER
  end subroutine ecef_to_aer

  !----------------------------------------------------------------------
  ! Subroutine: aer_to_ecef
  ! Purpose: Converts AER coordinates (azimuth, elevation, range) back to ECEF coordinates (x, y, z)
  !          via an intermediate conversion from AER to ENU.
  !----------------------------------------------------------------------
  elemental subroutine aer_to_ecef(az, el, r, ref_lat, ref_lon, ref_h, x, y, z)
    real(wp), intent(in) :: az, el, r
    real(wp), intent(in) :: ref_lat, ref_lon, ref_h
    real(wp), intent(out) :: x, y, z
    real(wp) :: east, north, up
    call aer_to_enu(az, el, r, east, north, up)  ! Convert from AER to ENU
    call enu_to_ecef(east, north, up, ref_lat, ref_lon, ref_h, x, y, z)  ! Convert from ENU to ECEF
  end subroutine aer_to_ecef

  !----------------------------------------------------------------------
  ! Subroutine: ECEFspher2ENU
  ! Purpose: Converts ECEF spherical coordinates (altitude, zenith, azimuth) into ENU coordinates
  !          using spherical Earth assumptions.
  !----------------------------------------------------------------------
  ! elemental subroutine ECEFspher2ENU(alt, theta, phi, theta1, phi1, x, y, z)
  !   real(wp), intent(in) :: alt, theta, phi  ! Spherical coordinates: altitude, zenith, and azimuth
  !   real(wp), intent(in) :: theta1, phi1  ! Reference spherical coordinates for ENU conversion
  !   real(wp), intent(out) :: x, y, z  ! Resulting ENU coordinates
  !   real(wp) :: gamma1, gamma2, xp, yp
  !   gamma1 = acos( cos(theta)*cos(theta1) + sin(theta)*sin(theta1)*cos(phi-phi1) )  ! Compute angular distance gamma1
  !   gamma2 = acos( cos(theta1)*cos(theta) + sin(theta1)*sin(theta)*cos(phi1-phi) )  ! Compute angular distance gamma2
  !   xp = Re * gamma1  ! Convert angular distance to linear displacement (x-direction)
  !   yp = Re * gamma2  ! Convert angular distance to linear displacement (y-direction)
  !   if (theta > theta1) then
  !     yp = -yp  ! Adjust sign of y displacement based on relative zenith angle
  !   end if
  !   if (phi < phi1) then
  !     xp = -xp  ! Adjust sign of x displacement based on relative azimuth
  !   end if
  !   x = xp
  !   y = yp
  !   z = alt  ! Altitude remains unchanged
  ! end subroutine ECEFspher2ENU

  !----------------------------------------------------------------------
  ! Subroutine: geog2geomag
  ! Purpose: Converts geographic coordinates (longitude, latitude in degrees)
  !          to geomagnetic coordinates (phi, theta in radians).
  !----------------------------------------------------------------------
  elemental subroutine geog2geomag(glon, glat, phi, theta)
    ! Converts geographic coordinates (glon, glat in degrees) to geomagnetic coordinates (phi, theta in radians).
    real(wp), intent(in) :: glon, glat
    real(wp), intent(out) :: phi, theta
    real(wp) :: glonwrap, thetag, argtmp, alpha
    glonwrap = mod(glon, 360._wp)  ! Wrap longitude to [0,360)
    thetag = pi/2 - glat*pi/180.0_wp  ! Convert geodetic latitude to zenith angle
    glonwrap = glonwrap * pi/180.0_wp  ! Convert longitude to radians
    theta = acos(cos(thetag)*cos(thetan) + sin(thetag)*sin(thetan)*cos(glonwrap - phin))  ! Compute angular difference
    argtmp = (cos(thetag) - cos(theta)*cos(thetan)) / (sin(theta)*sin(thetan))  ! Intermediate value for conversion
    alpha = acos( max( min(argtmp, 1._wp), -1._wp) )  ! Clamp and compute alpha
    if ((phin > glonwrap .and. phin - glonwrap > pi) .or. (phin < glonwrap .and. glonwrap - phin < pi)) then
      phi = pi - alpha  ! Adjust phi based on geomagnetic pole position
    else
      phi = alpha + pi  ! Adjust phi for alternate configuration
    end if
  end subroutine geog2geomag

  !----------------------------------------------------------------------
  ! Subroutine: geomag2geog
  ! Purpose: Converts geomagnetic coordinates (phi, theta in radians)
  !          back to geographic coordinates (longitude, latitude in degrees).
  !----------------------------------------------------------------------
  elemental subroutine geomag2geog(phi, theta, glon, glat)
    ! Converts geomagnetic coordinates (phi, theta in radians) to geographic coordinates (glon, glat in degrees).
    real(wp), intent(in) :: phi, theta
    real(wp), intent(out) :: glon, glat
    real(wp) :: phiwrap, thetag2p, beta, phig2, thetag2, argtmp
    phiwrap = mod(phi, 2.0_wp*pi)  ! Wrap phi to [0,2*pi)
    thetag2p = acos(cos(theta)*cos(thetan) - sin(theta)*sin(thetan)*cos(phiwrap))  ! Intermediate zenith angle
    argtmp = (cos(theta) - cos(thetag2p)*cos(thetan)) / (sin(thetag2p)*sin(thetan))  ! Intermediate value for beta
    beta = acos( max( min(argtmp, 1._wp), -1._wp) )  ! Clamp and compute beta
    if (phiwrap > pi) then
      phig2 = phin - beta  ! Adjust geographic longitude for phiwrap > pi
    else
      phig2 = phin + beta  ! Adjust geographic longitude for phiwrap <= pi
    end if
    phig2 = mod(phig2, 2.0_wp*pi)  ! Wrap geographic longitude to [0,2*pi)
    thetag2 = pi/2 - thetag2p  ! Convert zenith angle to geographic latitude
    glon = phig2 * 180.0_wp/pi  ! Convert longitude from radians to degrees
    glat = thetag2 * 180.0_wp/pi  ! Convert latitude from radians to degrees
  end subroutine geomag2geog

  !----------------------------------------------------------------------
  ! Function: r2alt
  ! Purpose: Converts a geocentric distance (r) into an altitude (alt)
  !----------------------------------------------------------------------
  elemental function r2alt(r) result(alt)
    real(wp), intent(in) :: r
    real(wp) :: alt
    alt = r - Re  ! Altitude is the difference between geocentric distance and Earth radius
  end function r2alt

  !----------------------------------------------------------------------
  ! Function: alt2r
  ! Purpose: Converts an altitude (alt) into a geocentric distance (r)
  !----------------------------------------------------------------------
  elemental function alt2r(alt) result(r)
    real(wp), intent(in) :: alt
    real(wp) :: r
    r = alt + Re  ! Geocentric distance is altitude plus Earth radius
  end function alt2r

  !----------------------------------------------------------------------
  ! Function: rotz
  ! Purpose: Constructs a rotation matrix for a rotation about the z-axis.
  !----------------------------------------------------------------------
  function rotz(alpha) result(Rz)
    real(wp), intent(in) :: alpha
    real(wp), dimension(3,3) :: Rz
    Rz = 0._wp  ! Initialize rotation matrix to zero
    Rz(1,1) = cos(alpha)  ! Rotation about z-axis (element 1,1)
    Rz(1,2) = -sin(alpha)  ! Rotation about z-axis (element 1,2)
    Rz(2,1) = sin(alpha)  ! Rotation about z-axis (element 2,1)
    Rz(2,2) = cos(alpha)  ! Rotation about z-axis (element 2,2)
    Rz(3,3) = 1._wp  ! No rotation for z-component
  end function rotz

  !----------------------------------------------------------------------
  ! Function: roty
  ! Purpose: Constructs a rotation matrix for a rotation about the y-axis.
  !----------------------------------------------------------------------
  function roty(alpha) result(Ry)
    real(wp), intent(in) :: alpha
    real(wp), dimension(3,3) :: Ry
    Ry = 0._wp  ! Initialize rotation matrix to zero
    Ry(1,1) = cos(alpha)  ! Rotation about y-axis (element 1,1)
    Ry(1,3) = sin(alpha)  ! Rotation about y-axis (element 1,3)
    Ry(2,2) = 1._wp  ! y-axis remains unchanged
    Ry(3,1) = -sin(alpha)  ! Rotation about y-axis (element 3,1)
    Ry(3,3) = cos(alpha)  ! Rotation about y-axis (element 3,3)
  end function roty

  !----------------------------------------------------------------------
  ! Function: rotgm2gg
  ! Purpose: Constructs a rotation matrix to convert from geomagnetic ECEF
  !          to geographic ECEF.
  !----------------------------------------------------------------------
  function rotgm2gg() result(Rgm2gg)
    real(wp), dimension(3,3) :: Rgm2gg, Rz, Ry
    Rz = rotz(phin)  ! Rotation about z by geomagnetic longitude
    Ry = roty(thetan)  ! Rotation about y by geomagnetic colatitude
    Rgm2gg = matmul(Rz, Ry)  ! Combined rotation from geomagnetic to geographic
  end function rotgm2gg

  !----------------------------------------------------------------------
  ! Function: rotgg2gm
  ! Purpose: Constructs a rotation matrix to convert from geographic ECEF
  !          to geomagnetic ECEF.
  !----------------------------------------------------------------------
  function rotgg2gm() result(Rgg2gm)
    real(wp), dimension(3,3) :: Rgg2gm, Rz, Ry
    Rz = rotz(phin)  ! Rotation about z for reverse conversion
    Ry = roty(thetan)  ! Rotation about y for reverse conversion
    Rgg2gm = matmul(transpose(Ry), transpose(Rz))  ! Combined reverse rotation from geographic to geomagnetic
  end function rotgg2gm

end module coord_trans