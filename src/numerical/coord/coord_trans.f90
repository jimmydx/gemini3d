module coord_trans
  implicit none
  public :: set_earth_model, geodetic_to_ecef, ecef_to_geodetic, ecef_to_enu, enu_to_ecef, &
            enu_to_aer, aer_to_enu, ecef_to_aer, aer_to_ecef, ECEFspher2ENU, &
            geog2geomag, geomag2geog, r2alt, alt2r, rotz, roty, rotgm2gg, rotgg2gm

  integer, parameter :: wp = selected_real_kind(15, 307)

  !-----------------------------------------------------------------------
  ! WGS84 ellipsoid parameters
  !-----------------------------------------------------------------------
  real(wp), parameter :: a  = 6378137.0_wp
  real(wp), parameter :: f  = 1.0_wp/298.257223563_wp
  real(wp), parameter :: b  = a*(1.0_wp - f)
  real(wp), parameter :: e2 = (a**2 - b**2)/a**2
  real(wp), parameter :: ep2= (a**2 - b**2)/b**2

  ! Spherical Earth parameter
  real(wp), parameter :: Re = 6371000.0_wp

  ! Earth model types
  integer, parameter :: SPHERICAL = 0, WGS84 = 1
  integer :: earth_model = WGS84

  ! Fundamental constants
  real(wp), parameter :: pi = 3.14159265358979323846_wp
  real(wp), parameter :: two_pi = 2.0_wp*pi

  !-----------------------------------------------------------------------
  ! Geomagnetic parameters
  !-----------------------------------------------------------------------
  real(wp), parameter :: thetan = 11.0_wp*pi/180.0_wp
  real(wp), parameter :: phin   = 289.0_wp*pi/180.0_wp

contains

  subroutine set_earth_model(model_choice)
    integer, intent(in) :: model_choice
    if (model_choice == SPHERICAL .or. model_choice == WGS84) then
      earth_model = model_choice
    else
      print*, "Invalid earth model choice. Defaulting to SPHERICAL."
      earth_model = SPHERICAL
    end if
  end subroutine set_earth_model

  elemental subroutine geodetic_to_ecef(lat, lon, h, x, y, z)
    real(wp), intent(in) :: lat, lon, h
    real(wp), intent(out) :: x, y, z
    if (earth_model == WGS84) then
      x = ( a / sqrt(1.0_wp - e2*sin(lat)**2) + h ) * cos(lat) * cos(lon)
      y = ( a / sqrt(1.0_wp - e2*sin(lat)**2) + h ) * cos(lat) * sin(lon)
      z = ( a / sqrt(1.0_wp - e2*sin(lat)**2)*(1.0_wp - e2) + h ) * sin(lat)
    else
      x = ( Re + h ) * cos(lat) * cos(lon)
      y = ( Re + h ) * cos(lat) * sin(lon)
      z = ( Re + h ) * sin(lat)
    end if
  end subroutine geodetic_to_ecef

  elemental subroutine ecef_to_geodetic(x, y, z, lat, lon, h)
    real(wp), intent(in) :: x, y, z
    real(wp), intent(out) :: lat, lon, h
    real(wp) :: p, theta, N_local
    if (earth_model == WGS84) then
      p = sqrt(x**2 + y**2)
      lon = atan2(y, x)
      if (p == 0.0_wp) then
        if (z >= 0.0_wp) then
          lat = pi/2.0_wp
        else
          lat = -pi/2.0_wp
        end if
        h = abs(z) - b
      else
        theta = atan2(z*a, p*b)
        lat = atan2(z + ep2*b*sin(theta)**3, p - e2*a*cos(theta)**3)
        N_local = a / sqrt(1.0_wp - e2*sin(lat)**2)
        h = p/cos(lat) - N_local
      end if
    else
      p = sqrt(x**2 + y**2)
      lon = atan2(y, x)
      if (p == 0.0_wp) then
        if (z >= 0.0_wp) then
          lat = pi/2.0_wp
        else
          lat = -pi/2.0_wp
        end if
        h = abs(z) - Re
      else
        lat = atan2(z, p)
        h = sqrt(x**2 + y**2 + z**2) - Re
      end if
    end if
  end subroutine ecef_to_geodetic

  elemental subroutine ecef_to_enu(x, y, z, ref_lat, ref_lon, ref_h, east, north, up)
    real(wp), intent(in) :: x, y, z
    real(wp), intent(in) :: ref_lat, ref_lon, ref_h
    real(wp), intent(out) :: east, north, up
    real(wp) :: dx, dy, dz
    dx = x - ((a / sqrt(1.0_wp - e2*sin(ref_lat)**2) + ref_h)*cos(ref_lat)*cos(ref_lon))
    dy = y - ((a / sqrt(1.0_wp - e2*sin(ref_lat)**2) + ref_h)*cos(ref_lat)*sin(ref_lon))
    dz = z - (((a*(1.0_wp - e2)) / sqrt(1.0_wp - e2*sin(ref_lat)**2) + ref_h)*sin(ref_lat))
    east = -sin(ref_lon)*dx + cos(ref_lon)*dy
    north = -sin(ref_lat)*cos(ref_lon)*dx - sin(ref_lat)*sin(ref_lon)*dy + cos(ref_lat)*dz
    up = cos(ref_lat)*cos(ref_lon)*dx + cos(ref_lat)*sin(ref_lon)*dy + sin(ref_lat)*dz
  end subroutine ecef_to_enu

  elemental subroutine enu_to_ecef(east, north, up, ref_lat, ref_lon, ref_h, x, y, z)
    real(wp), intent(in) :: east, north, up
    real(wp), intent(in) :: ref_lat, ref_lon, ref_h
    real(wp), intent(out) :: x, y, z
    real(wp) :: dx, dy, dz
    dx = -sin(ref_lon)*east - sin(ref_lat)*cos(ref_lon)*north + cos(ref_lat)*cos(ref_lon)*up
    dy = cos(ref_lon)*east - sin(ref_lat)*sin(ref_lon)*north + cos(ref_lat)*sin(ref_lon)*up
    dz = cos(ref_lat)*north + sin(ref_lat)*up
    x = ((a / sqrt(1.0_wp - e2*sin(ref_lat)**2) + ref_h)*cos(ref_lat)*cos(ref_lon)) + dx
    y = ((a / sqrt(1.0_wp - e2*sin(ref_lat)**2) + ref_h)*cos(ref_lat)*sin(ref_lon)) + dy
    z = (((a*(1.0_wp - e2)) / sqrt(1.0_wp - e2*sin(ref_lat)**2) + ref_h)*sin(ref_lat)) + dz
  end subroutine enu_to_ecef

  elemental subroutine enu_to_aer(east, north, up, az, el, r)
    real(wp), intent(in) :: east, north, up
    real(wp), intent(out) :: az, el, r
    real(wp) :: dist
    dist = sqrt(east**2 + north**2 + up**2)
    if (dist == 0.0_wp) then
      az = 0.0_wp
      el = 0.0_wp
      r = 0.0_wp
    else
      az = atan2(east, north)
      if (az < 0.0_wp) az = az + two_pi
      el = asin(up/dist)
      r = dist
    end if
  end subroutine enu_to_aer

  elemental subroutine aer_to_enu(az, el, r, east, north, up)
    real(wp), intent(in) :: az, el, r
    real(wp), intent(out) :: east, north, up
    east = r * cos(el) * sin(az)
    north = r * cos(el) * cos(az)
    up = r * sin(el)
  end subroutine aer_to_enu

  elemental subroutine ecef_to_aer(x, y, z, ref_lat, ref_lon, ref_h, az, el, r)
    real(wp), intent(in) :: x, y, z
    real(wp), intent(in) :: ref_lat, ref_lon, ref_h
    real(wp), intent(out) :: az, el, r
    real(wp) :: east, north, up
    call ecef_to_enu(x, y, z, ref_lat, ref_lon, ref_h, east, north, up)
    call enu_to_aer(east, north, up, az, el, r)
  end subroutine ecef_to_aer

  elemental subroutine aer_to_ecef(az, el, r, ref_lat, ref_lon, ref_h, x, y, z)
    real(wp), intent(in) :: az, el, r
    real(wp), intent(in) :: ref_lat, ref_lon, ref_h
    real(wp), intent(out) :: x, y, z
    real(wp) :: east, north, up
    call aer_to_enu(az, el, r, east, north, up)
    call enu_to_ecef(east, north, up, ref_lat, ref_lon, ref_h, x, y, z)
  end subroutine aer_to_ecef

  elemental subroutine ECEFspher2ENU(alt, theta, phi, theta1, phi1, x, y, z)
    real(wp), intent(in) :: alt, theta, phi
    real(wp), intent(in) :: theta1, phi1
    real(wp), intent(out) :: x, y, z
    real(wp) :: gamma1, gamma2, xp, yp
    gamma1 = acos( cos(theta)*cos(theta1) + sin(theta)*sin(theta1)*cos(phi-phi1) )
    gamma2 = acos( cos(theta1)*cos(theta) + sin(theta1)*sin(theta)*cos(phi1-phi) )
    xp = Re * gamma1
    yp = Re * gamma2
    if (theta > theta1) then
      yp = -yp
    end if
    if (phi < phi1) then
      xp = -xp
    end if
    x = xp
    y = yp
    z = alt
  end subroutine ECEFspher2ENU

  elemental subroutine geog2geomag(glon, glat, phi, theta)
    ! Converts geographic coordinates (glon, glat in degrees) to geomagnetic
    ! coordinates (phi, theta in radians).
    real(wp), intent(in) :: glon, glat
    real(wp), intent(out) :: phi, theta
    real(wp) :: glonwrap, thetag, argtmp, alpha
    glonwrap = mod(glon, 360._wp)
    thetag = pi/2 - glat*pi/180.0_wp
    glonwrap = glonwrap * pi/180.0_wp
    theta = acos(cos(thetag)*cos(thetan) + sin(thetag)*sin(thetan)*cos(glonwrap - phin))
    argtmp = (cos(thetag) - cos(theta)*cos(thetan)) / (sin(theta)*sin(thetan))
    alpha = acos( max( min(argtmp, 1._wp), -1._wp) )
    if ((phin > glonwrap .and. phin - glonwrap > pi) .or. (phin < glonwrap .and. glonwrap - phin < pi)) then
      phi = pi - alpha
    else
      phi = alpha + pi
    end if
  end subroutine geog2geomag

  elemental subroutine geomag2geog(phi, theta, glon, glat)
    ! Converts geomagnetic coordinates (phi, theta in radians) to geographic
    ! coordinates (glon, glat in degrees).
    real(wp), intent(in) :: phi, theta
    real(wp), intent(out) :: glon, glat
    real(wp) :: phiwrap, thetag2p, beta, phig2, thetag2, argtmp
    phiwrap = mod(phi, 2.0_wp*pi)
    thetag2p = acos(cos(theta)*cos(thetan) - sin(theta)*sin(thetan)*cos(phiwrap))
    argtmp = (cos(theta) - cos(thetag2p)*cos(thetan)) / (sin(thetag2p)*sin(thetan))
    beta = acos( max( min(argtmp, 1._wp), -1._wp) )
    if (phiwrap > pi) then
      phig2 = phin - beta
    else
      phig2 = phin + beta
    end if
    phig2 = mod(phig2, 2.0_wp*pi)
    thetag2 = pi/2 - thetag2p
    glon = phig2 * 180.0_wp/pi
    glat = thetag2 * 180.0_wp/pi
  end subroutine geomag2geog

  elemental function r2alt(r) result(alt)
    real(wp), intent(in) :: r
    real(wp) :: alt
    alt = r - Re
  end function r2alt

  elemental function alt2r(alt) result(r)
    real(wp), intent(in) :: alt
    real(wp) :: r
    r = alt + Re
  end function alt2r

  function rotz(alpha) result(Rz)
    real(wp), intent(in) :: alpha
    real(wp), dimension(3,3) :: Rz
    Rz = 0._wp
    Rz(1,1) = cos(alpha)
    Rz(1,2) = -sin(alpha)
    Rz(2,1) = sin(alpha)
    Rz(2,2) = cos(alpha)
    Rz(3,3) = 1._wp
  end function rotz

  function roty(alpha) result(Ry)
    real(wp), intent(in) :: alpha
    real(wp), dimension(3,3) :: Ry
    Ry = 0._wp
    Ry(1,1) = cos(alpha)
    Ry(1,3) = sin(alpha)
    Ry(2,2) = 1._wp
    Ry(3,1) = -sin(alpha)
    Ry(3,3) = cos(alpha)
  end function roty

  function rotgm2gg() result(Rgm2gg)
    real(wp), dimension(3,3) :: Rgm2gg, Rz, Ry
    Rz = rotz(phin)
    Ry = roty(thetan)
    Rgm2gg = matmul(Rz, Ry)
  end function rotgm2gg

  function rotgg2gm() result(Rgg2gm)
    real(wp), dimension(3,3) :: Rgg2gm, Rz, Ry
    Rz = rotz(phin)
    Ry = roty(thetan)
    Rgg2gm = matmul(transpose(Ry), transpose(Rz))
  end function rotgg2gm

end module coord_trans