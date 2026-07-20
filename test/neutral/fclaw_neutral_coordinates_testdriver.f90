program fclaw_neutral_coordinates_testdriver

use phys_consts, only : wp,pi
use gemini3d_config, only : gemini_cfg,neutral_coord_geomagnetic_spherical, &
  neutral_coord_geographic_loximuthal
use geomagnetic, only : ECEFspher2ENU,geog2loximuth
use neutraldata3Dobj_fclaw, only : set_fclaw_neutral_query_positions

implicit none (type, external)

integer, parameter :: n1=2,n2=2,n3=1
real(wp), parameter :: tol_m=1.e-8_wp
type(gemini_cfg) :: cfg
real(wp) :: alt(n1,n2,n3),theta(n1,n2,n3),phi(n1,n2,n3)
real(wp) :: glat(n1,n2,n3),glon(n1,n2,n3)
real(wp) :: x_m(n1,n2,n3),y_m(n1,n2,n3),z_m(n1,n2,n3)
real(wp) :: x_expected(n1,n2,n3),y_expected(n1,n2,n3),z_expected(n1,n2,n3)
real(wp) :: source_theta,source_phi

alt=reshape([12000._wp,87000._wp,250000._wp,400000._wp],shape(alt))
glat=reshape([33.3444_wp,34.1_wp,31.8_wp,40.2_wp],shape(glat))
glon=reshape([-97.5702_wp,-96.2_wp,-100.1_wp,-92.7_wp],shape(glon))
theta=reshape([1.1_wp,1.2_wp,1.3_wp,1.4_wp],shape(theta))
phi=reshape([5.0_wp,5.1_wp,5.2_wp,5.3_wp],shape(phi))

cfg%neutral_coord_mode=neutral_coord_geographic_loximuthal
cfg%neutral_ref_glat=33.3444_wp
cfg%neutral_ref_glon=-97.5702_wp
call set_fclaw_neutral_query_positions(cfg,alt,theta,phi,glat,glon,x_m,y_m,z_m)
call geog2loximuth(glat,glon,cfg%neutral_ref_glat,cfg%neutral_ref_glon,x_expected,y_expected)
if (maxval(abs(x_m-x_expected))>tol_m .or. maxval(abs(y_m-y_expected))>tol_m) then
  error stop 'ForestClaw geographic loximuthal query coordinates differ from the projection library'
end if
if (maxval(abs(z_m-alt))>tol_m) error stop 'ForestClaw geographic loximuthal z is not GEMINI altitude'

cfg%neutral_coord_mode=neutral_coord_geomagnetic_spherical
cfg%sourcemlat=43.01_wp
cfg%sourcemlon=329.27_wp
source_theta=pi/2._wp-cfg%sourcemlat*pi/180._wp
source_phi=cfg%sourcemlon*pi/180._wp
call set_fclaw_neutral_query_positions(cfg,alt,theta,phi,glat,glon,x_m,y_m,z_m)
call ECEFspher2ENU(alt,theta,phi,source_theta,source_phi,x_expected,y_expected,z_expected)
if (maxval(abs(x_m-x_expected))>tol_m .or. maxval(abs(y_m-y_expected))>tol_m .or. &
    maxval(abs(z_m-z_expected))>tol_m) then
  error stop 'ForestClaw legacy query coordinate path changed'
end if

print '(a)', 'ForestClaw neutral coordinate mapping tests: PASS'

end program fclaw_neutral_coordinates_testdriver
