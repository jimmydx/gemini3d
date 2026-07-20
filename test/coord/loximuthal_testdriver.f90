program loximuthal_testdriver

use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
use phys_consts, only : wp,Re,pi
use geomagnetic, only : geog2loximuth,loximuth2geog

implicit none (type, external)

integer, parameter :: n_axis=17,n_golden=9
real(wp), parameter :: ref_glat_deg=33.3444_wp
real(wp), parameter :: ref_glon_deg=-97.5702_wp
real(wp), parameter :: half_width_m=960000._wp
real(wp), parameter :: roundtrip_limit_m=1.e-6_wp
real(wp), parameter :: golden_limit_m=1._wp
real(wp), parameter :: golden_x_m(n_golden)=[ &
  -960000._wp,0._wp,960000._wp,-960000._wp,0._wp,960000._wp,-960000._wp,0._wp,960000._wp]
real(wp), parameter :: golden_y_m(n_golden)=[ &
  -960000._wp,-960000._wp,-960000._wp,0._wp,0._wp,0._wp,960000._wp,960000._wp,960000._wp]
real(wp), parameter :: golden_glat_deg(n_golden)=[ &
  24.7109125831802_wp,24.7109125831802_wp,24.7109125831802_wp, &
  33.344401_wp,33.344401_wp,33.344401_wp, &
  41.9778874168198_wp,41.9778874168198_wp,41.9778874168198_wp]
real(wp), parameter :: golden_glon_deg(n_golden)=[ &
  -107.459112735974_wp,-97.5702_wp,-87.6812872640256_wp, &
  -107.904982015618_wp,-97.5702_wp,-87.2354179843816_wp, &
  -108.498750054737_wp,-97.5702_wp,-86.6416499452627_wp]

integer :: i,j,k
real(wp) :: x_m,y_m,glat_deg,glon_deg,x_back_m,y_back_m
real(wp) :: center_lat,center_lon,east_lat,east_lon,north_lat,north_lon
real(wp) :: max_roundtrip_error_m,max_golden_error_m,golden_error_m
real(wp) :: x_wrap_a,x_wrap_b,y_dummy,expected_x_m,nan_a,nan_b

max_roundtrip_error_m=0._wp
do j=1,n_axis
  y_m=-half_width_m+2._wp*half_width_m*real(j-1,wp)/real(n_axis-1,wp)
  do i=1,n_axis
    x_m=-half_width_m+2._wp*half_width_m*real(i-1,wp)/real(n_axis-1,wp)
    call loximuth2geog(x_m,y_m,ref_glat_deg,ref_glon_deg,glat_deg,glon_deg)
    call geog2loximuth(glat_deg,glon_deg,ref_glat_deg,ref_glon_deg,x_back_m,y_back_m)
    if (.not. ieee_is_finite(glat_deg) .or. .not. ieee_is_finite(glon_deg)) then
      error stop 'loximuthal inverse produced a non-finite coordinate inside the MAGIC domain'
    end if
    max_roundtrip_error_m=max(max_roundtrip_error_m,hypot(x_back_m-x_m,y_back_m-y_m))
  end do
end do

call loximuth2geog(0._wp,0._wp,ref_glat_deg,ref_glon_deg,center_lat,center_lon)
call loximuth2geog(1000._wp,0._wp,ref_glat_deg,ref_glon_deg,east_lat,east_lon)
call loximuth2geog(0._wp,1000._wp,ref_glat_deg,ref_glon_deg,north_lat,north_lon)
if (abs(center_lat-ref_glat_deg)>1.e-12_wp .or. abs(center_lon-ref_glon_deg)>1.e-12_wp) then
  error stop 'loximuthal origin does not map to the configured reference'
end if
if (east_lon<=center_lon .or. abs(east_lat-center_lat)>1.e-12_wp) then
  error stop 'positive loximuthal x does not point east along the reference parallel'
end if
if (north_lat<=center_lat .or. abs(north_lon-center_lon)>1.e-12_wp) then
  error stop 'positive loximuthal y does not point north along the reference meridian'
end if

call geog2loximuth(ref_glat_deg,ref_glon_deg+1._wp,ref_glat_deg,ref_glon_deg,x_m,y_m)
expected_x_m=Re*pi/180._wp*cos(ref_glat_deg*pi/180._wp)
if (abs(x_m-expected_x_m)>1.e-8_wp .or. abs(y_m)>1.e-12_wp) then
  error stop 'loximuthal central-parallel limit is incorrect'
end if

call geog2loximuth(ref_glat_deg,ref_glon_deg+181._wp,ref_glat_deg,ref_glon_deg,x_wrap_a,y_dummy)
call geog2loximuth(ref_glat_deg,ref_glon_deg-179._wp,ref_glat_deg,ref_glon_deg,x_wrap_b,y_dummy)
if (abs(x_wrap_a-x_wrap_b)>1.e-8_wp) error stop 'loximuthal longitude wrapping is inconsistent'

call geog2loximuth(90._wp,ref_glon_deg,ref_glat_deg,ref_glon_deg,nan_a,nan_b)
if (ieee_is_finite(nan_a) .or. ieee_is_finite(nan_b)) then
  error stop 'loximuthal forward projection accepted a geographic pole'
end if
call loximuth2geog(0._wp,Re*pi,ref_glat_deg,ref_glon_deg,nan_a,nan_b)
if (ieee_is_finite(nan_a) .or. ieee_is_finite(nan_b)) then
  error stop 'loximuthal inverse projection accepted an out-of-range latitude'
end if

max_golden_error_m=0._wp
do k=1,n_golden
  call loximuth2geog(golden_x_m(k),golden_y_m(k),ref_glat_deg,ref_glon_deg,glat_deg,glon_deg)
  golden_error_m=hypot( &
    (glat_deg-golden_glat_deg(k))*pi/180._wp*Re, &
    longitude_delta_deg(glon_deg,golden_glon_deg(k))*pi/180._wp*Re* &
      cos((glat_deg+golden_glat_deg(k))*pi/360._wp))
  max_golden_error_m=max(max_golden_error_m,golden_error_m)
end do

if (max_roundtrip_error_m>roundtrip_limit_m) then
  error stop 'loximuthal round-trip tolerance exceeded'
end if
if (max_golden_error_m>golden_limit_m) then
  error stop 'loximuthal MATLAB golden-point tolerance exceeded'
end if

print '(a,i0)', 'Loximuthal domain points tested: ',n_axis*n_axis
print '(a,es12.4,a)', 'Maximum round-trip error: ',max_roundtrip_error_m,' m'
print '(a,es12.4,a)', 'Maximum MATLAB golden-point error: ',max_golden_error_m,' m'
print '(a,f12.1,a)', 'Repository Earth radius Re: ',Re,' m'
print '(a)', 'Loximuthal coordinate tests: PASS'

contains

elemental real(wp) function longitude_delta_deg(a_deg,b_deg) result(delta_deg)
  real(wp), intent(in) :: a_deg,b_deg

  delta_deg=modulo(a_deg-b_deg+180._wp,360._wp)-180._wp
end function longitude_delta_deg

end program loximuthal_testdriver
