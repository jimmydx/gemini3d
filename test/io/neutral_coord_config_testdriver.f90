program neutral_coord_config_testdriver

use phys_consts, only : wp
use gemini3d_config, only : read_configfile,gemini_cfg, &
  neutral_coord_geomagnetic_spherical,neutral_coord_geographic_loximuthal

implicit none (type, external)

type(gemini_cfg) :: cfg
character(1024) :: filename,expected_mode
integer :: stat

call get_command_argument(1,filename,status=stat)
if (stat/=0 .or. len_trim(filename)==0) error stop 'missing neutral coordinate test namelist'

cfg%infile=trim(filename)
cfg%outdir='.'
call read_configfile(cfg)

call get_command_argument(2,expected_mode,status=stat)
if (stat/=0 .or. len_trim(expected_mode)==0) then
  error stop 'invalid neutral coordinate configuration unexpectedly passed validation'
end if

if (trim(cfg%neutral_coord_mode)/=trim(expected_mode)) then
  error stop 'neutral coordinate mode was not parsed as expected'
end if

select case (trim(expected_mode))
case (neutral_coord_geomagnetic_spherical)
  if (cfg%neutral_ref_glat/=999._wp .or. cfg%neutral_ref_glon/=999._wp) then
    error stop 'legacy neutral coordinate defaults changed'
  end if
case (neutral_coord_geographic_loximuthal)
  if (abs(cfg%neutral_ref_glat-33.3444_wp)>epsilon(1._wp) .or. &
      abs(cfg%neutral_ref_glon+97.5702_wp)>epsilon(1._wp)) then
    error stop 'loximuthal geographic reference was not parsed correctly'
  end if
case default
  error stop 'unexpected neutral coordinate test mode'
end select

print '(a,a)', 'Neutral coordinate configuration test passed: ',trim(expected_mode)

end program neutral_coord_config_testdriver
