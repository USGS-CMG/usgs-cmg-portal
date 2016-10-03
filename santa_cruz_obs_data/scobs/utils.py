#!python
# coding=utf-8
import os
import json
import tempfile
from collections import namedtuple

import epic2cf
import numpy as np
from nco import Nco
import netCDF4 as nc4
from compliance_checker.runner import ComplianceChecker, CheckSuite

from scobs.mappings import (
    IGNORABLE_CODES,
    variable_name_overrides,
    long_name_overrides,
    special_map
)

import logging
logger = logging.getLogger(__name__)


def check_compliance(netcdf_file):
    check_suite = CheckSuite()
    check_suite.load_all_available_checkers()

    _, outfile = tempfile.mkstemp(prefix='scobs-cc-', suffix='.json')

    try:
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=netcdf_file,
            checker_names=['cf', 'acdd'],
            verbose=True,
            criteria='strict',
            output_format='json',
            output_filename=outfile
        )
        z = []
        score = namedtuple('ComplianceScore', 'name passed score possible')
        with open(outfile, 'rt') as f:
            for checker, results in json.load(f).items():
                z.append(
                    score(
                        name=checker,
                        passed=not errors,
                        score=results['scored_points'],
                        possible=results['possible_points']
                    )
                )
        return z

    except BaseException as e:
        logger.warning(e)
    finally:
        if os.path.isfile(outfile):
            os.remove(outfile)


def normalize_units(netcdf_file):
    with nc4.Dataset(netcdf_file, 'a') as nc:
        for v in nc.variables:
            nc_var = nc.variables.get(v)
            if hasattr(nc_var, 'units') and nc_var.units == "K":
                # Convert kelvin to Celsius
                nc_var[:] = nc_var[:] - 273.15
                nc_var.units = "degree_Celsius"
            elif hasattr(nc_var, 'standard_name') and nc_var.standard_name == 'sea_surface_wave_from_direction':
                # Convert "From" to "To" direction
                nc_var[:] = (nc_var[:] + 180) % 360
                nc_var.standard_name = 'sea_surface_wave_to_direction'
                nc_var.long_name = "Wave Direction (to TN)"


def normalize_time(netcdf_file):
    epoch_units       = 'seconds since 1970-01-01T00:00:00Z'
    millisecond_units = 'milliseconds since 1858-11-17T00:00:00Z'

    with nc4.Dataset(netcdf_file, 'a') as nc:
        # Signell said this works, any problems and we can all blame him!
        time_data = nc4.num2date(
            (
                np.int64(nc.variables['time'][:]) - 2400001
            ) * 3600 * 24 * 1000 + nc.variables['time2'][:].__array__(),
            units=millisecond_units
        )
        nc.renameVariable("time", "old_time")
        nc.sync()

        time = nc.createVariable('time', 'f8', ('time'))
        time.units          = epoch_units
        time.standard_name  = "time"
        time.long_name      = "time of measurement"
        time.calendar       = "gregorian"
        time.axis           = "T"
        time[:] = nc4.date2num(time_data, units=epoch_units).round()

    o = Nco()
    o.ncks(
        input=netcdf_file,
        output=netcdf_file,
        options=[
            '-O',
            '-h',
            '-x',
            '-v', 'time2,old_time'
        ]
    )


def normalize_epic_codes(netcdf_file):
    with nc4.Dataset(netcdf_file, 'a') as nc:
        for v in nc.variables:
            nc_var = nc.variables.get(v)
            if v in variable_name_overrides:
                ec = variable_name_overrides.get(v).get('epic_code', None)
                if ec is not None:
                    nc_var.epic_code = ec
                overrides = variable_name_overrides.get(v).get('overrides', dict())
                for k, d in overrides.items():
                    if k == 'convert':
                        nc_var[:] = d(nc_var[:])
                    elif k != 'original_units':
                        nc_var.setncattr(k, d)

            if hasattr(nc_var, 'long_name'):
                if not hasattr(nc_var, 'epic_code') or (hasattr(nc_var, 'epic_code') and nc_var.epic_code in IGNORABLE_CODES):
                    lookup_long_name = nc_var.long_name.lower().strip()
                    if lookup_long_name in long_name_overrides:
                        ec = long_name_overrides.get(lookup_long_name).get('epic_code', None)
                        if ec is not None:
                            nc_var.epic_code = ec
                        overrides = long_name_overrides.get(lookup_long_name).get('overrides', dict())
                        for k, d in overrides.items():
                            if k == 'convert':
                                nc_var[:] = d(nc_var[:])
                            elif k != 'original_units':
                                nc_var.setncattr(k, d)

            if hasattr(nc_var, "epic_code") and nc_var.epic_code:
                try:
                    int(nc_var.epic_code)
                except ValueError:
                    logger.debug("No EPIC code specified on {0}".format(v))
                else:

                    # Specialized cases for generic EPIC codes
                    if nc_var.epic_code in special_map:
                        attribs = special_map.get(int(nc_var.epic_code))(nc_var)
                    else:
                        attribs = epic2cf.mapping.get(int(nc_var.epic_code))

                    # Special case for 'Onset weather stations'.
                    # https://github.com/USGS-CMG/usgs-cmg-portal/issues/69
                    if int(nc_var.epic_code) in [905, 908] and 'hml' in netcdf_file.lower():
                        attribs.standard_name = 'surface_downwelling_photosynthetic_radiative_flux_in_air'

                    if attribs is not None and attribs.standard_name is not None:
                        # Convert data to CF units
                        nc_var[:] = attribs.convert(nc_var[:])
                        # Set attributes
                        nc_var.standard_name = attribs.standard_name
                        nc_var.long_name     = attribs.long_name
                        nc_var.units         = attribs.cf_units
                        if attribs.cell_methods is not None:
                            nc_var.cell_methods = attribs.cell_methods
                    else:
                        logger.debug("Could not find CF mapping for EPIC code {!s}".format(nc_var.epic_code))


def normalize_vectors(netcdf_file):
    with nc4.Dataset(netcdf_file, 'a') as nc:
        east  = None
        north = None
        for v in nc.variables:
            nc_var = nc.variables.get(v)
            if hasattr(nc_var, 'standard_name') and nc_var.standard_name == 'eastward_sea_water_velocity':
                east = nc_var
                continue
            if hasattr(nc_var, 'standard_name') and nc_var.standard_name == 'northward_sea_water_velocity':
                north = nc_var
                continue

        std_names = []
        for varname in nc.variables:
            var = nc.variables.get(varname)
            if hasattr(var, 'standard_name'):
                std_names.append(var.standard_name)

        # Only add the variables if they don't already exist
        if east is not None and north is not None and 'sea_water_speed' not in std_names and 'direction_of_sea_water_velocity' not in std_names:
            # We have vectors... create the speed and direction variables
            # astype avoids overflow errors on 32 bit floats
            speed = np.sqrt(np.square(east[:].astype(np.float64)) + np.square(north[:].astype(np.float64)))
            direction = np.degrees(np.arctan2(north[:], east[:]))

            east_fill_value = east._FillValue if hasattr(east, '_FillValue') else np.nan
            spd = nc.createVariable('CS_300', east.dtype, east.dimensions, fill_value=east_fill_value)
            spd.standard_name = 'sea_water_speed'
            spd.long_name = "Current speed"
            spd.units = 'm/s'
            spd.epic_code = 300
            spd[:] = speed

            drc = nc.createVariable('CD_310', east.dtype, east.dimensions, fill_value=east_fill_value)
            drc.standard_name = 'direction_of_sea_water_velocity'
            drc.long_name = "Current direction"
            drc.units = 'degree'
            drc.epic_code = 310
            drc[:] = direction


def normalize_depth_direction(netcdf_file):
    with nc4.Dataset(netcdf_file, 'a') as nc:
        # Get all depth variables
        depth_variables = []
        for dv in nc.variables:
            depth_variables += [ x for x in nc.variables.get(dv).dimensions if 'depth' in x ]
        depth_variables = sorted(list(set(depth_variables)))

        # Convert everything to positive up, unless it is specifically specified as "up" already
        for dv in depth_variables:
            dvar = nc.variables.get(dv)
            if hasattr(dvar, 'positive') and dvar.positive.lower() == 'up':
                pass
            else:
                dvar[:] = dvar[:] * -1.0


def normalize_netcdf4(netcdf_file):
    o = Nco()
    o.ncks(
        input=netcdf_file,
        output=netcdf_file,
        options=[
            '-O',
            '-h',
            '-4',
            '-L3'
        ]
    )
