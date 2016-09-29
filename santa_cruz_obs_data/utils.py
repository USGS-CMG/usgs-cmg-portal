#!python
# coding=utf-8
from datetime import datetime
from types import SimpleNamespace

import epic2cf
import numpy as np
from nco import Nco
import netCDF4 as nc4

from epic2cf.data import location_codes, time_codes, generic_codes, voltage_codes

import logging
logger = logging.getLogger()
logger.addHandler(logging.NullHandler())

IGNORABLE_CODES = location_codes + time_codes + generic_codes + voltage_codes

# Special case EPIC mapping for generic EPIC codes that are used
special_map = {
    # '20' can be  Air Temperature or Water Temperature.
    20 : lambda y: epic2cf.mapping.get(32) if hasattr(y, 'sensor_depth') and y.sensor_depth > 0 else epic2cf.mapping.get(21),

    56 : lambda x: SimpleNamespace(
        standard_name='backscatter_intensity',
        long_name='Backscatter Intensity',
        units='v',
        convert=lambda x: x,
        cf_units='v',
        cell_methods=None
    ),
}

variable_name_overrides = {
    'w_1204min' : dict(epic_code=1204, overrides=dict(cell_methods='time: minimum')),
    'u_1205min' : dict(epic_code=1205, overrides=dict(cell_methods='time: minimum')),
    'v_1206min' : dict(epic_code=1206, overrides=dict(cell_methods='time: minimum')),
    'w_1204max' : dict(epic_code=1204, overrides=dict(cell_methods='time: maximum')),
    'u_1205max' : dict(epic_code=1205, overrides=dict(cell_methods='time: maximum')),
    'v_1206max' : dict(epic_code=1206, overrides=dict(cell_methods='time: maximum')),
    'WG_402'    : dict(epic_code=402),
    'Turb'      : dict(epic_code=980),
    'Press'     : dict(epic_code=1301),
    'vspd_1'    : dict(epic_code=300),
    'vdir_1'    : dict(epic_code=310),
    'vspd_2'    : dict(epic_code=300),
    'vdir_2'    : dict(epic_code=310),
    'u_1'       : dict(epic_code=1205),
    'v_1'       : dict(epic_code=1206),
    'w_1'       : dict(epic_code=1204),
    'u_2'       : dict(epic_code=1205),
    'v_2'       : dict(epic_code=1206),
    'w_2'       : dict(epic_code=1204),
    'bearing'   : dict(epic_code=1411),
    'rotor'     : dict(epic_code=4006),
    'DO'        : dict(epic_code=None, overrides=dict(standard_name='mass_concentration_of_oxygen_in_sea_water',
                                                      convert=lambda x: x / 1000.,
                                                      units='kg/m^3')),
    'BGAPE'     : dict(epic_code=None, overrides=dict(standard_name='mass_concentration_of_phycoerythrin_expressed_as_chlorophyll_in_sea_water',
                                                      units='kg/m^3',
                                                      convert=lambda x: x / 1000000.)),
    'turbidity' : dict(epic_code=980),
    'Qs_133'    : dict(epic_code=None, overrides=dict(standard_name='net_downward_shortwave_flux_in_air',
                                                      original_units='w/milli-angstrom^2',
                                                      units='w/m^2',
                                                      convert=lambda x: x / 1e13)),
    'BPR_1301'  : dict(epic_code=1301),
    'WL_NAVD88' : dict(epic_code=18, overrides=dict(vertical_datum='NAVD88')),
    'CTDCON_4218': dict(epic_code=4218),
}

long_name_overrides = {
    'salinity 1':                           dict(epic_code=40),
    'salinity 2':                           dict(epic_code=40),
    'salinity 2 q':                         dict(epic_code=40),
    'ctd salinity, pss-78':                 dict(epic_code=4214),
    'salinity':                             dict(epic_code=40),
    'salinity (ppt)':                       dict(epic_code=40),
    'salinity (psu)':                       dict(epic_code=41),
    'northward velocity':                   dict(epic_code=1206),
    'north':                                dict(epic_code=1206),
    'mean northward velocity':              dict(epic_code=1206),
    'north lp':                             dict(epic_code=1206),
    'eastward velocity':                    dict(epic_code=1205),
    'east':                                 dict(epic_code=1205),
    'mean eastward velocity':               dict(epic_code=1205),
    'east lp':                              dict(epic_code=1205),
    'instrument transducer temp.':          dict(epic_code=1211),
    'temperature (c)':                      dict(epic_code=32),
    'fr temp':                              dict(epic_code=32),
    'adp transducer temp.':                 dict(epic_code=1211),
    'adcp transducer temp.':                dict(epic_code=1211),
    'transducer temp.':                     dict(epic_code=1211),
    'temp 1':                               dict(epic_code=32),
    'temp 2':                               dict(epic_code=32),
    'temperature':                          dict(epic_code=32),
    'internal temperature':                 dict(epic_code=32),
    'frtemp':                               dict(epic_code=32),
    'temp 2 q':                             dict(epic_code=32),
    'temp':                                 dict(epic_code=32),
    'temp lp':                              dict(epic_code=32),
    'sea surface temperature (degrees C)':  dict(epic_code=36),
    'conductivity':                         dict(epic_code=50),
    'attenuation':                          dict(epic_code=55),
    'sigma theta':                          dict(epic_code=70),
    'psdev':                                dict(epic_code=850),
    'pressure':                             dict(epic_code=9),
    'sp cond':                              dict(epic_code=48),
    'dissolved oxygen saturation (mg/l)':   dict(epic_code=None, overrides=dict(standard_name='mass_concentration_of_oxygen_in_sea_water',
                                                                                original_units='mg/l',
                                                                                units='kg/m^3',
                                                                                convert=lambda x: x / 1000.)),
    'raw aanderaa dissolved oxygen concentration (um/kg)': dict(epic_code=65),
    'standard deviation of inst pitch':     dict(epic_code=1219),
    'standard deviation of inst roll':      dict(epic_code=1220),
    'cond':                                 dict(epic_code=51),
    'cond 1':                               dict(epic_code=51),
    'cond 2':                               dict(epic_code=51),
    'rotor count':                          dict(epic_code=4006),
    'rotor speed':                          dict(epic_code=4005),
    'compass':                              dict(epic_code=1401),
    'vane':                                 dict(epic_code=1402),
    'barometric pressure':                  dict(epic_code=None, overrides=dict(standard_name='air_pressure',
                                                                                long_name='Air Pressure',
                                                                                original_units='PSI',
                                                                                units='dbar',
                                                                                convert=lambda x: x * 0.6894757293)),
}

global_attributes = {
    'naming_authority':         'gov.usgs.cmgp',
    'source':                   'USGS',
    'institution':              'USGS Coastal and Marine Geology Program',
    'project':                  'U.S. Geological Survey Oceanographic Time-Series Data',
    'keywords':                 'Oceans > Ocean Pressure > Water Pressure, Oceans > Ocean Temperature > Water Temperature, Oceans > Salinity/Density > Conductivity, Oceans > Salinity/Density > Salinity',
    'keywords_vocabulary':      'GCMD Science Keywords',
    'standard_name_vocabulary': 'CF-1.6',
    'creator_email':            'rsignell@usgs.gov',
    'creator_name':             'Rich Signell',
    'creator_phone':            '+1 (508) 548-8700',
    'creator_url':              'http://www.usgs.gov',
    'publisher_email':          'emontgomery@usgs.gov',
    'publisher_name':           'Ellyn Montgomery',
    'publisher_phone':          '+1 (508) 548-8700',
    'publisher_url':            'http://www.usgs.gov',
    'contributor_role':         'principalInvestigator',
    'Conventions':              'CF-1.6',
    'date_created':             datetime.utcnow().strftime("%Y-%m-%dT%H:%M:00Z")
}

coord_vars = [
    'feature_type_instance',
    'time',
    'time2',
    'time_cf',
    'old_time',
    'depth',
    'depth002',
    'depth003',
    'depth004',
    'depth005',
    'lat',
    'lon'
]


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
        time[:] = nc4.date2num(time_data, units=epoch_units).round()
        return time_data[0]

    o = Nco()
    o.ncks(
        input=netcdf_file,
        output=netcdf_file,
        options=[
            '-O',
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


def normalize_ctd_location(netcdf_file):
    with nc4.Dataset(netcdf_file) as nc:
        try:
            latitude  = nc.variables.get("lat")[0]
            longitude = nc.variables.get("lon")[0]
        except IndexError:
            latitude  = nc.variables.get("lat")[:]
            longitude = nc.variables.get("lon")[:]

    o = Nco()
    o.ncks(
        input=netcdf_file,
        output=netcdf_file,
        options=[
            '-O',
            '-x',
            '-v', 'lat,lon'
        ]
    )

    with nc4.Dataset(netcdf_file, 'a') as nc:
        lat = nc.createVariable('lat', 'f4')
        lat[:] = latitude

        lon = nc.createVariable('lon', 'f4')
        lon[:] = longitude


def normalize_ctd_depths(netcdf_file):
    with nc4.Dataset(netcdf_file) as nc:

        dvars = []
        #dvars += nc.get_variables_by_attributes(generic_name='depth')
        #dvars += nc.get_variables_by_attributes(long_name=lambda x: x and x.strip().lower() == 'depth (m)')
        dvars += nc.get_variables_by_attributes(name='depth')

        dvars = list(set(dvars))

        if len(dvars) > 1:
            raise ValueError("Uhhhh, more than 1 depth variable? {}".format(dvars))

        depths = dvars[0][:]
        dname = dvars[0].name

    o = Nco()
    o.ncks(
        input=netcdf_file,
        output=netcdf_file,
        options=[
            '-O',
            '-x',
            '-v', dname
        ]
    )

    with nc4.Dataset(netcdf_file, 'a') as nc:
        nc.createDimension('z', depths.size)
        d = nc.createVariable('z', 'f4', ('z',))
        d[:] = depths


def normalize_adcp_location(netcdf_file):
    pass


def normalize_adcp_depths(netcdf_file):
    pass
