#!python
# coding=utf-8

import numpy as np
from nco import Nco
import netCDF4 as nc4

from pyaxiom.netcdf.cf import CFDataset
from pyaxiom.netcdf.sensors.dsg import *

from scobs.utils import (
    normalize_time,
    normalize_epic_codes,
    normalize_vectors,
    normalize_units,
    normalize_netcdf4,
    normalize_locations,
    normalize_depth_locations,
    check_compliance
)
from scobs.mappings import global_attributes

import logging
logger = logging.getLogger(__name__)


def enhance_wvs(file):

    # Normalize
    normalize_netcdf4(file)
    normalize_time(file)
    normalize_epic_codes(file)
    normalize_vectors(file)
    normalize_units(file)
    normalize_locations(file)
    normalize_wvs_depths(file)
    normalize_wvs_datavariables(file)
    normalize_wvs(file)

    # Now check for compliance!
    check_compliance(file)

    # Load as a CF DSG
    CFDataset.load(file)

    yield file


def normalize_wvs_depths(netcdf_file):
    """
    Adjust some data variables to account for the height of the sensor above
    the sea floor
    """
    add_sensor_heights = [
        'sea_surface_height'
    ]
    with nc4.Dataset(netcdf_file, 'a') as nc:
        vs = nc.get_variables_by_attributes(
            standard_name=lambda x: x in add_sensor_heights,
            initial_sensor_height=lambda x: x is not None
        )
        for v in vs:
            v[:] = v[:] + v.initial_sensor_height
            v.NOTE = 'height relative to the sea floor'


def normalize_wvs_datavariables(netcdf_file):
    with nc4.Dataset(netcdf_file, 'a') as nc:
        for v in nc.get_variables_by_attributes(axis=lambda v: v is None):
            v.coordinates = 'time lat lon'


def normalize_wvs(netcdf_file):
    with CFDataset(netcdf_file, 'a') as nc:
        # Set the generic global attributes
        for n, v in global_attributes.items():
            nc.setncattr(n, v)

        nc.id = nc.MOORING
        nc.featureType = 'timeSeries'

        # Set the unique station ID to the MOORING id
        s = nc.createVariable('station', str)
        s[0] = nc.MOORING
        s_atts = {
            'cf_role': 'timeseries_id',
        }
        s.setncatts(s_atts)


def enhance_avg(file):
    # Normalize
    normalize_netcdf4(file)
    normalize_epic_codes(file)
    normalize_vectors(file)
    normalize_units(file)
    normalize_depth_locations(file)
    normalize_avg_datavariables(file)
    normalize_avg_depths(file)  # Changes dimension name from 'depth' to 'z'
    normalize_avg(file)
    normalize_time(file)

    # Now check for compliance!
    check_compliance(file)

    # Load as a CF DSG
    CFDataset.load(file)

    yield file


def normalize_avg_depths(netcdf_file):
    with CFDataset(netcdf_file, 'a') as nc:
        nc.renameDimension('depth', 'z')
        d = nc.variables['depth']
        d.axis = 'Z'

        # Add a special case variable for the non-depth measurements
        # that occur at the ADCP transducer.
        depth = nc.nominal_sensor_depth
        sdepth = nc.createVariable('station_depth', depth.dtype)
        sdepth.units = 'm'
        sdepth.standard_name = 'surface_altitude'
        sdepth.positive = 'down'
        sdepth.long_name = 'sensor depth below datum'
        sdepth[:] = depth


def normalize_avg_datavariables(netcdf_file):

    remove_attrs = [
        'minimum',
        'maximum',
        'valid_range',
    ]

    with nc4.Dataset(netcdf_file, 'a') as nc:
        for v in nc.get_variables_by_attributes(axis=lambda v: v is None):
            if 'depth' in v.dimensions:
                v.coordinates = 'time depth lat lon'
            else:
                v.coordinates = 'time lat lon'

            for ra in remove_attrs:
                if hasattr(v, ra):
                    v.delncattr(ra)


def normalize_avg(netcdf_file):
    with CFDataset(netcdf_file, 'a') as nc:
        # Set the generic global attributes
        for n, v in global_attributes.items():
            nc.setncattr(n, v)

        nc.id = nc.MOORING
        nc.featureType = 'timeSeriesProfile'

        # Set the unique station ID to the MOORING id
        s = nc.createVariable('station', str)
        s[0] = nc.MOORING
        s_atts = {
            'cf_role': 'timeseries_id',
        }
        s.setncatts(s_atts)

        # Create profile data
        p = nc.createVariable('profile', 'i4', ('time',))
        p[:] = np.arange(p.size)
        p_atts = {
            'cf_role': 'profile_id',
        }
        p.setncatts(p_atts)
