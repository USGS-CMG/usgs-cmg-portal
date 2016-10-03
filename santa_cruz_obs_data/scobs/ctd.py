#!python
# coding=utf-8

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
    check_compliance
)
from scobs.mappings import global_attributes

import logging
logger = logging.getLogger(__name__)


def enhance_ctd(file):
    # Normalize
    normalize_netcdf4(file)
    normalize_time(file)
    normalize_epic_codes(file)
    normalize_vectors(file)
    normalize_units(file)
    normalize_ctd_location(file)
    normalize_ctd_depths(file)
    normalize_ctd_datavariables(file)
    normalize(file)

    # Now check for compliance!
    scores = check_compliance(file)
    for s in scores:
        level = 'info'
        if not s.passed:
            level = 'error'
        getattr(logger, level)(s)

    # Load as a CF DSG
    CFDataset.load(file)


def normalize_ctd_location(netcdf_file):
    with CFDataset(netcdf_file) as nc:
        y = nc.variables.get("lat")
        y_data = y[:][0]
        y_atts = nc.vatts('lat')

        x = nc.variables.get("lon")
        x_data = x[:][0]
        x_atts = nc.vatts('lon')

    # Remove the old x/y variable so we can add new ones with no dimensions
    o = Nco()
    o.ncks(
        input=netcdf_file,
        output=netcdf_file,
        options=[
            '-O',
            '-h',
            '-x',
            '-v', 'lat,lon'
        ]
    )

    # Add the new X/Y variables
    with nc4.Dataset(netcdf_file, 'a') as nc:
        lat = nc.createVariable('lat', 'f4')
        lat[:] = y_data
        y_atts.update({
            'standard_name': 'latitude',
            'axis': 'Y'
        })
        lat.setncatts(y_atts)

        lon = nc.createVariable('lon', 'f4')
        lon[:] = x_data
        x_atts.update({
            'standard_name': 'longitude',
            'axis': 'X'
        })
        lon.setncatts(x_atts)


def normalize_ctd_depths(netcdf_file):
    with CFDataset(netcdf_file, 'a') as nc:
        depths = nc.variables['depth'][:][0]
        z_atts = nc.vatts('depth')

    # Remove the old depth variable so we can add a new one with no dimensions
    o = Nco()
    o.ncks(
        input=netcdf_file,
        output=netcdf_file,
        options=[
            '-O',
            '-h',
            '-x',
            '-v', 'depth'
        ]
    )

    # Add the new Z variable
    with nc4.Dataset(netcdf_file, 'a') as nc:
        z = nc.createVariable('z', 'f4')
        z[:] = depths
        z_atts.update({
            'standard_name': 'depth',
            'axis': 'Z',
            'positive': 'down'
        })
        z.setncatts(z_atts)


def normalize_ctd_datavariables(netcdf_file):
    with CFDataset(netcdf_file, 'a') as nc:

        nc.featureType = 'timeSeries'

        for v in nc.get_variables_by_attributes(axis=lambda v: v is None):
            v.coordinates = 'time lat lon z'


def normalize(netcdf_file):
    with CFDataset(netcdf_file, 'a') as nc:
        # Set the generic globall attributes
        for n, v in global_attributes.items():
            nc.setncattr(n, v)

        nc.id = nc.MOORING

        # Set the unique station ID to the MOORING id
        s = nc.createVariable('station', str)
        s[0] = nc.MOORING
        s_atts = {
            'cf_role': 'timeseries_id',
        }
        s.setncatts(s_atts)
