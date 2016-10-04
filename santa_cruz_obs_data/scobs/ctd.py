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
    normalize_locations,
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
    normalize_locations(file)
    normalize_ctd_depths(file)
    normalize_ctd_datavariables(file)
    normalize_ctd(file)

    # Now check for compliance!
    check_compliance(file)

    # Load as a CF DSG
    CFDataset.load(file)

    yield file


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
        z = nc.createVariable('depth', 'f4')
        z[:] = depths
        z_atts.update({
            'standard_name': 'depth',
            'axis': 'Z',
            'positive': 'down'
        })
        z.setncatts(z_atts)


def normalize_ctd_datavariables(netcdf_file):
    with CFDataset(netcdf_file, 'a') as nc:
        for v in nc.get_variables_by_attributes(axis=lambda v: v is None):
            v.coordinates = 'time depth lat lon'


def normalize_ctd(netcdf_file):
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
