#!python
# coding=utf-8

import os
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
    normalize_depth_locations,
    check_compliance
)
from scobs.mappings import global_attributes

import logging
logger = logging.getLogger(__name__)


def enhance_summary(file):
    # Normalize
    normalize_netcdf4(file)
    normalize_time(file)
    normalize_epic_codes(file)
    normalize_vectors(file)
    normalize_units(file)
    normalize_depth_locations(file)
    normalize_summary_depths(file)
    normalize_summary_datavariables(file)

    for nf in split_by_depths(file):
        normalize_summary(nf)

        # Now check for compliance!
        check_compliance(nf)

        # Load as a CF DSG
        CFDataset.load(nf)

        yield nf


def split_by_depths(netcdf_file):

    depthmap = {
        None: []
    }

    with CFDataset(netcdf_file) as nc:
        # Go through the depths of each variable and pull them out into
        # separate files
        for var in nc.get_variables_by_attributes(axis=lambda v: v is None):
            if hasattr(var, 'sensor_depth'):
                if var.sensor_depth not in depthmap:
                    depthmap[var.sensor_depth] = []
                depthmap[var.sensor_depth].append(var.name)
            else:
                depthmap[None].append(var.name)

    o = Nco()
    for i, (depth, vnames) in enumerate(depthmap.items()):

        if depth is None:
            continue

        base, ext = os.path.splitext(netcdf_file)
        outfile = '{}_z{}{}'.format(base, i, ext)

        allvars = vnames + ['lat', 'lon', 'depth', 'burst']

        o.ncks(
            input=netcdf_file,
            output=outfile,
            options=[
                '-O',
                '-h',
                '-v', ','.join(allvars)
            ]
        )

        # Now we open it up and add the correct value for the depth variable
        with CFDataset(outfile, 'a') as znc:
            z = znc.variables['depth']
            z[:] = depth

        yield outfile


def normalize_summary_depths(netcdf_file):

    redimension = []

    with CFDataset(netcdf_file, 'a') as nc:
        nc.renameDimension('depth', 'z')
        z_atts = nc.vatts('depth')

        # Get list of variables to re-write with different dimensions
        for vname, ov in nc.variables.items():
            if 'depth' in ov.dimensions:
                redimension.append(vname)

        for vname in redimension:
            ov = nc.variables[vname]
            vatts = nc.vatts(vname)
            nc.renameVariable(vname, '{}_old'.format(vname))
            nc.sync()

            v = nc.createVariable(vname, ov.dtype, ('time',))
            v.setncatts(vatts)
            v[:] = ov[:, 0]

    # Remove the old variables
    remove_vars = [ '{}_old'.format(vname) for vname in redimension ]
    remove_vars += ['depth']

    o = Nco()
    o.ncks(
        input=netcdf_file,
        output=netcdf_file,
        options=[
            '-O',
            '-h',
            '-x',
            '-v', ','.join(remove_vars)
        ]
    )

    # Add the new Z variable
    with nc4.Dataset(netcdf_file, 'a') as nc:
        z = nc.createVariable('depth', 'f4')
        z[:] = None
        z_atts.update({
            'standard_name': 'depth',
            'axis': 'Z',
            'positive': 'down'
        })
        z.setncatts(z_atts)


def normalize_summary_datavariables(netcdf_file):
    remove_attrs = [
        'minimum',
        'maximum',
        'valid_range',
    ]

    with nc4.Dataset(netcdf_file, 'a') as nc:
        for v in nc.get_variables_by_attributes(axis=lambda v: v is None):
            v.coordinates = 'time lat lon'

            for ra in remove_attrs:
                if hasattr(v, ra):
                    v.delncattr(ra)


def normalize_summary(netcdf_file):
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
