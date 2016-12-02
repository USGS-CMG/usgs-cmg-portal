#!python
# coding=utf-8

import os
import sys
import shutil
from glob import glob

import unittest

import numpy as np
import netCDF4 as nc4

import logging
from pyaxiom import logger
logger.level = logging.INFO

# Hack to get imports working with pytests
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from devices import main


class DeviceTests(unittest.TestCase):

    def setUp(self):
        self.output = os.path.join(os.path.dirname(__file__), 'device_output')
        self.download = os.path.join(os.path.dirname(__file__), 'cf_output')

    def download_and_process(self, project, ncfiles, feature, variables, since=None):
        lp = project.lower()
        output_files = []

        for n in ncfiles:
            downloaded_file = os.path.join(self.download, project, n)
            do_download = not os.path.isfile(downloaded_file)

        main(self.output, self.download, [project], do_download, filesubset=[ x.lower() for x in ncfiles ], since=since)

        for n in ncfiles:
            # CF file that was downloaded
            downloaded_file = os.path.join(self.download, project, n)
            assert os.path.isfile(downloaded_file)
            for v in variables:
                # Files that were created
                of = list(glob('{}/{}_{}/{}_{}_{}-*.nc'.format(self.output, lp, feature, lp, feature, v)))
                output_files += of

        return list(set(output_files))

    def test_combine_multiple_depth_files(self):
        project = 'MBAY_LT'
        ncfiles = ['4501spd-a_d1.nc', '4501spd-a_d2.nc', '4501spd-a_d3.nc', '4501spd-a_d4.nc']
        variables = ['u_1205']
        feature = 450
        output_files = self.download_and_process(project, ncfiles, feature, variables)

        assert len(output_files) == 1

        with nc4.Dataset(output_files[0]) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == '4501spd-a_d3.nc'  # The first file with 'u' in it
            assert nc.MOORING == feature
            assert nc.id == '{}_{}'.format(project.lower(), feature)
            assert 'eastward_sea_water_velocity' in nc.variables
            assert 'crs' in nc.variables
            assert 'latitude' in nc.variables
            assert 'longitude' in nc.variables
            assert 'height' in nc.variables
            assert nc.variables['height'].size == 2
            # Make sure it was sorted and converted to positive "down" (from "up")
            assert nc.variables['height'].positive == 'down'
            assert np.allclose(nc.variables['height'][:], np.asarray([31.3, 31.9], dtype=np.float32))

            assert np.isclose(nc.variables['longitude'][:], -70.78083)
            assert np.isclose(nc.variables['latitude'][:], 42.3745)

            assert 'time' in nc.variables['eastward_sea_water_velocity'].dimensions
            assert 'z' in nc.variables['eastward_sea_water_velocity'].dimensions

    def test_sensor_depth_attribute(self):
        project = 'EUROSTRATAFORM'
        ncfiles = ['7031adc-a.nc']
        variables = ['Tx_1211']
        feature = 703
        output_files = self.download_and_process(project, ncfiles, feature, variables)

        assert len(output_files) == 1

        with nc4.Dataset(output_files[0]) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == '7031adc-a.nc'  # The first file with 'u' in it
            assert nc.MOORING == feature
            assert nc.id == '{}_{}'.format(project.lower(), feature)
            assert 'sea_water_temperature' in nc.variables
            assert 'crs' in nc.variables
            assert 'latitude' in nc.variables
            assert 'longitude' in nc.variables
            assert 'height' in nc.variables
            assert nc.variables['height'].size == 1
            # Make sure it was sorted and converted to positive "down" (from "up")
            assert nc.variables['height'].positive == 'down'
            assert np.isclose(nc.variables['height'][0], 15.964893579483032)

            assert 'time' in nc.variables['sea_water_temperature'].dimensions
            assert 'z' not in nc.variables['sea_water_temperature'].dimensions
