#!python
# coding=utf-8

import os
import sys

import unittest

import numpy as np
import netCDF4 as nc4

import logging
from pyaxiom import logger
logger.level = logging.INFO

# Hack to get imports working with pytests
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from collect import main


class ConverterTests(unittest.TestCase):

    def setUp(self):
        self.output = os.path.join(os.path.dirname(__file__), 'testing_output')
        self.download = os.path.join(os.path.dirname(__file__), 'testing_download')
        self.csv = os.path.join(os.path.dirname(__file__), '..', 'project_metadata.csv')

    def download_and_process(self, project, ncfile):
        downloaded_file = os.path.join(self.download, ncfile)
        output_file = os.path.join(self.output, project, ncfile)
        do_download = not os.path.isfile(downloaded_file)
        main(self.output, self.download, do_download, [project.lower()], self.csv, filesubset=[ncfile.lower()])
        assert os.path.isfile(downloaded_file)
        assert os.path.isfile(output_file)
        return output_file

    def test_timeseries(self):
        project = 'BUZZ_BAY'
        ncfile = '2881-A.cdf'
        output_file = self.download_and_process(project, ncfile)

        with nc4.Dataset(output_file) as nc:
            assert 'u_1205' in nc.variables        # ADCP variable
            assert 'latitude' in nc.variables
            assert 'longitude' in nc.variables
            assert 'z' in nc.variables

            assert np.isclose(nc.variables['longitude'][:], -70.88355)
            assert np.isclose(nc.variables['latitude'][:], 41.558887)
            assert nc.variables['z'][:] == 9.0

            assert 'z' not in nc.variables['u_1205'].dimensions

    def test_timeseries_profile(self):
        project = 'EUROSTRATAFORM'
        ncfile = '7031adc-a.nc'
        output_file = self.download_and_process(project, ncfile)

        with nc4.Dataset(output_file) as nc:
            assert 'u_1205' in nc.variables
            assert 'latitude' in nc.variables
            assert 'longitude' in nc.variables
            assert 'z' in nc.variables

            assert np.isclose(nc.variables['longitude'][:], 13.8795)
            assert np.isclose(nc.variables['latitude'][:], 43.3331)
            assert np.allclose(nc.variables['z'][:],
                               np.asarray([1.4548931, 2.204893, 2.954893, 3.704894, 4.454894, 5.204894,
                                           5.954894, 6.704894, 7.454894, 8.204894, 8.954894, 9.704894,
                                           10.454894, 11.204893, 11.954893, 12.704893, 13.454893, 14.204893]))
            assert 'z' in nc.variables['u_1205'].dimensions
            assert np.allclose(nc.variables['u_1205'][0, :],
                               np.asarray([0.80603516, 1.24185622, 1.14698148, 1.21925867, 1.23519087, 1.13773572,
                                           1.2186166, 0.98565829, 1.0108161, 1.1536535, 1.08373952, 1.10589528,
                                           1.13214898, 1.12784564, 1.1278441, 1.13214445, 1.17817998, 0.98327768]))

    def test_timeseries_profile_with_detached_parameter(self):
        project = 'EUROSTRATAFORM'
        ncfile = '7031adc-a.nc'
        output_file = self.download_and_process(project, ncfile)

        with nc4.Dataset(output_file) as nc:
            assert 'u_1205' in nc.variables        # ADCP variable
            assert 'latitude' in nc.variables
            assert 'longitude' in nc.variables
            assert 'z' in nc.variables
            assert 'Tx_1211' in nc.variables       # Bottom temperature
            assert 'sensor_depth' in nc.variables  # Special NCJ variable

            assert np.isclose(nc.variables['longitude'][:], 13.8795)
            assert np.isclose(nc.variables['latitude'][:], 43.3331)

            assert 'z' not in nc.variables['Tx_1211'].dimensions
            assert np.isclose(nc.variables['Tx_1211'][0], 16.81)

    def test_combine_multiple_depth_dimensions(self):
        project = 'MOBILE_BAY'
        ncfile = '3951tct-a.cdf'
        output_file = self.download_and_process(project, ncfile)

        with nc4.Dataset(output_file) as nc:
            assert 'temp' in nc.variables
            assert 'cond' in nc.variables
            assert 'salinity' in nc.variables
            assert 'trans' in nc.variables
            assert 'att' in nc.variables
            assert 'latitude' in nc.variables
            assert 'longitude' in nc.variables
            assert 'z' in nc.variables
            assert nc.variables['z'].size == 5
            assert np.allclose(nc.variables['z'][:],
                               np.asarray([1.0, 1.1, 3.0, 3.3, 3.4]))

            assert np.isclose(nc.variables['longitude'][:], -88.00183)
            assert np.isclose(nc.variables['latitude'][:], 30.466167)

            assert 'z' in nc.variables['temp'].dimensions
            assert 'z' in nc.variables['cond'].dimensions
            assert 'z' in nc.variables['salinity'].dimensions
            assert 'z' in nc.variables['trans'].dimensions
            assert 'z' in nc.variables['att'].dimensions

    def test_split_multiple_z_dimensions(self):
        project = 'FI12'
        ncfile = '9205advs-cal.nc'
        self.download_and_process(project, ncfile)
        # This should produce extra files for some variables, namely NEP2_56 and Sed2_981
        nep_file = os.path.join(self.output, project, ncfile.replace('.nc', '_NEP2_56.nc'))
        sed_file = os.path.join(self.output, project, ncfile.replace('.nc', '_Sed2_981.nc'))

        assert os.path.isfile(nep_file)
        assert os.path.isfile(sed_file)

        with nc4.Dataset(nep_file) as nc:
            assert 'NEP2_56' in nc.variables
            assert 'z' in nc.variables
            assert 'z' not in nc.variables['NEP2_56'].dimensions

        with nc4.Dataset(sed_file) as nc:
            assert 'Sed2_981' in nc.variables
            assert 'z' in nc.variables
            assert 'z' not in nc.variables['Sed2_981'].dimensions
