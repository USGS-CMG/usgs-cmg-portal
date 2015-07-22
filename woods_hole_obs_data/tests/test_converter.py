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
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 288
            assert nc.id == os.path.splitext(ncfile)[0]
            assert 'u_1205' in nc.variables        # ADCP variable
            assert 'latitude' in nc.variables
            assert 'longitude' in nc.variables
            assert 'z' in nc.variables

            assert np.isclose(nc.variables['longitude'][:], -70.88355)
            assert np.isclose(nc.variables['latitude'][:], 41.558887)
            # Make sure it was converted to positive "up" (from "down")
            assert np.isclose(nc.variables['z'][:], -9.0)

            assert 'z' not in nc.variables['u_1205'].dimensions

    def test_timeseries_4_digit_mooring(self):
        project = 'FI14'
        ncfile = '10002sc-a.nc'
        output_file = self.download_and_process(project, ncfile)

        with nc4.Dataset(output_file) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 1000
            assert nc.id == os.path.splitext(ncfile)[0]
            assert 'T_28' in nc.variables
            assert 'latitude' in nc.variables
            assert 'longitude' in nc.variables
            assert 'z' in nc.variables

            assert np.isclose(nc.variables['longitude'][:], -73.14838)
            assert np.isclose(nc.variables['latitude'][:], 40.636913)
            # Make sure it was converted to positive "up" (from "down")
            assert np.isclose(nc.variables['z'][:], -9.81)

            assert 'z' not in nc.variables['T_28'].dimensions

    def test_timeseries_profile(self):
        project = 'EUROSTRATAFORM'
        ncfile = '7031adc-a.nc'
        output_file = self.download_and_process(project, ncfile)

        with nc4.Dataset(output_file) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 703
            assert nc.id == os.path.splitext(ncfile)[0]
            assert 'u_1205' in nc.variables
            assert 'latitude' in nc.variables
            assert 'longitude' in nc.variables
            assert 'z' in nc.variables
            assert nc.variables['z'].positive == 'up'

            assert np.isclose(nc.variables['longitude'][:], 13.8795)
            assert np.isclose(nc.variables['latitude'][:], 43.3331)
            # Make sure it was sorted and converted to positive "up" (from "down")
            assert np.allclose(nc.variables['z'][:],
                               np.asarray([14.204893, 13.454893, 12.704893, 11.954893, 11.204893, 10.454894,
                                           9.704894,  8.954894,  8.204894,  7.454894,  6.704894,  5.954894,
                                           5.204894,  4.454894,  3.704894,  2.954893,  2.204893,  1.4548931]) * -1)
            assert 'z' in nc.variables['u_1205'].dimensions
            assert np.allclose(nc.variables['u_1205'][0, :],
                               np.asarray([0.09832777, 0.11781799, 0.113214445, 0.1127844, 0.11278457, 0.11321489,
                                           0.11058952, 0.10837395, 0.11536535, 0.10108161, 0.09856583, 0.12186166,
                                           0.113773575, 0.12351909, 0.12192587, 0.11469815, 0.12418562, 0.08060351]))

    def test_timeseries_profile_with_detached_parameter(self):
        project = 'EUROSTRATAFORM'
        ncfile = '7031adc-a.nc'
        output_file = self.download_and_process(project, ncfile)

        with nc4.Dataset(output_file) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 703
            assert nc.id == os.path.splitext(ncfile)[0]
            assert 'u_1205' in nc.variables        # ADCP variable
            assert 'latitude' in nc.variables
            assert 'longitude' in nc.variables
            assert 'z' in nc.variables
            assert nc.variables['z'].positive == 'up'
            assert 'Tx_1211' in nc.variables       # Bottom temperature
            assert 'sensor_depth' in nc.variables  # Special NCJ variable

            assert np.isclose(nc.variables['longitude'][:], 13.8795)
            assert np.isclose(nc.variables['latitude'][:], 43.3331)

            assert 'z' not in nc.variables['Tx_1211'].dimensions
            # Make sure it was converted to positive "up" (from "down")
            assert np.isclose(nc.variables['sensor_depth'][:], -15.964893579483032)
            assert np.isclose(nc.variables['Tx_1211'].sensor_depth, -15.964893579483032)

            assert np.isclose(nc.variables['Tx_1211'][0], 16.81)

    def test_combine_multiple_depth_dimensions(self):
        project = 'MOBILE_BAY'
        ncfile = '3951tct-a.cdf'
        output_file = self.download_and_process(project, ncfile)

        with nc4.Dataset(output_file) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 395
            assert nc.id == os.path.splitext(ncfile)[0]
            assert 'temp' in nc.variables
            assert 'cond' in nc.variables
            assert 'salinity' in nc.variables
            assert 'trans' in nc.variables
            assert 'att' in nc.variables
            assert 'latitude' in nc.variables
            assert 'longitude' in nc.variables
            assert 'z' in nc.variables
            assert nc.variables['z'].size == 5
            # Make sure it was sorted and converted to positive "up" (from "down")
            assert nc.variables['z'].positive == 'up'
            assert np.allclose(nc.variables['z'][:],
                               np.asarray([3.4, 3.3, 3.0, 1.1, 1.0]) * -1)

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
        output_file = self.download_and_process(project, ncfile)
        # This should produce extra files for some variables, namely NEP2_56 and Sed2_981
        nep_file = os.path.join(self.output, project, ncfile.replace('.nc', '_NEP2_56.nc'))
        sed_file = os.path.join(self.output, project, ncfile.replace('.nc', '_Sed2_981.nc'))

        assert os.path.isfile(nep_file)
        assert os.path.isfile(sed_file)

        with nc4.Dataset(nep_file) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 920
            assert nc.id == os.path.splitext(ncfile)[0]
            assert 'NEP2_56' in nc.variables
            assert 'z' in nc.variables
            # Make sure it was converted to positive "up" (from "down")
            assert nc.variables['z'].positive == 'up'
            assert np.isclose(nc.variables['z'][:], -18.571000007152556)
            assert 'z' not in nc.variables['NEP2_56'].dimensions

        with nc4.Dataset(sed_file) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 920
            assert nc.id == os.path.splitext(ncfile)[0]
            assert 'Sed2_981' in nc.variables
            # Make sure it was converted to positive "up" (from "down")
            assert nc.variables['z'].positive == 'up'
            assert np.isclose(nc.variables['z'][:], -18.571000007152556)
            assert 'z' in nc.variables
            assert 'z' not in nc.variables['Sed2_981'].dimensions

        with nc4.Dataset(output_file) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 920
            assert nc.id == os.path.splitext(ncfile)[0]
            assert 'w_1204min' in nc.variables
            # Make sure it was converted to positive "up" (from "down")
            assert nc.variables['z'].positive == 'up'
            assert np.isclose(nc.variables['z'][:], -19.840353)
            assert 'z' in nc.variables
            # Translated to m/s from cm/s
            assert np.isclose(nc.variables['w_1204min'][0], -0.06763188)
