#!python
# coding=utf-8

import os
import sys
import shutil

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

    def tearDown(self):
        return
        try:
            shutil.rmtree(self.output)
        except OSError:
            pass

    def download_and_process(self, project, ncfile):
        downloaded_file = os.path.join(self.download, project, ncfile)
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

    def test_metadata_variables(self):
        project = 'FI14'
        ncfile = '10021wh-a.nc'
        output_file = self.download_and_process(project, ncfile)

        with nc4.Dataset(output_file) as nc:
            assert 'bindist' in nc.variables
            assert 'z' in nc.variables['bindist'].dimensions

    def test_split_multiple_z_dimensions(self):
        project = 'FI12'
        ncfile = '9205advs-cal.nc'
        output_file = self.download_and_process(project, ncfile)
        # This should produce extra files for some variables, namely NEP2_56 and Sed2_981

        z1_file = os.path.join(self.output, project, ncfile.replace('.nc', '_z1.nc'))
        z2_file = os.path.join(self.output, project, ncfile.replace('.nc', '_z2.nc'))
        z3_file = os.path.join(self.output, project, ncfile.replace('.nc', '_z3.nc'))

        assert os.path.isfile(z1_file)
        assert os.path.isfile(z2_file)
        assert os.path.isfile(z3_file)

        with nc4.Dataset(z1_file) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 920
            assert nc.id == os.path.splitext(os.path.basename(z1_file))[0]
            assert 'ATTN1_55' in nc.variables
            assert 'tran1_4010' in nc.variables
            assert 'z' in nc.variables
            # Make sure it was converted to positive "up" (from "down")
            assert nc.variables['z'].positive == 'up'
            assert np.isclose(nc.variables['z'][:], -17.7310000333786)
            assert 'z' not in nc.variables['ATTN1_55'].dimensions

        with nc4.Dataset(z2_file) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 920
            assert nc.id == os.path.splitext(os.path.basename(z2_file))[0]
            assert 'NEP2_56' in nc.variables
            assert 'Sed2_981' in nc.variables
            # Make sure it was converted to positive "up" (from "down")
            assert nc.variables['z'].positive == 'up'
            assert np.isclose(nc.variables['z'][:], -18.5710000071526)
            assert 'z' in nc.variables
            assert 'z' not in nc.variables['NEP2_56'].dimensions

        with nc4.Dataset(z3_file) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 920
            assert nc.id == os.path.splitext(os.path.basename(z3_file))[0]
            assert 'P_4023' in nc.variables
            assert 'SDP_850' in nc.variables
            # Make sure it was converted to positive "up" (from "down")
            assert nc.variables['z'].positive == 'up'
            assert np.isclose(nc.variables['z'][:], -17.2110000524521)
            assert 'z' in nc.variables
            assert 'z' not in nc.variables['P_4023'].dimensions

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

    def test_split_multiple_z_psc_cal_file(self):
        project = 'MVCO_11'
        ncfile = '9104pcs-cal.cdf'
        output_file = self.download_and_process(project, ncfile)

        # This should produce extra files for some variables
        z1_file = os.path.join(self.output, project, ncfile.replace('.cdf', '_z1.cdf'))
        z2_file = os.path.join(self.output, project, ncfile.replace('.cdf', '_z2.cdf'))

        assert os.path.isfile(z1_file)
        assert os.path.isfile(z2_file)

        with nc4.Dataset(z1_file) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 910
            assert nc.id == os.path.splitext(os.path.basename(z1_file))[0]
            assert 'Tx_1211' in nc.variables
            assert 'z' in nc.variables
            # Make sure it was converted to positive "up" (from "down")
            assert nc.variables['z'].positive == 'up'
            assert np.isclose(nc.variables['z'][:], -10.97)
            assert 'z' not in nc.variables['Tx_1211'].dimensions

        with nc4.Dataset(z2_file) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 910
            assert nc.id == os.path.splitext(os.path.basename(z2_file))[0]
            assert 'NEP2_56' in nc.variables
            # Make sure it was converted to positive "up" (from "down")
            assert nc.variables['z'].positive == 'up'
            assert np.isclose(nc.variables['z'][:], -10.8159999847412)
            assert 'z' in nc.variables
            assert 'z' not in nc.variables['NEP2_56'].dimensions

        with nc4.Dataset(output_file) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 910
            assert nc.id == os.path.splitext(ncfile)[0]
            assert 'u_1205' in nc.variables
            assert 'ATTN1_55' in nc.variables
            assert 'tran1_4010' in nc.variables
            # Make sure it was converted to positive "up" (from "down")
            assert nc.variables['z'].positive == 'up'
            assert np.allclose(nc.variables['z'][:],
                               np.asarray([-12.65049, -12.58649, -12.52249, -12.45849, -12.39449, -12.33049,
                                           -12.26649, -12.20249, -12.13849, -12.07449, -12.01049, -11.94649,
                                           -11.88249, -11.81849, -11.75449, -11.69049, -11.62649, -11.56249,
                                           -11.49849, -11.43449, -11.37049, -11.30649]))
            assert np.allclose(nc.variables['sensor_depth'][:],
                               -9.25)
            assert 'z' in nc.variables

    def test_variable_with_extra_dimensions(self):
        project = 'FI12'
        ncfile = '9261awWvs-cal.nc'
        output_file = self.download_and_process(project, ncfile)

        with nc4.Dataset(output_file) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 926
            assert nc.id == os.path.splitext(ncfile)[0]
            assert 'direction' in nc.variables
            assert 'direction' in nc.dimensions

            assert 'frequency' in nc.variables
            assert 'frequency' in nc.dimensions
            assert 'burst' in nc.dimensions
            assert 'SpecAmp' in nc.variables
            assert 'burst' in nc.variables['SpecAmp'].dimensions
            assert 'frequency' in nc.variables['SpecAmp'].dimensions

    def test_variable_with_extra_dimensions_size_1(self):
        project = 'WFAL'
        ncfile = '9411aqd-cal.nc'
        output_file = self.download_and_process(project, ncfile)

        with nc4.Dataset(output_file) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 941
            assert nc.id == os.path.splitext(ncfile)[0]
            assert 'bindist' in nc.variables
            assert 'z' not in nc.variables['bindist'].dimensions
            assert 'z' not in nc.dimensions

    def test_inconsistent_dimensions(self):
        project = 'PV_SHELF07'
        ncfile = '8482ls-cal.nc'
        output_file = self.download_and_process(project, ncfile)

        with nc4.Dataset(output_file) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 848
            assert nc.id == os.path.splitext(ncfile)[0]
            assert 'CDF' in nc.variables
            assert 'time' in nc.variables['CDF'].dimensions
            assert 'CDF' in nc.variables['CDF'].dimensions
            assert 'z' not in nc.variables['CDF'].dimensions

    def test_vector_conversion(self):
        project = 'MBAY_LT'
        ncfile = '4801spd-a_d1.nc'
        output_file = self.download_and_process(project, ncfile)

        with nc4.Dataset(output_file) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 480
            assert nc.id == os.path.splitext(ncfile)[0]
            assert nc.variables['platform'].type == 'fixed'
            assert 'CS_300' in nc.variables
            assert 'CD_310' in nc.variables

    def test_fillvalues(self):
        project = 'BW2011'
        ncfile = '9041ysi-a.nc'
        output_file = self.download_and_process(project, ncfile)

        with nc4.Dataset(output_file) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 904
            assert nc.id == os.path.splitext(ncfile)[0]

    def test_vertical_direction(self):
        project = 'MVCO_11'
        ncfile = '9111aqd-a.nc'
        output_file = self.download_and_process(project, ncfile)

        with nc4.Dataset(output_file) as nc:
            assert nc.original_folder == project
            assert nc.original_filename == ncfile
            assert nc.MOORING == 911
            assert nc.id == os.path.splitext(ncfile)[0]

            # Data under water should be negative (positive "up")
            for v in nc.variables['z'][:]:
                assert v < 0
