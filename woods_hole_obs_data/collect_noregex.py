#!python
# coding=utf-8

import os
import csv
import shutil
import logging
import tempfile
import requests
import argparse
from copy import copy
from glob import glob
from datetime import datetime

import pytz
import epic2cf
from cf_units import Unit
from epic2cf.data import *
import netCDF4
import numpy as np

from thredds_crawler.crawl import Crawl
from pyaxiom.netcdf.sensors import TimeSeries
from pyaxiom.netcdf.dataset import EnhancedDataset
from pyaxiom.utils import DotDict

import coloredlogs

# Log to stdout
logger = logging.getLogger()
coloredlogs.install(level=logging.DEBUG)
logger.setLevel(logging.DEBUG)

# Log to a file
df = logging.FileHandler('cf_errors.log', mode='w', encoding='utf-8')
df.setLevel(logging.ERROR)
df.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(df)

# Don't show the HTTP connection spam
requests_log = logging.getLogger("requests").setLevel(logging.WARNING)
crawler_log = logging.getLogger("thredds_crawler").setLevel(logging.INFO)
epic_log = logging.getLogger("epic2cf").setLevel(logging.INFO)
pyaxiom_log = logging.getLogger("pyaxiom").setLevel(logging.ERROR)

IGNORABLE_CODES = location_codes + time_codes + generic_codes + voltage_codes


# Special case EPIC mapping for generic EPIC codes that are used
# EG '20' can be  Air Temperature or Water Temperature.
def correct_temperature(var, filename):
    # https://github.com/USGS-CMG/usgs-cmg-portal/issues/170#issuecomment-189485296
    for x in ['met', 'hlm', 'alm', 'hwlb']:
        if x in filename:
            return epic2cf.mapping.get(21)
    return epic2cf.mapping.get(28)


def correct_backscatter(var, filename):
    return DotDict(standard_name='backscatter_intensity', long_name='Backscatter Intensity', units='v', convert=lambda x: x, cf_units='v', cell_methods=None)


def correct_airpressure(var, filename):
    z = dict(
        standard_name='air_pressure',
        long_name='Air Pressure',
        cf_units='Pa',
        convert=lambda x: Unit(var.units.lower()).convert(x, Unit('Pa')),
        units=var.units.lower(),
        cell_methods=None
    )
    return DotDict(**z)


special_map = {
    20   : correct_temperature,
    56   : correct_backscatter,
    915  : correct_airpressure
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
    'fDOMQSU'    : dict(epic_code=None, overrides=dict(standard_name='concentration_of_colored_dissolved_organic_matter_in_sea_water_expressed_as_equivalent_mass_fraction_of_quinine_sulfate_dihydrate',
                                                       units='1')),
    'BPR_1301'  : dict(epic_code=1301),
    'WL_NAVD88' : dict(epic_code=18, overrides=dict(vertical_datum='NAVD88')),
    'CTDCON_4218': dict(epic_code=4218),
    'BP_915'     : dict(epic_code=915),
    'BPR_915'    : dict(epic_code=915),
    'Baro'       : dict(epic_code=915),
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
    'barometric pressure':                  dict(epic_code=915)
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

coord_vars      = ['feature_type_instance', 'time', 'time2', 'time_cf', 'old_time', 'depth', 'depth002', 'depth003', 'depth004', 'depth005', 'lat', 'lon']


def download(folder, project_metadata, filesubset, since):

    # Use thredds_crawler to find DAP endpoints of the RAW data.
    total_datasets = []
    skips = Crawl.SKIPS + ['.*OTHER.*', '.*ancillary.*', '.*OLD_VERSIONS.*']

    try:
        for k, v in project_metadata.items():
            # http://regexr.com/3conn
            datasets = Crawl(v['catalog_xml'],
                             select=['(.*.nc$)'],
                             skip=skips,
                             after=since).datasets
            logger.info("Found {0} datasets in {1}!".format(len(datasets), k))
            total_datasets += datasets
        logger.info("Found {0} TOTAL datasets!".format(len(total_datasets)))
    except KeyboardInterrupt:
        logger.info("Breaking out of crawling loop.")
        total_datasets = []

    try:
        os.makedirs(folder)
    except OSError:
        pass

    # Save datasets to download directory
    saved_files = []
    for num, d in enumerate(total_datasets):

        if filesubset and d.name.lower() not in filesubset:
            continue

        try:
            http_url = next(s["url"] for s in d.services if s["service"].lower() == "httpserver")
            project_name = http_url.split("/")[-2]
        except StopIteration:
            logger.error("No HTTPServer endpoint found, skipping")
            continue

        # Make download folder
        save_file = os.path.join(folder, project_name, d.name)
        if not os.path.isdir(os.path.dirname(save_file)):
            os.makedirs(os.path.dirname(save_file))
        logger.info("Downloading {0}".format(http_url))
        try:
            with open(save_file, "wb") as f:
                r = requests.get(http_url, stream=True)
                if not r.ok:
                    logger.error("Could not download '{!s}' from '{!s}', skipping".format(d.name, http_url))
                    break
                for block in r.iter_content(1024):
                    if not block:
                        break
                    f.write(block)
        except KeyboardInterrupt:
            logger.info("Breaking out of download loop.")
            raise
        except BaseException:
            logger.error("Could not download... error with HTTP endpoint.  Skipping.")
            continue

        # Try to open file, if it fails, writing failed.
        try:
            with EnhancedDataset(save_file, 'a') as nc:
                name, _ = os.path.splitext(d.name)
                nc.id = "{0}/{1}".format(project_name, name)
        except BaseException:
            os.remove(save_file)
            raise
        else:
            logger.info("{!s} saved ({!s}/{!s})".format(d.name, num + 1, len(total_datasets)))
            saved_files.append(save_file)

    return saved_files


def convert_attributes(ncvar, convert_function):

    attrs = [
        'valid_min',
        'valid_max',
        'valid_range',
        'actual_min',
        'actual_max',
        'actual_range',
        'minimum',
        'maximum',
    ]

    vatts = { k : getattr(ncvar, k) for k in ncvar.ncattrs() }
    for k, v in vatts.items():
        if k in attrs:
            attvalue = getattr(ncvar, k)
            vs = ncvar.dtype.type(attvalue)
            calc = convert_function(vs)
            setattr(ncvar, k, calc)


def normalize_epic_codes(netcdf_file, original_filename):
    with EnhancedDataset(netcdf_file, 'a') as nc:
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
                        convert_attributes(nc_var, d)
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
                                convert_attributes(nc_var, d)
                            elif k != 'original_units':
                                nc_var.setncattr(k, d)

            if hasattr(nc_var, "epic_code") and nc_var.epic_code:
                try:
                    epic_code = int(nc_var.epic_code)
                except ValueError:
                    logger.debug("No EPIC code specified on {0}".format(v))
                else:

                    # Specialized cases for generic EPIC codes
                    if epic_code in special_map:
                        attribs = special_map.get(epic_code)(nc_var, original_filename)
                    else:
                        attribs = epic2cf.mapping.get(epic_code)

                    # Special case for 'Onset weather stations'.
                    # https://github.com/USGS-CMG/usgs-cmg-portal/issues/69
                    if epic_code in [905, 908] and 'hml' in netcdf_file.lower():
                        attribs.standard_name = 'surface_downwelling_photosynthetic_radiative_flux_in_air'

                    if attribs is not None and attribs.standard_name is not None:
                        # Convert data to CF units
                        nc_var[:] = attribs.convert(nc_var[:])
                        convert_attributes(nc_var, attribs.convert)
                        # Set attributes
                        nc_var.standard_name = attribs.standard_name
                        nc_var.long_name     = attribs.long_name
                        nc_var.units         = attribs.cf_units
                        nc_var.epic_code     = epic_code  # Set it again to be sure it is an int
                        if attribs.cell_methods is not None:
                            nc_var.cell_methods = attribs.cell_methods
                    else:
                        logger.debug("Could not find CF mapping for EPIC code {!s}".format(epic_code))


def normalize_vectors(netcdf_file):
    with EnhancedDataset(netcdf_file, 'a') as nc:
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
            speed = np.sqrt(np.square(east[:]) + np.square(north[:]))
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


def normalize_units(netcdf_file):
    with EnhancedDataset(netcdf_file, 'a') as nc:
        for v in nc.variables:
            nc_var = nc.variables.get(v)
            if hasattr(nc_var, 'units') and nc_var.units == "K":
                # Convert kelvin to Celsius
                def d(x):
                    return x - 273.15
                nc_var[:] = d(nc_var[:])
                nc_var.units = "degree_Celsius"
                convert_attributes(nc_var, d)
            elif hasattr(nc_var, 'standard_name') and nc_var.standard_name == 'sea_surface_wave_from_direction':
                # Convert "From" to "To" direction
                def d(x):
                    return (x + 180) % 360
                nc_var[:] = d(nc_var[:])
                nc_var.standard_name = 'sea_surface_wave_to_direction'
                nc_var.long_name = "Wave Direction (to TN)"
                convert_attributes(nc_var, d)


def normalize_time(netcdf_file):
    epoch_units       = 'seconds since 1970-01-01T00:00:00Z'
    millisecond_units = 'milliseconds since 1858-11-17T00:00:00Z'

    with EnhancedDataset(netcdf_file, 'a') as nc:
        # Signell said this works, any problems and we can all blame him!
        time_data = netCDF4.num2date((np.int64(nc.variables['time'][:])-2400001)*3600*24*1000 + nc.variables['time2'][:].__array__(), units=millisecond_units)  # noqa
        nc.renameVariable("time", "old_time")
        nc.sync()

        time = nc.createVariable('time', 'f8', ('time'))
        time.units          = epoch_units
        time.standard_name  = "time"
        time.long_name      = "time of measurement"
        time.calendar       = "gregorian"
        time[:] = netCDF4.date2num(time_data, units=epoch_units).round()
        return time_data[0]


def normalize_variable_attribute_types(netcdf_file):

    attrs = [
        'valid_min',
        'valid_max',
        'valid_range',
        'actual_min',
        'actual_max',
        'actual_range',
        'minimum',
        'maximum',
    ]

    with EnhancedDataset(netcdf_file, 'a') as nc:
        for ncvarname, ncvar in nc.variables.items():
            vatts = { k : getattr(ncvar, k) for k in ncvar.ncattrs() }
            for k, v in vatts.items():
                if k in attrs:
                    attvalue = getattr(ncvar, k)
                    if isinstance(attvalue, np.ndarray):
                        if attvalue.dtype.kind in ['U', 'S']:
                            ncvar.delncattr(k)
                            continue
                    elif not attvalue or isinstance(attvalue, str):
                        ncvar.delncattr(k)
                        continue

                    vs = ncvar.dtype.type(attvalue)
                    setattr(ncvar, k, vs)


def main(output, download_folder, do_download, projects, csv_metadata_file, filesubset=None, since=None):
    project_metadata = dict()
    with open(csv_metadata_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            project_name = row['project_name']
            if isinstance(project_name, str) and project_name[0] == '#':
                continue
            if projects and project_name.lower() not in projects:
                # Skip projects if a subset was defined
                continue
            project_metadata[project_name] = dict()
            for k, v in row.items():
                project_metadata[project_name][k] = v

    if do_download:
        try:
            downloaded_files = download(download_folder, project_metadata, filesubset, since)
        except KeyboardInterrupt:
            logger.exception('Error downloading datasets from THREDDS')
            downloaded_files = []
    else:
        downloaded_files = glob(os.path.join(download_folder, '**', '*'))
        if since is not None:
            def should_keep(d):
                modt = datetime.utcfromtimestamp(os.path.getmtime(d)).replace(tzinfo=pytz.utc)
                return modt >= since
            downloaded_files = [ dl for dl in downloaded_files if should_keep(dl) ]

    for down_file in sorted(downloaded_files):

        temp_fd, temp_file = tempfile.mkstemp(prefix='cmg_collector', suffix='nc')
        try:

            if filesubset is not None:
                if os.path.basename(down_file).lower() not in filesubset:
                    # aka "9631ecp-a.nc"
                    # Skip this file!
                    continue

            project_name = os.path.basename(os.path.dirname(down_file))
            if projects:
                if project_name.lower() not in projects:
                    # Skip this project!
                    continue
            shutil.copy(down_file, temp_file)

            # Cleanup to CF-1.6
            try:
                first_time = normalize_time(temp_file)
            except (TypeError, ValueError, IndexError):
                logger.exception("Could not normalize the time variable. Skipping {0}.".format(down_file))
                continue
            except OverflowError:
                logger.error("Dates out of range. Skipping {0}.".format(down_file))
                continue

            normalize_variable_attribute_types(temp_file)
            normalize_epic_codes(temp_file, down_file)
            normalize_vectors(temp_file)
            normalize_units(temp_file)

            # Create list of variables that we want to save.
            mooring_id   = None
            latitude     = None
            longitude    = None

            fname = os.path.basename(down_file)
            feature_name, file_ext = os.path.splitext(os.path.basename(down_file))
            try:
                if int(fname[0]) <= 9 and int(fname[0]) >= 2:
                    # 1.) everything with first char between 2-9 is 3-digit
                    mooring_id = int(fname[0:3])
                elif int(fname[0]) == 1:
                    # 2.) if MOORING starts with 1, and data is newer than 2014, it's 4 digit, otherwise 3 digit.
                    if first_time > datetime(2014, 1, 1, 0):
                        # 4 digit if after Jan 1, 2014
                        mooring_id = int(fname[0:4])
                    else:
                        # 3 digit if before
                        mooring_id = int(fname[0:3])
            except ValueError:
                logger.exception("Could not create a suitable station_id. Skipping {0}.".format(down_file))
                continue

            file_name = os.path.basename(down_file)
            output_directory = os.path.join(output, project_name)
            logger.info("Translating {0} into CF1.6 format: {1}".format(down_file, os.path.abspath(os.path.join(output_directory, file_name))))

            with EnhancedDataset(temp_file) as nc:

                try:
                    latitude  = nc.variables.get("lat")[0]
                    longitude = nc.variables.get("lon")[0]
                except IndexError:
                    latitude  = nc.variables.get("lat")[:]
                    longitude = nc.variables.get("lon")[:]
                except TypeError:
                    logger.error("Could not find lat/lon variables. Skipping {0}.".format(down_file))

                file_global_attributes = { k : getattr(nc, k) for k in nc.ncattrs() }
                file_global_attributes.update(global_attributes)
                file_global_attributes['id'] = feature_name
                file_global_attributes['MOORING'] = mooring_id
                file_global_attributes['original_filename'] = fname
                file_global_attributes['original_folder'] = project_name

                no_override = ['id', 'MOORING', 'original_filename', 'original_folder', 'catalog_xml', 'project_name']
                if project_name in project_metadata:
                    for k, v in project_metadata[project_name].items():
                        if v and k.lower() not in no_override:
                            file_global_attributes[k] = v

                if 'summary' in file_global_attributes:
                    # Save the original summary
                    file_global_attributes['WHOI_Buoy_Group_summary'] = file_global_attributes['summary']

                # Better title/summary for discovery via catalogs
                project_title = file_global_attributes.get('project_title', project_name).strip()
                project_summary = file_global_attributes.get('project_summary', '').strip()
                file_global_attributes['title'] = 'USGS-CMG time-series data: {0} - {1} - {2}'.format(project_name, mooring_id, feature_name)
                file_global_attributes['summary'] = 'USGS-CMG time-series data from the {} project, mooring {} and package {}. {}'.format(project_title, mooring_id, feature_name, project_summary).strip()

                times  = nc.variables.get('time')[:]

                # Get all depth values
                depth_variables = []
                for dv in nc.variables:
                    depth_variables += [ x for x in nc.variables.get(dv).dimensions if 'depth' in x ]
                depth_variables = sorted(list(set(depth_variables)))

                try:
                    assert depth_variables
                    depth_values = np.asarray([ nc.variables.get(x)[:] for x in depth_variables ]).flatten()
                except (AssertionError, TypeError):
                    logger.warning("No depth variables found in {}, skipping.".format(down_file))
                    continue

                # Convert everything to positive up, unless it is specifically specified as "up" already
                depth_conversion = -1.0
                if depth_variables:
                    pull_positive = nc.variables.get(depth_variables[0])
                    if hasattr(pull_positive, 'positive') and pull_positive.positive.lower() == 'up':
                        depth_conversion = 1.0
                depth_values = depth_values * depth_conversion

                if not os.path.isdir(output_directory):
                    os.makedirs(output_directory)
                ts = TimeSeries(output_directory, latitude, longitude, feature_name, file_global_attributes, times=times, verticals=depth_values, output_filename=file_name, vertical_positive='up')

                # Set the platform type from the global attribute 'platform_type', defaulting to 'fixed'
                with EnhancedDataset(ts.out_file, 'a') as onc:
                    platform_type = getattr(onc, 'platform_type', 'fixed').lower()
                    onc.variables['platform'].setncattr('type', platform_type)
                    onc.variables['platform'].setncattr('nodc_name', "FIXED PLATFORM, MOORINGS")
                    # Add ERDDAP variables
                    onc.cdm_data_type = "TimeSeries"
                    onc.cdm_timeseries_variables = "latitude,longitude,z,feature_type_instance"

                v = []
                depth_files = []
                for other in sorted(nc.variables):  # Sorted for a reason... don't change!
                    try:
                        if other in coord_vars:
                            continue

                        ovsd = None  # old var sensor depth
                        old_var = nc.variables.get(other)
                        variable_attributes = { k : getattr(old_var, k) for k in old_var.ncattrs() }
                        # Remove/rename some attributes
                        # https://github.com/USGS-CMG/usgs-cmg-portal/issues/67
                        if 'valid_range' in variable_attributes:
                            del variable_attributes['valid_range']
                        if 'minimum' in variable_attributes:
                            variable_attributes['actual_min'] = variable_attributes['minimum']
                            del variable_attributes['minimum']
                        if 'maximum' in variable_attributes:
                            variable_attributes['actual_max'] = variable_attributes['maximum']
                            del variable_attributes['maximum']
                        if 'sensor_depth' in variable_attributes:
                            # sensor_depth is ALWAYS positive "down", so don't convert!
                            # This is contrary to the "positive" attribute on the Z axis.
                            # variable_attributes['sensor_depth'] = variable_attributes['sensor_depth'] * -1
                            # Round the sensor_depth attribute
                            variable_attributes['sensor_depth'] = np.around(variable_attributes['sensor_depth'], decimals=4)
                            ovsd = np.around(old_var.sensor_depth * depth_conversion, decimals=4)

                        fillvalue = None
                        if hasattr(old_var, "_FillValue"):
                            fillvalue = old_var._FillValue

                        # Figure out if this is a variable that is repeated at different depths
                        # as different variable names.   Assumes sorted.
                        new_var_name = other.split('_')[0]
                        if new_var_name in ts.ncd.variables:
                            # Already in new file (processed when the first was encountered in the loop below)
                            continue

                        # Get the depth index
                        depth_variable = [ x for x in old_var.dimensions if 'depth' in x ]
                        if depth_variable and len(old_var.dimensions) > 1 and 'time' in old_var.dimensions:
                            depth_index = np.squeeze(np.where(depth_values == (nc.variables.get(depth_variable[0])[:] * depth_conversion)))

                            # Find other variable names like this one
                            depth_indexes = [(other, depth_index)]
                            for search_var in sorted(nc.variables):
                                # If they have different depth dimension names we need to combine them into one variable
                                if search_var != other and search_var.split('_')[0] == new_var_name and \
                                   depth_variable[0] != [ x for x in nc.variables[search_var].dimensions if 'depth' in x ][0]:
                                    # Found a match at a different depth
                                    search_depth_variable = [ x for x in nc.variables.get(search_var).dimensions if 'depth' in x ]
                                    depth_index = np.squeeze(np.where(depth_values == (nc.variables.get(search_depth_variable[0])[:] * depth_conversion)))
                                    depth_indexes.append((search_var, depth_index))
                                    logger.info("Combining '{}' with '{}' as '{}' (different variables at different depths but are the same parameter)".format(search_var, other, new_var_name))

                            values = np.ma.empty((times.size, len(depth_values)), dtype=old_var.dtype)
                            values.fill_value = values.dtype.type(fillvalue)
                            fillvalue = values.dtype.type(fillvalue)
                            values.mask = True
                            inconsistent = False
                            for nm, index in depth_indexes:
                                try:
                                    values[:, index] = np.squeeze(nc.variables.get(nm)[:])
                                except ValueError:
                                    inconsistent = True
                                    break

                            # If we just have one index we want to use the original name
                            if len(depth_indexes) == 1:
                                # Just use the original variable name
                                new_var_name = other

                            if inconsistent is True:
                                # Incorrect array size, most likely a strange variable
                                ts.add_variable_object(old_var, dimension_map=dict(depth='z'), reduce_dims=True)
                            else:
                                # Create this one, should be the first we encounter for this type
                                ts.add_variable(new_var_name, values=values, times=times, fillvalue=fillvalue, attributes=variable_attributes)

                        elif len(old_var.dimensions) == 1 and old_var.dimensions[0] == 'time':
                            # A single time dimensioned variable, like pitch, roll, record count, etc.
                            ts.add_variable(other, values=old_var[:], times=times, unlink_from_profile=True, fillvalue=fillvalue, attributes=variable_attributes)
                        elif old_var.ndim <= 3 and ovsd and \
                                ((depth_values.size == 1 and not depth_variable and 'time' in old_var.dimensions) or
                                 (depth_values.size  > 1 and not depth_variable and 'time' in old_var.dimensions and 'sensor_depth' in ts.ncd.variables)):

                            if 'sensor_depth' in ts.ncd.variables and np.isclose(ts.ncd.variables['sensor_depth'][:], ovsd):
                                ts.add_variable(other, values=old_var[:], times=times, unlink_from_profile=True, verticals=[ovsd], fillvalue=fillvalue, attributes=variable_attributes)
                            else:
                                # Search through secondary files that have been created for detached variables at a certain depth and
                                # try to match this variable with one of the depths.
                                found_df = False
                                for dfts in depth_files:
                                    if isinstance(ovsd, np.ndarray):
                                        # Well, this is a bad file.
                                        raise ValueError("The sensor_depth attribute has more than one value, please fix the source NetCDF: {}".format(down_file))
                                    if np.isclose(dfts.ncd.variables[ts.vertical_axis_name][:], ovsd):
                                        dfts.add_variable(other, values=old_var[:], times=times, unlink_from_profile=True, verticals=[ovsd], fillvalue=fillvalue, attributes=variable_attributes)
                                        found_df = True
                                        break

                                # If we couldn't match the current or one of the existing secondary depth files, create a new one.
                                if found_df is False:
                                    new_file_name = file_name.replace(file_ext, '_z{}{}'.format(len(depth_files) + 1, file_ext))
                                    fga = copy(file_global_attributes)
                                    fga['id'] = os.path.splitext(new_file_name)[0]
                                    new_ts = TimeSeries(output_directory, latitude, longitude, feature_name, fga, times=times, verticals=[ovsd], output_filename=new_file_name, vertical_positive='up')
                                    new_ts.add_variable(other, values=old_var[:], times=times, verticals=[ovsd], fillvalue=fillvalue, attributes=variable_attributes)
                                    depth_files.append(new_ts)
                        elif old_var.ndim <= 3 and (depth_values.size > 1 and not depth_variable and 'time' in old_var.dimensions):
                            if ovsd:
                                # An ADCP or profiling dataset, but this variable is measued at a single depth.
                                # Example: Bottom Temperature on an ADCP
                                # Skip things with a dimension over 3 (some beam variables like `brange`)
                                ts.add_variable(other, values=old_var[:], times=times, unlink_from_profile=True, verticals=[ovsd], fillvalue=fillvalue, attributes=variable_attributes)
                            else:
                                ts.add_variable(other, values=old_var[:], times=times, unlink_from_profile=True, fillvalue=fillvalue, attributes=variable_attributes)
                        else:
                            if 'time' in old_var.dimensions and old_var.ndim <= 3:
                                ts.add_variable(other, values=old_var[:], times=times, fillvalue=fillvalue, attributes=variable_attributes)
                            else:
                                ts.add_variable_object(old_var, dimension_map=dict(depth='z'), reduce_dims=True)

                    except BaseException:
                        logger.exception("Error processing variable {0} in {1}. Skipping it.".format(other, down_file))
        except KeyboardInterrupt:
            logger.info("Breaking out of Translate loop!")
            break
        except BaseException:
            logger.exception("Error. Skipping {0}.".format(down_file))
            continue
        finally:
            try:
                for df in depth_files:
                    del df
            except NameError:
                pass
            try:
                del ts
            except NameError:
                pass
            os.close(temp_fd)
            if os.path.isfile(temp_file):
                os.remove(temp_file)


if __name__ == "__main__":

    def valid_date(s):
        try:
            return datetime.strptime(s, "%Y-%m-%dT%H:%M").replace(tzinfo=pytz.utc)
        except ValueError:
            msg = "Not a valid date: '{0}'.".format(s)
            raise argparse.ArgumentTypeError(msg)

    parser = argparse.ArgumentParser()

    parser.add_argument('-o', '--output',
                        required=True,
                        help="Directory to output NetCDF files to",
                        nargs='?')
    parser.add_argument('-d', '--download',
                        action='store_true',
                        default=False,
                        help="Should we download a new set of files or use the files that have already been downloaded? \
                             Useful for debugging.  Downloaded files are never altered in any way so you can rerun \
                             the processing over and over without having to redownload any files.")
    parser.add_argument('-f', '--folder',
                        default=os.path.abspath(os.path.join(".", "download")),
                        help="Specify the folder location of NetCDF files you wish to translate. If this is used along with '--download', the files \
                              will be downloaded into this folder and then processed.  If used without the '--download' option, this is the \
                              location of the root folder you wish to translate into NetCDF files."),
    parser.add_argument('-p', '--projects',
                        help="Specific projects to process (optional).",
                        nargs='*')
    parser.add_argument('-l', '--files',
                        help="Specific files to process (optional).",
                        nargs='*')
    parser.add_argument('-c', '--csv_metadata_file',
                        help="CSV file to load metadata about each project from.  Defaults to 'project_metadata.csv'.",
                        default='project_metadata.csv',
                        nargs='?')
    parser.add_argument('-s', '--since',
                        help="Only process data modified since (utc) - format YYYY-MM-DDTHH:II",
                        type=valid_date)
    args = parser.parse_args()

    projects = args.projects
    if projects:
        projects = [ x.lower() for x in args.projects ]

    files = args.files
    if files:
        files = [ x.lower() for x in args.files ]
    main(args.output, args.folder, args.download, projects, os.path.realpath(args.csv_metadata_file), files, args.since)
