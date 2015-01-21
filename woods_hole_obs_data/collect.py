#!python
# coding=utf-8

import os
import csv
import shutil
import logging
import requests
import argparse
from glob import glob
from datetime import datetime

import epic2cf
from epic2cf.data import *
import netCDF4
import numpy as np

from thredds_crawler.crawl import Crawl
from pytools.netcdf.sensors.create import create_timeseries_file

import coloredlogs

# Log to stdout
logger = logging.getLogger()
coloredlogs.install(level=logging.INFO)
logger.setLevel(logging.INFO)
fh = logging.FileHandler('usgs_cmg.log')
fh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - [%(levelname)s] %(message)s')
fh.setFormatter(formatter)
logger.addHandler(fh)

# Don't show the HTTP connection spam
requests_log = logging.getLogger("requests").setLevel(logging.WARNING)
crawler_log = logging.getLogger("thredds_crawler").setLevel(logging.INFO)
epic_log = logging.getLogger("epic2cf").setLevel(logging.INFO)

IGNORABLE_CODES = location_codes + time_codes + generic_codes + voltage_codes

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
                                                      convert=lambda x: x/1000.,
                                                      units='kg/m^3')),
    'BGAPE'     : dict(epic_code=None, overrides=dict(standard_name='mass_concentration_of_phycoerythrin_expressed_as_chlorophyll_in_sea_water',
                                                      units='kg/m^3',
                                                      convert=lambda x: x/1000000.)),
    'turbidity' : dict(epic_code=980),
    'Qs_133'    : dict(epic_code=None, overrides=dict(standard_name='net_downward_shortwave_flux_in_air',
                                                      original_units='w/milli-angstrom^2',
                                                      units='w/m^2',
                                                      convert=lambda x: x/1e13)),
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
                                                                                convert=lambda x: x/1000.)),
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
                                                                                convert=lambda x: x*0.6894757293)),
}

global_attributes = {
    'naming_authority':         'gov.usgs.cmgp',
    'source':                   'USGS',
    'institution':              'USGS Woods Hole Coastal and Marine Science Center',
    'project':                  'Coastal and Marine Geology Program',
    'keywords':                 'Oceans > Ocean Pressure > Water Pressure, Oceans > Ocean Temperature > Water Temperature, Oceans > Salinity/Density > Conductivity, Oceans > Salinity/Density > Salinity',
    'keywords_vocabulary':      'CMD Science Keywords',
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

coord_vars      = ['time', 'time2', 'time_cf', 'old_time', 'depth', 'depth002', 'depth003', 'depth004', 'depth005', 'lat', 'lon']


# Has depths: EUROSTRATAFORM, GLOBEC_GB
# Has mutiple depth variables (depth1, depth2, depth3, etc.): MOBILE_BAY /home/kwilcox/Development/usgs_cmg_woods_hole/download/3951tct-a.cdf
# ADCP with a bottom temperature variable that isn't varying by depth: /home/kwilcox/Development/usgs_cmg_woods_hole/download/6551adc-a.nc


def nc_close(nc):
    if nc is not None:
        try:
            nc.sync()
            nc.close()
        except RuntimeError:
            pass


def download(folder, project_metadata):

    # Use thredds_crawler to find DAP endpoints of the RAW data.
    total_datasets = []
    skips = Crawl.SKIPS + ['.*OTHER.*', '.*ancillary.*']

    try:
        for k, v in project_metadata.items():
            datasets = Crawl(v['catalog_xml'], select=['.*-[A|a]+\..*'], skip=skips).datasets
            logger.info("Found {0} datasets in {1}!".format(len(datasets), k))
            total_datasets += datasets
        logger.info("Found {0} TOTAL datasets!".format(len(total_datasets)))
    except KeyboardInterrupt:
        logger.info("Breaking out of crawling loop.")
        total_datasets = []

    shutil.rmtree(folder, ignore_errors=True)
    os.makedirs(folder)

    # Save datasets to download directory
    saved_files = []
    for num, d in enumerate(total_datasets):
        try:
            http_url = next(s["url"] for s in d.services if s["service"].lower() == "httpserver")
        except StopIteration:
            logger.error("No HTTPServer endpoint found, skipping")
            continue

        # Make download folder
        save_file = os.path.join(folder, d.name)
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
            logger.info("Could not download... error with HTTP endpoint.  Skipping.")
            continue

        # Try to open file, if it fails, writing failed.
        try:
            nc = netCDF4.Dataset(save_file, 'a')
            name, _ = os.path.splitext(d.name)
            project_name = http_url.split("/")[-2]
            nc.id = "{0}/{1}".format(project_name, name)
        except BaseException:
            os.remove(save_file)
            raise
        else:
            logger.info("{!s} saved ({!s}/{!s})".format(d.name, num + 1, len(total_datasets)))
            saved_files.append(save_file)
        finally:
            nc_close(nc)

    return saved_files


def normalize_epic_codes(netcdf_file):
    nc = netCDF4.Dataset(netcdf_file, 'a')
    try:
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
                        setattr(nc_var, k, d)

            if hasattr(nc_var, 'long_name'):
                if not hasattr(nc_var, 'epic_code') or hasattr(nc_var, 'epic_code') and nc_var.epic_code in IGNORABLE_CODES:
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
                                setattr(nc_var, k, d)

            if hasattr(nc_var, "epic_code") and nc_var.epic_code:
                try:
                    int(nc_var.epic_code)
                except ValueError:
                    logger.debug("No EPIC code specified on {0}".format(v))
                else:
                    attribs = epic2cf.mapping.get(int(nc_var.epic_code))
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
    except BaseException:
        logger.exception("Error.")
        raise
    finally:
        nc_close(nc)


def normalize_vectors(netcdf_file):
    nc = netCDF4.Dataset(netcdf_file, 'a')
    try:
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
        if east is not None and north is not None and 'eastward_sea_water_velocity' not in std_names and 'northward_sea_water_velocity' not in std_names:
            # We have vectors... create the speed and direction variables
            speed = np.sqrt(np.square(east[:]) + np.square(north[:]))
            direction = np.degrees(np.arctan2(north[:], east[:]))

            east_fill_value = east._FillValue if hasattr(east, '_FillValue') else np.nan
            spd = nc.createVariable('CS_300', 'f4', east.dimensions, fill_value=east_fill_value)
            spd.standard_name = 'sea_water_speed'
            spd.epic_code     = 300
            spd[:] = speed

            drc = nc.createVariable('CD_310', 'f4', east.dimensions, fill_value=east_fill_value)
            drc.standard_name = 'direction_of_sea_water_velocity'
            drc.epic_code     = 310
            drc[:] = direction

    except BaseException:
        logger.exception("Error")
        raise
    finally:
        nc_close(nc)


def normalize_units(netcdf_file):
    nc = netCDF4.Dataset(netcdf_file, 'a')
    try:
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
    except BaseException:
        logger.exception("Error")
        raise
    finally:
        nc_close(nc)


def normalize_time(netcdf_file):
    epoch_units       = 'seconds since 1970-01-01T00:00:00Z'
    millisecond_units = 'milliseconds since 1858-11-17T00:00:00Z'

    try:
        nc = netCDF4.Dataset(netcdf_file, 'a')
        # Signell said this works, any problems and we can all blame him!
        time_data = netCDF4.num2date((np.int64(nc.variables['time'][:])-2400001)*3600*24*1000 + nc.variables['time2'][:].__array__(), units=millisecond_units)
        nc.renameVariable("time", "old_time")
        nc.sync()

        time = nc.createVariable('time', 'f8', ('time'))
        time.units          = epoch_units
        time.standard_name  = "time"
        time.long_name      = "time of measurement"
        time.calendar       = "gregorian"
        time[:] = netCDF4.date2num(time_data, units=epoch_units).round()
    except TypeError:
        raise TypeError("The TIME variable can not be converted.  Number of timesteps in original file: {!s}".format(nc.variables.get('time').size))
    except BaseException:
        raise
    finally:
        nc_close(nc)


def main(output, download_folder, do_download, output_format, projects, csv_metadata_file):

    project_metadata = dict()
    with open(csv_metadata_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            project_name = row['project_name']
            if projects and project_name.lower() not in projects:
                # Skip projects if a subset was defined
                continue
            project_metadata[project_name] = dict()
            for k, v in row.items():
                project_metadata[project_name][k] = v

    if do_download:
        try:
            downloaded_files = download(download_folder, project_metadata)
        except KeyboardInterrupt:
            downloaded_files = []
    else:
        downloaded_files = glob(os.path.join(download_folder, "*"))

    temp_folder = os.path.abspath(os.path.join(".", "temp"))
    shutil.rmtree(temp_folder, ignore_errors=True)
    os.makedirs(temp_folder)

    for down_file in downloaded_files:

        # For debugging
        #if os.path.basename(down_file) != "5881adc-a.nc":
        #    continue

        temp_file = os.path.join(temp_folder, os.path.basename(down_file))
        shutil.copy(down_file, temp_file)

        if projects:
            tmpnc = netCDF4.Dataset(temp_file)
            project_name, _ = tmpnc.id.split("/")
            nc_close(tmpnc)
            if project_name.lower() not in projects:
                # Skip this project!
                continue

        nc = None
        try:
            # Cleanup to CF-1.6
            normalize_time(temp_file)
            normalize_epic_codes(temp_file)
            normalize_vectors(temp_file)
            normalize_units(temp_file)

            # Create list of variables that we want to save.
            station_id   = None
            station_urn  = None
            station_name = None
            latitude     = None
            longitude    = None
            starting     = None
            ending       = None

            nc = netCDF4.Dataset(temp_file)

            # Default station_id
            project_name, _ = nc.id.split("/")
            # Now try to come up with a better one.
            if hasattr(nc, 'MOORING') and hasattr(nc, 'id'):
                mooring_id = str(nc.MOORING).replace(':', '').strip()
                station_id = "{0}_{1}".format(project_name, mooring_id[0:3]).lower()
                station_name = "{0} ({1})".format(project_name, mooring_id[0:3])
            else:
                try:
                    # Mooring ID is the first three numbers of the file
                    station_id = int(os.path.basename(down_file)[0:3])
                    station_id = "{0}_mooring_{0}".format(project_name, station_id)
                    station_name = "{0} Mooring ({0})".format(project_name, station_id)
                except BaseException:
                    logger.error("Could not create a suitable station_id. Skipping {0}.".format(down_file))
                    continue

            try:
                latitude  = nc.variables.get("lat")[0]
                longitude = nc.variables.get("lon")[0]
            except IndexError:
                latitude  = nc.variables.get("lat")[:]
                longitude = nc.variables.get("lon")[:]
            starting  = datetime.utcfromtimestamp(nc.variables.get("time")[0])
            ending    = datetime.utcfromtimestamp(nc.variables.get("time")[-1])

            station_urn = "urn:ioos:station:{0}:{1}".format('gov.usgs.cmgp', station_id).lower()

            logger.info("FILE: {0}".format(down_file))
            #logger.info("STATION: {0}".format(station_urn))

            def create_axiom():
                data_variables  = list()
                other_variables = list()
                for v in nc.variables:
                    if v in coord_vars:
                        continue  # Skip coordinate variables`
                    nc_var = nc.variables.get(v)
                    try:
                        var_epic_code = int(nc_var.epic_code)
                    except (AttributeError, ValueError):
                        var_epic_code = None
                    if hasattr(nc_var, "cf_role"):
                        other_variables.append(v)
                    elif hasattr(nc_var, "standard_name") and nc_var.standard_name.lower() == "time":
                        pass  # Skip time variables
                    elif hasattr(nc_var, "standard_name") and nc_var.standard_name.lower() in ["latitude", "longitude"]:
                        pass  # Skip lat/lon variables
                    elif hasattr(nc_var, "standard_name") and nc_var.standard_name.lower() == "surface_altitude":
                        other_variables.append(v)
                    elif hasattr(nc_var, "axis"):
                        other_variables.append(v)
                    elif var_epic_code and var_epic_code in metadata_codes:
                        other_variables.append(v)  # Add metadata variables to the file.  This later checks the dimensions match the data variable.
                    elif var_epic_code and var_epic_code in IGNORABLE_CODES:
                        continue  # Skip ignorable variables
                    elif v == "bindist":
                        other_variables.append(v)
                    elif hasattr(nc_var, "standard_name") and nc_var.standard_name.lower() in ['northward_sea_water_velocity', 'eastward_sea_water_velocity']:
                        continue  # We created speed/direction variables, so skip these
                    elif hasattr(nc_var, "standard_name"):
                        if hasattr(nc_var, "cell_methods"):
                            data_variables.append((v, nc_var.standard_name, nc_var.cell_methods, ))
                        else:
                            data_variables.append((v, nc_var.standard_name, None))
                    else:
                        if hasattr(nc_var, 'long_name'):
                            if var_epic_code:
                                logger.warning("Skipping {0}, no standard_name attribute.  'long_name' is {1}.  epic_code is {2}.".format(v, nc_var.long_name.strip(), nc_var.epic_code))
                            else:
                                logger.warning("Skipping {0}, no standard_name attribute.  'long_name' is {1}.".format(v, nc_var.long_name.strip()))
                        else:
                            logger.warning("Skipping {0}, no standard_name or long_name attribute.".format(v))
                        continue

                for dv, std, cm in data_variables:
                    try:
                        logger.debug("Exporting Axiom format: {0}".format(dv))
                        #logger.info("Creating file with the following variables: {!s}".format(other_variables + [dv]))
                        file_name = "{0}_{1}_TO_{2}.nc".format(dv, starting.strftime("%Y-%m-%dT%H:%MZ"), ending.strftime("%Y-%m-%dT%H:%MZ"))
                        sensor_urn = "{0}:{1}".format(station_urn.replace("station", "sensor"), std)
                        if cm is not None:
                            split_cms = cm.split(":")
                            cms = [ "{0}:{1}".format(c[0].strip(), c[1].strip()) for c in zip(split_cms[0::2], split_cms[1::2]) ]
                            sensor_urn_with_cellmethods = "{0}#cell_methods={1}".format(sensor_urn, ",".join(cms))
                            output_directory = os.path.join(output, sensor_urn_with_cellmethods)
                        else:
                            output_directory = os.path.join(output, sensor_urn)

                        if not os.path.isdir(output_directory):
                            os.makedirs(output_directory)

                        file_global_attributes = { k : getattr(nc, k) for k in nc.ncattrs() }
                        file_global_attributes.update(global_attributes)
                        file_global_attributes['id']               = station_id
                        file_global_attributes['title']            = station_name
                        if project_name in project_metadata:
                            for k, v in project_metadata[project_name].items():
                                if v and k.lower() not in ['id', 'title', 'catalog_xml', 'project_name']:
                                    file_global_attributes[k] = v

                        ts = nc.variables.get("time")[:]

                        depth_variable = [ x for x in nc.variables.get(dv).dimensions if 'depth' in x ]
                        if len(depth_variable) >= 1:
                            depth_variable = depth_variable[0]
                        else:
                            depth_variable = None

                        if depth_variable:
                            zs = nc.variables.get(depth_variable)[:]
                            times     = np.ma.repeat(ts, zs.size)
                            verticals = np.ma.ravel(np.ma.repeat([zs], ts.size, axis=0))
                            values    = np.squeeze(nc.variables.get(dv)[:]).flatten()
                        else:
                            sensor_depth = 0
                            if hasattr(nc.variables.get(dv), 'sensor_depth'):
                                sensor_depth = nc.variables.get(dv).sensor_depth
                            times     = ts
                            verticals = np.ma.repeat(sensor_depth, ts.size)
                            values    = np.squeeze(nc.variables.get(dv)[:]).flatten()

                        assert values.size == verticals.size == times.size
                        create_timeseries_file(output_directory=output_directory, latitude=latitude, longitude=longitude, full_station_urn=station_urn, full_sensor_urn=sensor_urn, global_attributes=file_global_attributes, attributes=nc.variables.get(dv).__dict__, output_filename=file_name, times=times, verticals=verticals, values=values)

                        new_nc = netCDF4.Dataset(os.path.join(output_directory, file_name), 'a')
                        new_nc.setncattrs(file_global_attributes)
                        for other in other_variables:
                            old_var = nc.variables.get(other)

                            # Get new variable name
                            variable_name = sensor_urn.split(":")[-1]
                            new_var = new_nc.variables.get(variable_name)

                            other_var = new_nc.createVariable(other, old_var.dtype, new_var.dimensions)
                            for k in old_var.ncattrs():
                                other_var.setncattr(k, old_var.getncattr(k))
                            try:
                                other_var[:] = nc.variables.get(other)[:]
                                logger.debug("Adding metadata variable: {0}".format(other))
                            except BaseException:
                                logger.debug("Skipping: {0}.  It had more dimensions than the core variable '{1}'.".format(other, dv))
                                continue
                            new_nc.sync()
                        nc_close(new_nc)

                    except BaseException:
                        logger.exception("Error. Skipping {0} in {1}.".format(down_file, dv))
                        if os.path.exists(os.path.join(output_directory, file_name)):
                            os.remove(os.path.join(output_directory, file_name))
                        continue

            def create_cf16():
                file_name = os.path.basename(down_file)
                logger.debug("Exporting CF1.6 format: {0}".format(file_name))
                #logger.info("Creating file with the following variables: {!s}".format(other_variables + [dv]))
                file_name = os.path.basename(down_file)
                output_directory = os.path.join(output, project_name)

                if not os.path.isdir(output_directory):
                    os.makedirs(output_directory)

                file_global_attributes = { k : getattr(nc, k) for k in nc.ncattrs() }
                file_global_attributes.update(global_attributes)
                file_global_attributes['id']               = station_id
                file_global_attributes['title']            = station_name
                if project_name in project_metadata:
                    for k, v in project_metadata[project_name].items():
                        if v and k.lower() not in ['id', 'title', 'catalog_xml', 'project_name']:
                            file_global_attributes[k] = v

                # Compute additional metadata
                file_global_attributes['time_coverage_start'] = starting.isoformat()
                file_global_attributes['time_coverage_end'] = ending.isoformat()
                # duration (ISO8601 format)
                file_global_attributes['time_coverage_duration'] = "P%sS" % unicode(int(round((ending - starting).total_seconds())))
                # resolution (ISO8601 format)
                # subtract adjacent times to produce an array of differences, then get the most common occurance
                unique_times = np.unique(nc.variables.get('time')[:])
                diffs = unique_times[1:] - unique_times[:-1]
                uniqs, inverse = np.unique(diffs, return_inverse=True)
                time_diffs = diffs[np.bincount(inverse).argmax()]
                file_global_attributes['time_coverage_resolution'] = "P%sS" % unicode(int(round(time_diffs)))

                new_nc = netCDF4.Dataset(os.path.join(output_directory, file_name), 'w')
                new_nc.createDimension("time", nc.variables.get('time').size)

                # Location
                lat = new_nc.createVariable("latitude", "f8")
                lat.units           = "degrees_north"
                lat.standard_name   = "latitude"
                lat.long_name       = "sensor latitude"
                lat[:] = latitude
                file_global_attributes['geospatial_lat_min'] = latitude
                file_global_attributes['geospatial_lat_max'] = latitude

                lon = new_nc.createVariable("longitude", "f8")
                lon.units           = "degrees_east"
                lon.standard_name   = "longitude"
                lon.long_name       = "sensor longitude"
                lon[:] = longitude
                file_global_attributes['geospatial_lon_min'] = longitude
                file_global_attributes['geospatial_lon_max'] = longitude

                # Station name
                feature_name, _ = os.path.splitext(os.path.basename(down_file))
                new_nc.createDimension("site_name", len(feature_name))
                site = new_nc.createVariable("site", "S1", ("site_name",))
                site.cf_role = "timeseries_id"
                site.standard_site = "station_id"
                site.long_name = "Identifier for each site"
                site[:] = list(feature_name)

                # Metadata variables
                crs = new_nc.createVariable("crs", "i4")
                crs.long_name           = "http://www.opengis.net/def/crs/EPSG/0/4326"
                crs.grid_mapping_name   = "latitude_longitude"
                crs.epsg_code           = "EPSG:4326"
                crs.semi_major_axis     = float(6378137.0)
                crs.inverse_flattening  = float(298.257223563)

                platform = new_nc.createVariable("platform", "i4")
                platform.ioos_code  = station_urn
                platform.short_name = global_attributes.get("title", station_urn)
                platform.long_name  = global_attributes.get("description", station_urn)

                # Get all depth values
                depth_variables = []
                for dv in nc.variables:
                    depth_variables += [ x for x in nc.variables.get(dv).dimensions if 'depth' in x ]
                depth_variables = sorted(list(set(depth_variables)))
                depth_values = np.asarray([ nc.variables.get(x)[:] for x in depth_variables ]).flatten()

                depth_var = None
                if depth_values.size == 1:
                    depth_var = new_nc.createVariable('depth', 'f4')
                elif depth_values.size > 0:
                    new_nc.createDimension('depth', depth_values.size)
                    depth_var = new_nc.createVariable('depth', 'f4', ('depth', ))

                if depth_var is not None:
                    for k in nc.variables.get(depth_variables[0]).ncattrs():
                        if k not in ['_FillValue', 'missing_value']:
                            depth_var.setncattr(k, nc.variables.get(depth_variables[0]).getncattr(k))
                    depth_var.axis = "Z"
                    depth_var.positive = "down"
                    depth_var[:] = depth_values

                def create_variable_from_old(netcdf, old_var, new_var_name=None, new_dims=None):
                    if hasattr(old_var, "_FillValue"):
                        new_var = netcdf.createVariable(new_var_name or old_var.name, old_var.dtype, new_dims or old_var.dimensions, fill_value=old_var._FillValue, zlib=True)
                    else:
                        new_var = netcdf.createVariable(new_var_name or old_var.name, old_var.dtype, new_dims or old_var.dimensions, zlib=True)
                    for k in old_var.ncattrs():
                        try:
                            if k not in ['_FillValue', 'missing_value']:
                                new_var.setncattr(k, old_var.getncattr(k))
                        except AttributeError:
                            # Already has this attribute.  Bug in netCDF C library < 4.2
                            pass
                    # Remove whitespace from long_name
                    if hasattr(new_var, 'long_name'):
                        new_var.long_name = new_var.long_name.strip()
                    return new_var

                v = []
                for other in nc.variables:
                    old_var = nc.variables.get(other)
                    new_dims = []
                    if 'time' in old_var.dimensions:
                        new_dims += ['time']
                    if depth_values.size > 1 and (old_var.ndim > 3 or 'depth' in old_var.dimensions):
                        new_dims += ['depth']

                    new_coords = None
                    if other in ['lat', 'lon', 'old_time', 'time2', 'time_cf', 'depth', 'depth002', 'depth003', 'depth004', 'depth005']:
                        continue
                    elif 'time' in old_var.dimensions and 'depth' in old_var.dimensions:
                        new_coords = 'time depth latitude longitude'
                    elif depth_values.size > 1 and 'time' in old_var.dimensions and 'depth' not in old_var.dimensions and other != 'time':
                        # Special case bottom transducer temp
                        new_coords = 'time depth latitude longitude'
                    else:
                        # We just copy these variables...
                        logger.info('Copying {}'.format(other))
                        new_var = create_variable_from_old(new_nc, old_var, other, new_dims)
                        new_var[:] = old_var[:]
                        continue

                    old_depth_variable = None
                    if len(depth_variables) > 1:
                        old_depth_variable = [ x for x in old_var.dimensions if 'depth' in x ][0]
                        # Normalize the variable name
                        new_var_name = other.split("_")[0]
                        if not nc.variables.get(new_var_name):
                            # Try to match
                            new_var_name = other.split("_")[0] + "_1"  # sorted([ x for x in nc.variables if hasattr(nc.variables.get(x), 'long_name') and nc.variables.get(x).long_name == ])[0]
                            if not nc.variables.get(new_var_name):
                                # Fall back to the original name.
                                new_var_name = other
                                logger.info('Could not match "{}" to any additional variables at varying depths'.format(other))

                        if new_nc.variables.get(new_var_name):
                            new_var = new_nc.variables.get(new_var_name)
                        else:
                            new_var = create_variable_from_old(new_nc, old_var, new_var_name, new_dims)
                    else:
                        new_var_name = other
                        new_var = create_variable_from_old(new_nc, old_var, new_var_name, new_dims)

                    # Set 'crs' and 'coordinates' attributes on data variables
                    new_var.coordinates = new_coords
                    new_var.grid_mapping = 'crs'

                    # Figure out the depth index for this variable
                    if old_depth_variable:
                        if old_depth_variable == 'depth':
                            depth_index = 0
                        else:
                            try:
                                depth_index = int(old_depth_variable[-3:]) - 1
                            except ValueError:
                                logger.exception("Could not pull the correct depth index!!")
                                continue

                        logger.info('Copying {} {} into {} {} at depth index {}'.format(other, old_var.dimensions, new_var_name, new_var.dimensions, depth_index))
                        new_var[:, depth_index] = old_var[:]
                    else:
                        logger.info("Copying {} into {}".format(other, new_var_name))
                        try:
                            new_var[:] = old_var[:]
                        except BaseException:
                            logger.exception("Error!  Skipping {}".format(new_var_name))
                    new_nc.sync()

                # Handle temperature at ADCP transducer heads
                sensor_depth = None
                if hasattr(nc, 'inst_depth'):
                    sensor_depth = nc.inst_depth
                elif hasattr(nc, 'nominal_sensor_depth'):
                    sensor_depth = nc.nominal_sensor_depth

                if sensor_depth is not None:
                    inst_depth = new_nc.createVariable('sensor_depth', 'f4')
                    inst_depth.units = 'm'
                    inst_depth.standard_name = 'surface_altitude'
                    inst_depth.long_name = 'sensor depth below datum'
                    inst_depth.positive = 'up'
                    inst_depth.datum = 'NAVD88'
                    inst_depth[:] = sensor_depth

                unique_verticals = np.unique(new_nc.variables.get('depth')[:])
                file_global_attributes['geospatial_vertical_positive'] = "down"
                file_global_attributes['geospatial_vertical_min']      = np.min(unique_verticals)
                file_global_attributes['geospatial_vertical_max']      = np.max(unique_verticals)
                if unique_verticals.size > 1:
                    vertical_diffs = unique_verticals[1:] - unique_verticals[:-1]
                    file_global_attributes['geospatial_vertical_resolution'] = " ".join(map(unicode, list(vertical_diffs)))
                    file_global_attributes['featureType'] = 'timeSeriesProfile'
                elif unique_verticals.size == 1:
                    file_global_attributes['featureType'] = 'timeSeries'

                new_nc.setncatts(file_global_attributes)
                new_nc.close()

            if output_format.lower() == 'axiom':
                create_axiom()
            elif output_format.lower() == 'cf16':
                create_cf16()

        except BaseException:
            logger.exception("Error. Skipping {0}.".format(down_file))
            continue
        finally:
            nc_close(nc)
            os.remove(temp_file)

    shutil.rmtree(temp_folder, ignore_errors=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('format', action='store', choices=('axiom', 'cf16',), help="Which type of file to produce. You most likely want 'cf16'.")
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
    parser.add_argument('-c', '--csv_metadata_file',
                        help="CSV file to load metadata about each project from.  Defaults to 'project_metadata.csv'.",
                        default='project_metadata.csv',
                        nargs='?')
    args = parser.parse_args()

    projects = args.projects
    if projects:
        projects = map(lambda x: x.lower(), args.projects)
    main(args.output, args.folder, args.download, args.format, projects, os.path.realpath(args.csv_metadata_file))
