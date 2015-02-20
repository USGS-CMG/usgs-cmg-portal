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
from pyaxiom.netcdf.sensors import TimeSeries

import coloredlogs

# Log to stdout
logger = logging.getLogger()
coloredlogs.install(level=logging.DEBUG)
logger.setLevel(logging.DEBUG)

# Don't show the HTTP connection spam
requests_log = logging.getLogger("requests").setLevel(logging.WARNING)
crawler_log = logging.getLogger("thredds_crawler").setLevel(logging.INFO)
epic_log = logging.getLogger("epic2cf").setLevel(logging.INFO)
pyaxiom_log = logging.getLogger("pyaxiom").setLevel(logging.WARNING)

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


# Has depths and bottom variables: EUROSTRATAFORM, GLOBEC_GB
# Has mutiple depth variables (depth1, depth2, depth3, etc.): MOBILE_BAY /home/kwilcox/Development/usgs_cmg_woods_hole/download/3951tct-a.cdf


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
            datasets = Crawl(v['catalog_xml'], select=['.*-[A|a]{1}(?!1h)\.*.*'], skip=skips).datasets
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

    nc = None
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


def main(output, download_folder, do_download, projects, csv_metadata_file):

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
    try:
        os.makedirs(temp_folder)
    except OSError:
        pass  # Exists

    for down_file in downloaded_files:

        # For debugging
        #if os.path.basename(down_file) != "8451met-a.nc":
        #    continue

        nc = None
        try:
            temp_file = os.path.join(temp_folder, os.path.basename(down_file))
            shutil.copy(down_file, temp_file)

            if projects:
                tmpnc = netCDF4.Dataset(temp_file)
                project_name, _ = tmpnc.id.split("/")
                nc_close(tmpnc)
                if project_name.lower() not in projects:
                    # Skip this project!
                    continue

            # Cleanup to CF-1.6
            normalize_time(temp_file)
            normalize_epic_codes(temp_file)
            normalize_vectors(temp_file)
            normalize_units(temp_file)

            # Create list of variables that we want to save.
            station_id   = None
            station_name = None
            latitude     = None
            longitude    = None

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

            file_name = os.path.basename(down_file)
            output_directory = os.path.join(output, project_name)
            logger.info("Translating {0} into CF1.6 format: {1}".format(down_file, os.path.abspath(os.path.join(output_directory, file_name))))

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

            times           = nc.variables.get('time')[:]
            feature_name, _ = os.path.splitext(os.path.basename(down_file))
            # Get all depth values
            depth_variables = []
            for dv in nc.variables:
                depth_variables += [ x for x in nc.variables.get(dv).dimensions if 'depth' in x ]
            depth_variables = sorted(list(set(depth_variables)))
            depth_values = np.asarray([ nc.variables.get(x)[:] for x in depth_variables ]).flatten()

            ts = TimeSeries(output_directory, latitude, longitude, feature_name, file_global_attributes, times=times, verticals=depth_values, output_filename=file_name)

            v = []
            for other in sorted(nc.variables):  # Sorted for a reason... don't change!
                if other in coord_vars:
                    continue

                old_var = nc.variables.get(other)
                variable_attributes = { k : getattr(old_var, k) for k in old_var.ncattrs() }

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
                    depth_index = np.squeeze(np.where(depth_values == nc.variables.get(depth_variable[0])[:]))

                    # Find other variable names like this one
                    depth_indexes = [(other, depth_index)]
                    for search_var in sorted(nc.variables):
                        # If they have different depth dimension names we need to combine them into one variable
                        if search_var != other and search_var.split('_')[0] == new_var_name and \
                           depth_variable[0] != [ x for x in nc.variables[search_var].dimensions if 'depth' in x ][0]:
                            # Found a match at a different depth
                            search_depth_variable = [ x for x in nc.variables.get(search_var).dimensions if 'depth' in x ]
                            depth_index = np.squeeze(np.where(depth_values == nc.variables.get(search_depth_variable[0])[:]))
                            depth_indexes.append((search_var, depth_index))
                            logger.info("Combining '{}' with '{}' as '{}' (different variables at different depths but are the same parameter)".format(search_var, other, new_var_name))

                    values = np.ma.empty((times.size, len(depth_values)))
                    values.fill_value = fillvalue
                    values.mask = True
                    for nm, index in depth_indexes:
                        values[:, index] = np.squeeze(nc.variables.get(nm)[:])

                    # If we just have one index we want to use the original name
                    if len(depth_indexes) == 1:
                        # Just use the original variable name
                        new_var_name = other

                    # Create this one, should be the first we encounter for this type
                    ts.add_variable(new_var_name, values=values, times=times, fillvalue=fillvalue, attributes=variable_attributes)
                elif depth_variable and 'time' not in old_var.dimensions:
                    # elif (depth_variable and len(old_var.dimensions) == 1 and 'depth' == old_var.dimensions[0]) or \
                    # Metadata variable like bin distance
                    meta_var = ts.ncd.createVariable(other, old_var.dtype, ('z',), fill_value=fillvalue)
                    for k, v in variable_attributes.iteritems():
                        if k != '_FillValue':
                            setattr(meta_var, k, v)
                else:
                    values = old_var[:]
                    if len(old_var.dimensions) == 1 and old_var.dimensions[0] == 'time':
                        # Metadata variables like pitch, roll, record count, etc.
                        ts.add_variable(other, values=values, times=times, unlink_from_profile=True, fillvalue=fillvalue, attributes=variable_attributes)
                    elif depth_values.size > 1:
                        # No Z variables in a profile dataset, aka Bottom Temperature
                        ts.add_variable(other, values=values, times=times, verticals=[old_var.sensor_depth], unlink_from_profile=True, fillvalue=fillvalue, attributes=variable_attributes)
                    else:
                        ts.add_variable(other, values=values, times=times, fillvalue=fillvalue, attributes=variable_attributes)

                ts.ncd.sync()
            ts.ncd.close()

        except BaseException:
            logger.exception("Error. Skipping {0}.".format(down_file))
            continue
        finally:
            nc_close(nc)
            os.remove(temp_file)

    shutil.rmtree(temp_folder, ignore_errors=True)


if __name__ == "__main__":
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
    parser.add_argument('-c', '--csv_metadata_file',
                        help="CSV file to load metadata about each project from.  Defaults to 'project_metadata.csv'.",
                        default='project_metadata.csv',
                        nargs='?')
    args = parser.parse_args()

    projects = args.projects
    if projects:
        projects = map(lambda x: x.lower(), args.projects)
    main(args.output, args.folder, args.download, projects, os.path.realpath(args.csv_metadata_file))
