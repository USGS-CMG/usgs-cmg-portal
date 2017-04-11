#!python
# coding=utf-8

import os
import logging
import argparse
from glob import glob
from datetime import datetime

import pytz
import requests
import netCDF4
from epic2cf.data import metadata_codes, voltage_codes, location_codes, time_codes
from thredds_crawler.crawl import Crawl

from pyaxiom.netcdf.sensors.timeseries import TimeSeries, get_dataframe_from_variable
from pyaxiom.netcdf.dataset import EnhancedDataset
from pyaxiom.utils import urnify
from pyaxiom.urn import IoosUrn

import coloredlogs

# Log to stdout
logger = logging.getLogger()
coloredlogs.install(level=logging.INFO)
logger.setLevel(logging.INFO)

# Log to a file
df = logging.FileHandler('dv_errors.log', mode='w', encoding='utf-8')
df.setLevel(logging.WARNING)
df.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(df)

# Don't show the HTTP connection spam
pyaxiom_log = logging.getLogger("pyaxiom").setLevel(logging.INFO)
requests_log = logging.getLogger("requests").setLevel(logging.WARNING)
crawler_log = logging.getLogger("thredds_crawler").setLevel(logging.INFO)


def nc_close(nc):
    if nc is not None:
        try:
            nc.close()
        except RuntimeError:
            pass


def download(folder, projects, filesubset, since):

    # Use thredds_crawler to find DAP endpoints of the CF-1.6 data.
    skips = Crawl.SKIPS
    if projects:
        skips += ['^(?!{}).*^(?!.*\.(cdf|nc)).*$'.format('|'.join(projects))]

    catalog = 'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/CF-1.6/catalog.html'

    try:
        datasets = Crawl(catalog, select=['.*\.(cdf|nc)'], skip=skips, after=since).datasets
        logger.info("Found {0} TOTAL datasets!".format(len(datasets)))
    except KeyboardInterrupt:
        logger.info("Breaking out of crawling loop.")
        datasets = []

    try:
        os.makedirs(folder)
    except OSError:
        pass

    # Save datasets to download directory
    saved_files = []
    for num, d in enumerate(datasets):

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
            nc = netCDF4.Dataset(save_file)
        except BaseException:
            os.remove(save_file)
        else:
            logger.info("{!s} saved ({!s}/{!s})".format(d.name, num + 1, len(datasets)))
            saved_files.append(save_file)
        finally:
            nc_close(nc)

    return saved_files


def main(output, download_folder, projects, do_download, filesubset=None, since=None):

    if do_download:
        try:
            downloaded_files = download(download_folder, projects, filesubset, since)
        except KeyboardInterrupt:
            logger.exception('Interpted downloading datasets from THREDDS')
            downloaded_files = []
    else:
        downloaded_files = glob(os.path.join(download_folder, '**', '*'))
        if since is not None:
            def should_keep(d):
                modt = datetime.utcfromtimestamp(os.path.getmtime(d)).replace(tzinfo=pytz.utc)
                return modt >= since
            downloaded_files = [ dl for dl in downloaded_files if should_keep(dl) ]

    epic_skips = metadata_codes + voltage_codes + location_codes + time_codes

    downloaded_files = sorted(downloaded_files)

    # Take the downloaded_files and split them up into arrays of related files.
    # This is needed because some of the files need to be combined... they
    # represent that same station/mooring/variables, but at different depths
    i = 0
    combinations = []
    while i < len(downloaded_files):

        combo = []
        nc_file = os.path.abspath(downloaded_files[i])
        combo.append(nc_file)

        try:

            if filesubset is not None:
                if os.path.basename(nc_file).lower() not in filesubset:
                    # aka "9631ecp-a.nc"
                    # Skip this file!
                    continue

            with EnhancedDataset(nc_file) as tmpnc:

                if projects:
                    if hasattr(tmpnc, 'original_folder') and tmpnc.original_folder.upper() not in projects:
                        continue

                logger.info("Scanned {}".format(nc_file))
                # Now search for files that are of the same var, but with a different depth
                thisbase = os.path.basename(nc_file).lower().split("_d")[0]

                for j in range(i + 1, len(downloaded_files)):
                    nextout = os.path.abspath(downloaded_files[j])
                    nextbase = os.path.basename(nextout).lower().split("_d")[0]
                    if thisbase == nextbase:
                        # Found a match
                        logger.info("Scanned {}".format(nextout))
                        combo.append(nextout)
                        # Now skip file because we added it already
                        i += 1
                    else:
                        # Doesn't match, move on to the next outfile
                        break
                # Add to combinations so it is processed
                combinations.append(combo)
        except BaseException:
            logger.exception("Error. Skipping {0}.".format(nc_file))
            continue
        finally:
            i += 1

    # Now iterate over each set of files and combine as necessary
    for c in combinations:
        dataframes_to_create = dict()
        for f in c:
            try:
                with EnhancedDataset(f) as ncd:
                    for var in ncd.get_variables_by_attributes(coordinates=lambda v: v is not None):
                        if not hasattr(var, 'standard_name'):
                            logger.warning("{}: Skipping variable {} because it has no standard_name".format(f, var.name))
                            continue
                        if hasattr(var, 'epic_code') and var.epic_code in epic_skips:
                            logger.warning("{}: Skipping metadata variable {}".format(f, var.standard_name))
                            continue

                        df = get_dataframe_from_variable(ncd, var)

                        if var.name in dataframes_to_create:
                            logger.info("Combining variable {}".format(var.name))
                            old_df = dataframes_to_create[var.name]['frame']
                            dataframes_to_create[var.name]['frame'] = old_df.combine_first(df)
                        else:
                            logger.info("New variable {}".format(var.name))
                            df_dict = dict(frame=df,
                                           varname=var.name,
                                           ncfile=f)
                            dataframes_to_create[var.name] = df_dict
            except BaseException:
                logger.exception("Error. Skipping {0}.".format(f))
                continue

        # Create a file for each dataframe
        for varname, creation in dataframes_to_create.items():
            create_file(output, creation['ncfile'], creation['varname'], creation['frame'])


def create_file(output, ncfile, varname, df):
    with EnhancedDataset(ncfile) as ncd:
        var = ncd[varname]

        latitude  = ncd.get_variables_by_attributes(standard_name='latitude')[0][:]
        longitude = ncd.get_variables_by_attributes(standard_name='longitude')[0][:]
        project = ncd.original_folder
        feature_name = '{}_{}'.format(project, ncd.MOORING).lower()

        station_urn = IoosUrn(authority=ncd.naming_authority, label=feature_name, asset_type='station').urn

        discriminant = ncd.id.replace('-', '_')
        output_filename = '{0}_{1}-{2}_{3}_TO_{4}.nc'.format(feature_name, var.name, discriminant, df['time'].min().strftime("%Y%m%dT%H%M%SZ"), df['time'].max().strftime("%Y%m%dT%H%M%SZ"))
        output_directory = os.path.join(output, feature_name)

        if not os.path.isdir(output_directory):
            os.makedirs(output_directory)

        file_global_attributes = { k : getattr(ncd, k) for k in ncd.ncattrs() }
        # original_folder is the project name
        file_global_attributes.update(dict(title='{} - {}'.format(project, ncd.MOORING),
                                           id=feature_name))

        variable_attributes = { k : getattr(var, k) for k in var.ncattrs() }
        # Add the specific sensor as a discriminant
        variable_attributes.update(dict(discriminant=discriminant))

        fillvalue = -9999.9
        if hasattr(var, "_FillValue"):
            fillvalue = var._FillValue

        vertical_datum = None
        if 'crs' in ncd.variables and hasattr(ncd.variables['crs'], 'vertical_datum'):
            vertical_datum = ncd.variables['crs'].vertical_datum

        ts = TimeSeries.from_dataframe(df, output_directory, output_filename, latitude, longitude, station_urn, file_global_attributes, var.standard_name, variable_attributes, sensor_vertical_datum=vertical_datum, fillvalue=fillvalue, vertical_axis_name='height', vertical_positive='down')
        ts.add_instrument_variable(variable_name=var.standard_name)


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
    parser.add_argument('-f', '--folder',
                        required=True,
                        help="Folder that the CF1.6 NetCDF files from USGS CMG THREDDS server were downloaded to",
                        nargs='?')
    parser.add_argument('-d', '--download',
                        action='store_true',
                        default=False,
                        help="Should we download a new set of files or use the files that have already been downloaded? \
                             Useful for debugging.  Downloaded files are never altered in any way so you can rerun \
                             the processing over and over without having to redownload any files.")
    parser.add_argument('-p', '--projects',
                        help="Specific projects to process (optional).",
                        nargs='*')
    parser.add_argument('-l', '--files',
                        help="Specific files to process (optional).",
                        nargs='*')
    parser.add_argument('-s', '--since',
                        help="Only process data modified since (utc) - format YYYY-MM-DDTHH:II",
                        type=valid_date)
    args = parser.parse_args()

    projects = args.projects
    if projects:
        projects = [ x.upper() for x in args.projects ]

    files = args.files
    if files:
        files = [ x.lower() for x in args.files ]
    main(args.output, args.folder, projects, args.download, files, args.since)
