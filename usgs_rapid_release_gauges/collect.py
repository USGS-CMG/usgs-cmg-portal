#!python
# coding=utf-8

import os
import logging
import argparse
import calendar
from datetime import datetime

import pytz
import requests
from bs4 import BeautifulSoup

from pyaxiom.netcdf.sensors import TimeSeries

from dateutil.parser import parse

import numpy as np
import pandas as pd

import coloredlogs

requests_log = logging.getLogger("requests.packages.urllib3")
requests_log.setLevel(logging.WARN)

# Log to stdout
logger = logging.getLogger()
coloredlogs.install(level=logging.INFO)
logger.setLevel(logging.INFO)


def main(output_format, output):

    waf = "http://ga.water.usgs.gov/flood/hurricane/sandy/datafiles/"

    r = requests.get(waf)
    soup = BeautifulSoup(r.text)
    for link in soup.find_all('a'):

        # Skip non .txt files
        site_id, ext = os.path.splitext(link['href'])
        if ext != ".txt":
            continue

        csv_link = waf + link['href']
        logger.info("Downloading '{}'".format(csv_link))
        d = requests.get(csv_link)
        try:
            d.raise_for_status()
        except requests.exceptions.HTTPError:
            logger.error("Could not download: {!s}, skipping. Status code: {!s}".format(csv_link, d.status_code))
            continue

        contents = d.text
        for line in contents.split("\n"):
            if "agency_cd" in line:
                parse_type_1(output_format, site_id, contents, output, csv_link)
                break
            elif "date_time_GMT" in line:
                parse_type_2(output_format, site_id, contents, output, csv_link)
                break
            else:
                continue


def split_file(c, data_key):
    comments = list()
    data     = list()
    headers  = list()

    lines    = c.split("\n")
    for line in lines:
        if line.startswith("#") or not line:
            comments.append(line)
        elif data_key in line:
            headers = map(lambda x: x.strip(), line.split("\t"))
            continue
        else:
            row = map(lambda x: x.strip(), line.split("\t"))
            if filter(None, row):
                assert len(row) == len(headers)
                data.append(row)

    return comments, headers, data


def parse_type_1(output_format, site_id, contents, output, csv_link):
    """
    # ---------------------------------- WARNING ----------------------------------------
    # The data you have obtained from this automated U.S. Geological Survey database
    # have not received Director's approval and as such are provisional and subject to
    # revision.  The data are released on the condition that neither the USGS nor the
    # United States Government may be held liable for any damages resulting from its use.
    # Additional info: http://waterdata.usgs.gov/ga/nwis/help/?provisional
    #
    # File-format description:  http://waterdata.usgs.gov/nwis/?tab_delimited_format_info
    # Automated-retrieval info: http://waterdata.usgs.gov/nwis/?automated_retrieval_info
    #
    # Contact:   gs-w_support_nwisweb@usgs.gov
    # retrieved: 2012-11-20 12:05:22 EST       (caww01)
    #
    # Data for the following 1 site(s) are contained in this file
    #    USGS 395740074482628 South Branch Rancocas Cr at S Main St nr Lumberton
    # -----------------------------------------------------------------------------------
    #
    # Data provided for site 395740074482628
    #    DD parameter   Description
    #    03   00035     Wind speed, miles per hour
    #    07   00025     Barometric pressure, millimeters of mercury
    #    09   00045     Precipitation, total, inches
    #    19   63160     Stream water level elevation above NAVD 1988, in feet
    #
    # Data-value qualification codes included in this output:
    #     P  Provisional data subject to revision.
    #
    agency_cd   site_no datetime    tz_cd   03_00035    03_00035_cd 07_00025    07_00025_cd 09_00045    09_00045_cd 19_63160    19_63160_cd
    5s  15s 20d 6s  14n 10s 14n 10s 14n 10s 14n 10s
    USGS    395740074482628 2012-10-28 13:00    EST 4.2 P   755 P           3.22    P
    USGS    395740074482628 2012-10-28 13:15    EST 6.4 P   754 P   0.00    P   3.36    P
    USGS    395740074482628 2012-10-28 13:30    EST 3.6 P   754 P   0.00    P   3.50    P
    USGS    395740074482628 2012-10-28 13:45    EST 3.2 P   754 P   0.00    P   3.63    P
    USGS    395740074482628 2012-10-28 14:00    EST 7.0 P   754 P   0.00    P   3.76    P
    USGS    395740074482628 2012-10-28 14:15    EST 4.0 P   754 P   0.00    P   3.87    P
    ...
    """
    # lat/lon point: http://waterservices.usgs.gov/nwis/site/?sites=395740074482628

    variable_map = {
        '01_00065' : {'long_name' : 'Gage height', 'geoid_name' : 'NAVD88', 'vertical_datum' : 'NAVD88', 'water_surface_reference_datum' : 'NAVD88', 'standard_name' : 'water_surface_height_above_reference_datum', 'units': 'feet'},
        '03_00035' : {'long_name' : 'Wind Speed', 'standard_name' : 'wind_speed', 'units': 'mph'},
        '04_00035' : {'long_name' : 'Wind Gust', 'standard_name' : 'wind_speed_of_gust', 'units': 'mph'},
        '05_00035' : {'long_name' : 'Wind Speed', 'standard_name' : 'wind_speed', 'units': 'mph'},
        '06_00035' : {'long_name' : 'Wind Gust', 'standard_name' : 'wind_speed_of_gust', 'units': 'mph'},
        '04_00036' : {'long_name' : 'Wind Direction', 'standard_name' : 'wind_from_direction', 'units': 'degrees'},
        '02_00036' : {'long_name' : 'Wind Direction', 'standard_name' : 'wind_from_direction', 'units': 'degrees'},
        '05_00025' : {'long_name' : 'Air Pressure', 'standard_name' : 'air_pressure', 'units': 'mm of mercury'},
        '07_00025' : {'long_name' : 'Air Pressure', 'standard_name' : 'air_pressure', 'units': 'mm of mercury'},
        '09_00025' : {'long_name' : 'Air Pressure', 'standard_name' : 'air_pressure', 'units': 'mm of mercury'},
        '03_00045' : {'long_name' : 'Total Precipitation', 'standard_name' : 'lwe_thickness_of_precipitation_amount', 'units': 'inches'},
        '08_00045' : {'long_name' : 'Total Precipitation', 'standard_name' : 'lwe_thickness_of_precipitation_amount', 'units': 'inches'},
        '09_00045' : {'long_name' : 'Total Precipitation', 'standard_name' : 'lwe_thickness_of_precipitation_amount', 'units': 'inches'},
        '06_00052' : {'long_name' : 'Relative Humidity', 'standard_name' : 'relative_humidity', 'units': 'percent'},
        '07_00052' : {'long_name' : 'Relative Humidity', 'standard_name' : 'relative_humidity', 'units': 'percent'},
        '08_00052' : {'long_name' : 'Relative Humidity', 'standard_name' : 'relative_humidity', 'units': 'percent'},
        '05_00020' : {'long_name' : 'Air Temperature', 'standard_name' : 'air_temperature', 'units': 'degrees_Celsius'},
        '06_00020' : {'long_name' : 'Air Temperature', 'standard_name' : 'air_temperature', 'units': 'degrees_Celsius'},
        '07_00020' : {'long_name' : 'Air Temperature', 'standard_name' : 'air_temperature', 'units': 'degrees_Celsius'},
        '19_63160' : {'long_name' : 'Water Surface Height Above Reference Datum (NAVD88)', 'geoid_name' : 'NAVD88', 'vertical_datum' : 'NAVD88', 'water_surface_reference_datum' : 'NAVD88', 'standard_name' : 'water_surface_height_above_reference_datum', 'units': 'feet'},
        '01_63160' : {'long_name' : 'Water Surface Height Above Reference Datum (NAVD88)', 'geoid_name' : 'NAVD88', 'vertical_datum' : 'NAVD88', 'water_surface_reference_datum' : 'NAVD88', 'standard_name' : 'water_surface_height_above_reference_datum', 'units': 'feet'},
    }

    # Get metadata from a seperate endpoint.
    d = requests.get("http://waterservices.usgs.gov/nwis/site/?sites={!s}".format(site_id))
    try:
        d.raise_for_status()
    except requests.exceptions.HTTPError:
        logger.error("Could not find lat/lon endpoint for station {!s}, skipping. Status code: {!s}".format(site_id, d.status_code))
        return
    _, hz, dz = split_file(d.text, "agency_cd")
    # Strip off the one line after the headers
    dz = dz[1:]
    dfz  = pd.DataFrame(dz, columns=hz)
    lat  = float(dfz["dec_lat_va"][0])
    lon  = float(dfz["dec_long_va"][0])
    sensor_vertical_datum = dfz["alt_datum_cd"][0] or "NAVD88"
    try:
        z = float(dfz["alt_va"][0])
    except ValueError:
        z = 0.
    loc  = "POINT({!s} {!s} {!s})".format(lon, lat, z)
    name = dfz["station_nm"][0]

    comments, headers, data = split_file(contents, "agency_cd")
    df = pd.DataFrame(data, columns=headers)

    fillvalue = -9999.9

    # Combine date columns
    dates = df["datetime"]
    tz = df["tz_cd"]
    new_dates = list()
    for i in range(len(dates)):
        try:
            new_dates.append(parse(dates[i] + " " + tz[i]).astimezone(pytz.utc))
        except BaseException:
            # Remove row.  Bad date.
            df.drop(i, axis=0, inplace=True)
            continue
    df['time'] = new_dates
    df['depth'] = [ z for x in range(len(df['time'])) ]

    # Strip out "_cd" columns (quality checks for USGS)
    for h in headers:
        if "_cd" in h:
            df.drop(h, axis=1, inplace=True)

    # Add global attributes to appear in the resulting NetCDF file
    global_attributes = dict(
        title=name,
        summary='USGS Hurricane Sandy Rapid Response Stations.  Data acquired from "http://ga.water.usgs.gov/flood/hurricane/sandy/datafiles/.',
        keywords="usgs, waterdata, elevation, water, waterlevel, sandy, hurricane, rapid, response, %s" % site_id,
        keywords_vocaublary="None",
        naming_authority='gov.usgs',
        id=site_id,
        cdm_data_type="Station",
        history="NetCDF file generated from {!s}".format(csv_link),
        creator="USGS",
        creator_url="http://waterdata.usgs.gov",
        creator_institution="USGS",
        creator_urn="gov.usgs",
        publisher="Axiom Data Science",
        publisher_uri="http://axiomdatascience.com",
        processing_level="None",
        acknowledgement="None",
        geospatial_bounds=loc,
        geospatial_lat_min=lat,
        geospatial_lat_max=lat,
        geospatial_lon_min=lon,
        geospatial_lon_max=lon,
        license="Freely Distributed",
        date_created=datetime.utcnow().replace(second=0, microsecond=0).isoformat()
    )

    def to_floats(x):
        try:
            return float(x)
        except ValueError:
            return fillvalue

    min_time = df['time'].min()
    max_time = df['time'].max()

    full_station_urn = "urn:ioos:station:{!s}:{!s}".format(global_attributes["naming_authority"], site_id)
    if output_format == 'cf16':
        output_filename = '{}_{}-{}.nc'.format(site_id, min_time.strftime('%Y%m%dT%H%M%S'), max_time.strftime('%Y%m%dT%H%M%S'))
        times = [ calendar.timegm(x.timetuple()) for x in df["time"] ]
        verticals = df['depth'].values
        ts = TimeSeries(output, latitude=lat, longitude=lon, station_name=full_station_urn, global_attributes=global_attributes, output_filename=output_filename, times=times, verticals=verticals, vertical_axis_name='z')

    for var in df.columns:
        if var in ['datetime', 'time', 'depth', 'tz_cd', 'site_no', 'agency_cd']:
            continue

        try:
            var_meta = variable_map[var]
        except KeyError:
            logger.error("Variable {!s} was not found in variable map!".format(var))
            continue

        # Convert to floats
        df[var] = df[var].map(to_floats)

        # Change feet to meters
        if var_meta["units"] in ["feet", "ft"]:
            df[var] = np.asarray([ v * 0.3048 if v != fillvalue else v for v in df[var] ])
            var_meta["units"] = "meters"

        if output_format == 'axiom':
            full_sensor_urn = "urn:ioos:sensor:{!s}:{!s}:{!s}".format(global_attributes["naming_authority"], site_id, var_meta["standard_name"])
            output_directory = os.path.join(output, full_sensor_urn)
            output_filename = '{}_{}-{}.nc'.format(var, min_time.strftime('%Y%m%dT%H%M%S'), max_time.strftime('%Y%m%dT%H%M%S'))
            ts = TimeSeries.from_dataframe(df, output_directory, output_filename, lat, lon, full_station_urn, global_attributes, var_meta["standard_name"], var_meta, sensor_vertical_datum=sensor_vertical_datum, fillvalue=fillvalue, data_column=var, vertical_axis_name='height')
            ts.add_instrument_metadata(urn=full_sensor_urn)
            ts.close()
        elif output_format == 'cf16':
            # Variable names shouldn't start with a number
            try:
                int(var[0])
                variable_name = 'v_{}'.format(var)
            except:
                variable_name = var
            ts.add_variable(variable_name, values=df[var].values, attributes=var_meta, fillvalue=fillvalue, sensor_vertical_datum=sensor_vertical_datum)

    if output_format == 'cf16':
        ts.close()


def parse_type_2(output_format, site_id, contents, output, csv_link):
    """
    # These data are provisional and subject to revision.
    # Data processed as of 12/05/2012 11:54:29.
    # Data collected as part of Hurricane Sandy (2012) Storm Tide project.
    # Data are archived at http://water.usgs.gov/floods/events/2012/isaac/index.php
    # Elevation determined from GPS surveys (NAVD 88).
    # Time datum is GMT (Greenwich Mean Time).
    # Water density estimated on basis of sensor location
    #   where saltwater = 63.989 lb/ft3       (Saltwater = dissolved solids concentration greater than 20000 milligrams per liter)
    #   where brackish water = 63.052 lb/ft3  (Brackish water = dissolved solids concentration between 1000 and 20000 milligrams per liter)
    #   where freshwater = 62.428 lb/ft3      (Freshwater = dissolved solids concentration less than 1000 milligrams per liter)
    # The equation used to compute elevation from recorded pressure is
    #  (((sp-bp)*144)/d)+e
    # Where sp = surge pressure in psi; bp = barometric pressure in psi;
    #  d = water density in lb/ft3; and e = elevation of sensor in ft above NAVD 88.
    # Barometric data from nearest pressure sensor. Location for the barometric sensor is listed below.
    # Elevation is computer-rounded to two decimal places.
    #      Sensor information
    # Site id = SSS-NY-WES-001WL
    # Site type = water level
    # Horizontal datum used is NAD 83
    # Sensor location latitude 40.942755
    # Sensor location longitude -73.719828
    # Sensor elevation above NAVD 88 = -3.97 ft
    # Lowest recordable water elevation is -3.90 ft
    # Water density value used = 63.989 lb/ft3
    # Barometric sensor site (source of bp) = SSS-NY-WES-002BP
    # Barometric sensor location latitude 40.90754368
    # Barometric sensor location longitude -73.8692184

    date_time_GMT   elevation   nearest_barometric_sensor_psi
    10-28-2012 06:00:00 0.88    14.5145
    10-28-2012 06:00:30 0.86    14.5145
    10-28-2012 06:01:00 0.85    14.5170
    10-28-2012 06:01:30 0.85    14.5145
    10-28-2012 06:02:00 0.84    14.5170
    10-28-2012 06:02:30 0.81    14.5145
    10-28-2012 06:03:00 0.76    14.5145
    ...
    """

    variable_map = {
        'elevation' : {'long_name' : 'Water Level Elevation above Reference Datum (NAVD88)', 'geoid_name' : 'NAVD88', 'vertical_datum' : 'NAVD88', 'water_surface_reference_datum' : 'NAVD88', 'standard_name' : 'water_surface_height_above_reference_datum', 'units': 'feet'},
    }

    def to_floats(x):
        try:
            return float(x)
        except ValueError:
            return fillvalue

    comments, headers, data = split_file(contents, "date_time_GMT")
    df = pd.DataFrame(data, columns=headers)
    fillvalue = -9999.9

    lat     = None
    lon     = None
    z       = 0
    name    = site_id
    sensor_vertical_datum = "NAVD88"

    for c in comments:
        if "Sensor location latitude" in c:
            lat = float(filter(None, map(lambda x: x.strip(), c.split(" ")))[-1])
        elif "Sensor location longitude" in c:
            lon = float(filter(None, map(lambda x: x.strip(), c.split(" ")))[-1])
        elif "Site id" in c:
            site_id = filter(None, map(lambda x: x.strip(), c.split(" ")))[-1]
            name = site_id
        elif "Sensor elevation" in c:
            sensor_vertical_datum = "".join(c.split("=")[0].split(" ")[4:6])
            l = filter(None, map(lambda x: x.strip(), c.split(" ")))
            z = float(l[-2])
            if l[-1] in ["feet", "ft"]:
                z *= 0.3048

    loc = "POINT({!s} {!s} {!s})".format(lon, lat, z)
    df['time'] = df["date_time_GMT"].map(lambda x: parse(x + " UTC"))
    df['depth'] = [ z for x in range(len(df['time'])) ]

    # Add global attributes to appear in the resulting NetCDF file
    global_attributes = dict(
        title=name,
        summary='USGS Hurricane Sandy Rapid Response Stations.  Data acquired from http://ga.water.usgs.gov/flood/hurricane/sandy/datafiles/.',
        keywords="usgs, waterdata, elevation, water, waterlevel, sandy, hurricane, rapid, response, %s" % site_id,
        keywords_vocaublary="None",
        naming_authority='gov.usgs',
        id=site_id,
        cdm_data_type="Station",
        history="NetCDF file generated from {!s}".format(csv_link),
        creator="USGS",
        creator_url="http://waterdata.usgs.gov",
        creator_institution="USGS",
        creator_urn="gov.usgs",
        publisher="Axiom Data Science",
        publisher_uri="http://axiomdatascience.com",
        processing_level="None",
        acknowledgement="None",
        geospatial_bounds=loc,
        geospatial_lat_min=lat,
        geospatial_lat_max=lat,
        geospatial_lon_min=lon,
        geospatial_lon_max=lon,
        license="Freely Distributed",
        date_created=datetime.utcnow().replace(second=0, microsecond=0).isoformat()
    )

    full_station_urn = "urn:ioos:station:{!s}:{!s}".format(global_attributes["naming_authority"], site_id)
    min_time = df["time"].min()
    max_time = df["time"].max()

    if output_format == 'cf16':
        times = [ calendar.timegm(x.timetuple()) for x in df['time'] ]
        verticals = df['depth'].values
        output_filename = '{}_{}-{}.nc'.format(site_id, min_time.strftime('%Y%m%dT%H%M%S'), max_time.strftime('%Y%m%dT%H%M%S'))
        ts = TimeSeries(output, latitude=lat, longitude=lon, station_name=full_station_urn, global_attributes=global_attributes, output_filename=output_filename, times=times, verticals=verticals)

    for var in df.columns:
        if var in ['date_time_GMT', 'time', 'depth']:
            continue
        try:
            int(var[0])
            variable_name = 'v_{}'.format(var)
        except:
            variable_name = var

        try:
            var_meta = variable_map[var]
        except KeyError:
            logger.error("Variable {!s} was not found in variable map!".format(var))
            continue

        # Convert to floats
        df[var] = df[var].map(to_floats)
        if var_meta["units"] in ["feet", "ft"]:
            df[var] = [ v * 0.3048 if v != fillvalue else v for v in df[var] ]
            var_meta["units"] = "meters"

        if output_format == 'axiom':
            full_sensor_urn = "urn:ioos:sensor:{!s}:{!s}:{!s}".format(global_attributes["naming_authority"], site_id, var_meta["standard_name"])
            output_directory = os.path.join(output, full_sensor_urn)
            output_filename = '{}_{}-{}.nc'.format(var, min_time.strftime('%Y%m%dT%H%M%S'), max_time.strftime('%Y%m%dT%H%M%S'))
            ts = TimeSeries.from_dataframe(df, output_directory, output_filename, lat, lon, full_station_urn, global_attributes, var_meta["standard_name"], var_meta, sensor_vertical_datum=sensor_vertical_datum, fillvalue=fillvalue, data_column=var)
            ts.add_instrument_metadata(urn=full_sensor_urn)
            ts.close()
        elif output_format == 'cf16':
            ts.add_variable(variable_name, values=df[var].values, attributes=var_meta, fillvalue=fillvalue, sensor_vertical_datum=sensor_vertical_datum)

    if output_format == 'cf16':
        ts.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('format', action='store', choices=('axiom', 'cf16',), help="Which type of file to produce. You most likely want 'cf16'.")
    parser.add_argument('-o', '--output',
                        required=True,
                        help="Directory to output NetCDF files to",
                        nargs='?')
    args          = parser.parse_args()
    main(args.format, args.output)
