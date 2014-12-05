#!python
# coding=utf-8

import os
import sys
import shutil
import logging
import requests
import argparse
from copy import copy
from glob import glob
from datetime import datetime

import epic
import netCDF4
import numpy as np

from thredds_crawler.crawl import Crawl
from pytools.netcdf.sensors.create import create_timeseries_file

# Log to stdout
logger = logging.getLogger()
logger.setLevel(logging.INFO)
ch = logging.StreamHandler(sys.stdout)
fh = logging.FileHandler('usgs_cmg.log')
ch.setLevel(logging.INFO)
fh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - [%(levelname)s] %(message)s')
ch.setFormatter(formatter)
fh.setFormatter(formatter)
logger.addHandler(ch)
logger.addHandler(fh)

# Don't show the HTTP connection spam
requests_log = logging.getLogger("requests").setLevel(logging.WARNING)
crawler_log = logging.getLogger("thredds_crawler").setLevel(logging.INFO)

variable_name_overrides = {
    'w_1204min' : dict(epic_code=1204, overrides=dict(cell_methods='time: minimum')),
    'u_1205min' : dict(epic_code=1205, overrides=dict(cell_methods='time: minimum')),
    'v_1206min' : dict(epic_code=1206, overrides=dict(cell_methods='time: minimum')),
    'w_1204max' : dict(epic_code=1204, overrides=dict(cell_methods='time: maximum')),
    'u_1205max' : dict(epic_code=1205, overrides=dict(cell_methods='time: maximum')),
    'v_1206max' : dict(epic_code=1206, overrides=dict(cell_methods='time: maximum')),
    'WG_402'    : dict(epic_code=401),
    'Turb'      : dict(epic_code=980),
    'Press'     : dict(epic_code=1301),
    'vspd_1'    : dict(epic_code=300),
    'vdir_1'    : dict(epic_code=310),
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
}

long_name_overrides = {
    'salinity 2 q':                         dict(epic_code=40),
    'salinity 1':                           dict(epic_code=40),
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
    'dissolved oxygen saturation (mg/l)':   dict(epic_code=None, overrides=dict(standard_name='mass_concentration_of_oxygen_in_sea_water',
                                                                                original_units='mg/l',
                                                                                units='kg/m^3',
                                                                                convert=lambda x: x/1000.)),
    'raw aanderaa dissolved oxygen concentration (um/kg)': dict(epic_code=65),
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
}

project_metadata = {
    'ARGO_MERCHANT':  'B. Butman:Argo Merchant Experiment:A moored array deployed after the ARGO MERCHANT ran aground onNantucket Shoals designed to help understand the fate of the spilled oil.',
    'BARNEGAT':       'N. Ganju:Light attenuation and sediment resuspension in Barnegat Bay New Jersey: Light attenuation is a critical parameter governing the ecological function of shallow estuaries.  Near-bottom and mid-water observations of currents, pressure, chlorophyll, and fDOM were collected at three pairs of sites sequentially at different locations in the estuary to characterize the conditions.',
    'BUZZ_BAY':       'B. Butman:Currents and Sediment Transport in Buzzards Bay:Investigation of the near-bottom circulation in Buzzards Bay and consequent transport of fine-grained sediments that may be contaminated with PCBs from inner New Bedford Harbor.',
    'BW2011':         'N. Ganju: Blackwater 2011: Oceanographic and Water-Quality Measurements made at several sites in 2 watersheds in Blackwater National Wildlife Refuge.',
    'CAMP':           'B. Butman:California Area Monitoring Program (CAMP):A four-year multi-disciplinary field and laboratory study to investigate the sediment transport regime in the vicinity of production drilling rigs in the Santa Barbara Basin',
    'CAPE_COD_BAY':   'B. Butman:Currents and Sediment Transport in Cape Cod Bay:A pilot study to determine the effect of winter storms on sediment movement at two potential dredge spoil disposal areas.',
    'CC_MISC':        'B. Butman:Transport studies - Nauset Inlet:Part of a collaborative study of sediment movement in Nauset Inlet.',
    'CHANDELEUR':     'C. Sherwood:Chandeleur Islands Oceanographic Measurements:A program to measure waves water levels and currents near the Chandeleur Islands Louisiana and adjacent berm construction site.',
    'DEEP_REEF':      'J. Lacey:Gulf of Mexico - Pinnacles:Pressure data from the Gulf of Mexico',
    'DIAMONDSHOALS':  'J. Warner:Cape Hatteras- Diamond Shoals:This experiment was designed to investigate the ocean circulation and sediment transport dynamics at Diamond Shoals NC.',
    'DWDS_106':       'B. Butman:Sediment Transport at Deep Water Dump Site 106:Near-bottom current measurements to understand the fate and transport of sludge from the New York Metropolitan region discharged at the sea surface.',
    'ECOHAB_I':       'R. Signell:Ecology of Harmful Algal Blooms (ECOHAB-I):A field program to study the transport and fate of toxic dinoflagellate blooms in the western Gulf of Maine.',
    'ECOHAB_II':      'R. Signell:Ecology of Harmful Algal Blooms (ECOHAB-II):A field program to continue investigating the transport and fate of toxic dinoflagellate blooms in the western Gulf of Maine.',
    'EUROSTRATAFORM': 'C. Sherwood:EuroSTRATAFORM:The EuroSTRATAFORM Po and Apennine Sediment Transport and Accumulation (PASTA) experiment was an international study of sediment-transport processes and formation of geological strata in the Adriatic Sea.',
    'FARALLONES':     'M. Noble:Farallons:Program to measure the currents and circulation on the continental slope off San Francisco CA and thus infer the transport of dredged materialat the newly-established deep-water disposal site.',
    'FI12':           'J. Warner:Fire Island NY - Offshore: Oceanographic and meteorological observations were made at 7 sites on and around the sand ridges offshore of Fire Island NY in winter 2012 to study coastal processes.',
    'GB_SED':         'B. Butman:Georges Bank Current and Sediment Transport Studies:A series of studies to assess environmental hazards to petroleum development in the Georges Bank and New England Shelf region',
    'GLOBEC_GB':      'R. Schlitz:GLOBEC Georges Bank Program:A moored array program to investigate the circulation and mixing of plankton on Georges Bank.',
    'GLOBEC_GSC':     'R. Schlitz:GLOBEC Great South Channel Circulation Experiment:A moored array program to investigate the recirculation of water and plankton around Georges Bank.',
    'GULF_MAINE':     'B. Butman:Deep Circulation in the Gulf of Maine:A two-year field study to investigate the deep flow between the major basins in the Gulf of Maine and the effects on the distribution of suspended sediments.',
    'HUDSON_SVALLEY': 'B. Butman:Circulation and Sediment Transport in the Hudson Shelf Valley:Field experiments have been carried out to understand the transport of sediments and associated contaminants in the Hudson Shelf Valley offshore of New York.',
    'HURRIRENE_BB':   'B. Butman: Observations in Buzzards Bay during and after a Hurricane: Oceanographic data collected in Buzzards Bay MA during Hurricane Irene August 2011.',
    'KARIN_RIDGE':    'M. Noble:Karin Ridge Experiment:Current measurements collected at 2 sites in Karin Ridge Seamount.',
    'LYDONIA_C':      'B. Butman:Lydonia Canyon Dynamics Experiment:A major field experiment to determine the importance of submarine canyons in sediment transport along and across the continental margin.',
    'MAB_SED':        'B. Butman:Sediment Transport Observations in the Middle Atlantic Bight:A series of studies to assess environmental hazards to petroleum development in the Middle Atlantic Bight.',
    'MAMALA_BAY':     'D. Cacchione:Mamala bay Experiment:Current measurements collected at 350-450 meters in Mamala Bay near Waikiki Beach.',
    'MBAY_CIRC':      'R. Signell: Massachusetts Bay Circulation Experiment:Current measurements collected at 6 sites in Massachusetts Bay throughout the year to map the tidal wind and density driven currents.',
    'MBAY_IWAVE':     'B. Butman:Massachusetts Bay Internal Wave Experiment:A 1-month 4-element moored array experiment to measure the currents associated with large-amplitude internal waves generated by tidal flow across Stellwagen Bank.',
    'MBAY_LT':        'B. Butman:Long-term observations in Massachusetts Bay; Site A-Boston Harbor:Measurements of currents and other oceanographic properties were made to assess the impact of sewage discharge from the proposed outfall site.',
    'MBAY_LTB':       'B. Butman:Long-term observations in Massachusetts Bay; Site B-Scituate:Measurements of currents and other oceanographic properties were made to assess the impact of sewage discharge from the proposed outfall site.',
    'MBAY_STELL':     'R. Signell:Monitoring on Stellwagen Bank:A year-long series of current measurements on the eastern flank of Stellwagen Bank to document the currents at the mouth of Massachusetts Bay driven by the Maine Coastal current.',
    'MBAY_WEST':      'B. Butman:Currents and Sediment Transport in Western Massachusetts Bay:A pilot winter-time experiment to investigate circulation and sediment transport. Designed to provide information to aid in citing the new ocean outfall for the Boston sewer system.',
    'MOBILE_BAY':     'B. Butman:Mobile Bay Study:Measure currents and transport out of Mobile Bay.',
    'MONTEREY_BAY':   'M. Noble:Monterey Bay National Marine Sanctuary Program:Part of a large multi-disciplinary experiment to characterize the geologic environment and to generate a sediment budget.',
    'MONTEREY_CAN':   'M. Noble:Monterey Canyon Experiment: A program to determine the mechanisms that govern the circulation within and the transport of sediment and water through Monterey Submarine Canyon.',
    'MVCO_11':        'C. Sherwood: OASIS MVCO 2011: Near-seabed Oceanographic Observations made as part of the 2011 OASIS Project at the MVCO.',
    'MYRTLEBEACH':    'J. Warner:Myrtle Beach Experiment SC:Measurements collected as part of a larger study to understand the physical processes that control the transport of sediments in Long Bay South Carolina.',
    'NE_SLOPE':       'B. Butman:Currents on the New England Continental Slope:A study designed to describe the currents and to investigate the transport of sediment from the shelf to the slope.',
    'OCEANOG_C':      'B. Butman:Oceanographer Canyon Dynamics Experiment:A field experiment to determine the importance of submarine canyons in sediment transport along and across the continental margin.',
    'ORANGE_COUNTY':  'M. Noble:Orange County Sanitation District Studies:Observations to monitor coastal ocean process that transport suspended material and associated comtaminants across the shelf.',
    'PONCHARTRAIN':   'R. Signell:Lake Ponchartrain Project:A series of moored array studies to investigate the circulation and particle transport in Lake Pontchartrain.',
    'PV_SHELF':       'M. Noble:Palos Verdes Shelf Study:Initial observations of currents and circulation near the White Point ocean outfalls determine how often coastal ocean processes move the DDT contaminated sediments in this region.',
    'PV_SHELF04':     'M. Noble:Palos Verdes Shelf 2004:Additional observations to estimate the quantity and direction of sediment erosion and transport on the shelf near the White Point ocean outfalls.',
    'PV_SHELF07':     'M. Noble:Palos Verdes Shelf 2007:Follow-up observations to evaluate how often coastal ocean processes move the DDT contaminated sediments near the White Point ocean outfalls.',
    'SAB_SED':        'B. Butman:Sediment Transport Observations in the Southern Atlantic Bight:A series of studies to assess environmental hazards to petroleum development in the South Atlantic Bight.',
    'SOUTHERN_CAL':   'M. Noble:Southern California Project:A series of moorings were deployed to understand how coastal ocean processes that move sediments change with location on the shelf.',
    'STRESS':         'B. Butman:Sediment Transport on Shelves and Slopes (STRESS):Experiment on the California continental margin to investigate storm-driven sediment transport.',
    'WFAL':           'N. Ganju:West Falmouth Harbor Fluxes:Oceanographic and water-quality observations made at six locations in West Falmouth Harbor and Buzzards Bay.',
    'WRIGHTSVILLE':   'R. Thieler:Wrightsville Beach Study: Measurements of bottom currents and waves to investigate the flow field and sediment transport in a rippled scour depression offshore of Wrightsville Beach NC.',
}


def nc_close(nc):
    if nc is not None:
        try:
            nc.sync()
            nc.close()
        except RuntimeError:
            pass


def download(folder):

    full_catalog = 'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/catalog.xml'
    adcp_test = 'http://geoport.whoi.edu/thredds/catalog/usgs/data2/rsignell/data/adcp/catalog.html'

    catalogs = {
        'ARGO_MERCHANT':  'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/ARGO_MERCHANT/catalog.xml',
        'BARNEGAT':       'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/BARNEGAT/catalog.xml',
        'BUZZ_BAY':       'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/BUZZ_BAY/catalog.xml',
        'BW2011':         'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/BW2011/catalog.xml',
        'CAMP':           'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/CAMP/catalog.xml',
        'CAPE_COD_BAY':   'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/CAPE_COD_BAY/catalog.xml',
        'CC_MISC':        'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/CC_MISC/catalog.xml',
        'CHANDELEUR':     'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/CHANDELEUR/catalog.xml',
        'DEEP_REEF':      'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/DEEP_REEF/catalog.xml',
        'DIAMONDSHOALS':  'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/DIAMONDSHOALS/catalog.xml',
        'DWDS_106':       'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/DWDS_106/catalog.xml',
        'ECOHAB_I':       'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/ECOHAB_I/catalog.xml',
        'ECOHAB_II':      'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/ECOHAB_II/catalog.xml',
        'EUROSTRATAFORM': 'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/EUROSTRATAFORM/catalog.xml',
        'FARALLONES':     'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/FARALLONES/catalog.xml',
        'FI12':           'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/FI12/catalog.xml',
        'GB_SED':         'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/GB_SED/catalog.xml',
        'GLOBEC_GB':      'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/GLOBEC_GB/catalog.xml',
        'GLOBEC_GSC':     'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/GLOBEC_GSC/catalog.xml',
        'GULF_MAINE':     'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/GULF_MAINE/catalog.xml',
        'HUDSON_SVALLEY': 'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/HUDSON_SVALLEY/catalog.xml',
        'HURRIRENE_BB':   'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/HURRIRENE_BB/catalog.xml',
        'KARIN_RIDGE':    'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/KARIN_RIDGE/catalog.xml',
        'LYDONIA_C':      'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/LYDONIA_C/catalog.xml',
        'MAB_SED':        'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/MAB_SED/catalog.xml',
        'MAMALA_BAY':     'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/MAMALA_BAY/catalog.xml',
        'MBAY_CIRC':      'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/MBAY_CIRC/catalog.xml',
        'MBAY_IWAVE':     'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/MBAY_IWAVE/catalog.xml',
        'MBAY_LT':        'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/MBAY_LT/catalog.xml',
        'MBAY_LTB':       'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/MBAY_LTB/catalog.xml',
        'MBAY_STELL':     'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/MBAY_STELL/catalog.xml',
        'MBAY_WEST':      'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/MBAY_WEST/catalog.xml',
        'MOBILE_BAY':     'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/MOBILE_BAY/catalog.xml',
        'MONTEREY_BAY':   'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/MONTEREY_BAY/catalog.xml',
        'MONTEREY_CAN':   'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/MONTEREY_CAN/catalog.xml',
        'MVCO_11':        'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/MVCO_11/catalog.xml',
        'MYRTLEBEACH':    'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/MYRTLEBEACH/catalog.xml',
        'NE_SLOPE':       'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/NE_SLOPE/catalog.xml',
        'OCEANOG_C':      'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/OCEANOG_C/catalog.xml',
        'ORANGE_COUNTY':  'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/ORANGE_COUNTY/catalog.xml',
        'PONCHARTRAIN':   'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/PONCHARTRAIN/catalog.xml',
        'PV_SHELF':       'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/PV_SHELF/catalog.xml',
        'PV_SHELF04':     'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/PV_SHELF04/catalog.xml',
        'PV_SHELF07':     'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/PV_SHELF07/catalog.xml',
        'SAB_SED':        'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/SAB_SED/catalog.xml',
        'SOUTHERN_CAL':   'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/SOUTHERN_CAL/catalog.xml',
        'STRESS':         'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/STRESS/catalog.xml',
        'WFAL':           'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/WFAL/catalog.xml',
        'WRIGHTSVILLE':   'http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/WRIGHTSVILLE/catalog.xml',
    }

    # Use thredds_crawler to find DAP endpoints of the RAW data.
    total_datasets = []
    skips = Crawl.SKIPS + ['.*OTHER.*', '.*ancillary.*']

    try:
        for k, v in catalogs.items():
            datasets = Crawl(v, select=['.*-[A|a]+\..*'], skip=skips).datasets
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
    for d in total_datasets:
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
            logger.info("{!s} saved".format(d.name, folder))
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
                    logger.info("No EPIC code specified on {0}".format(v))
                else:
                    attribs = epic.mapping.get(int(nc_var.epic_code), None)
                    if attribs is not None and attribs["standard_name"] is not None:
                        # Convert data to CF units
                        nc_var[:] = attribs["convert"](nc_var[:])
                        # Set attributes
                        nc_var.standard_name = attribs['standard_name']
                        nc_var.long_name     = attribs['long_name']
                        nc_var.units         = attribs['cf_units']
                        if attribs['cell_methods'] is not None:
                            nc_var.cell_methods = attribs['cell_methods']
                    else:
                        logger.warning("Could not find CF mapping for EPIC code {!s}".format(nc_var.epic_code))
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

        if east is not None and north is not None:
            # We have vectors... create the speed and direction variables
            speed = np.sqrt(np.square(east[:]) + np.square(north[:]))
            direction = np.degrees(np.arctan2(north[:], east[:]))

            east_fill_value = east._FillValue if hasattr(east, '_FillValue') else np.nan
            spd = nc.createVariable('spd_300', 'f4', east.dimensions, fill_value=east_fill_value)
            spd.standard_name = 'sea_water_speed'
            spd.epic_code     = 300
            spd[:] = speed

            drc = nc.createVariable('dir_310', 'f4', east.dimensions, fill_value=east_fill_value)
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
        time_data = netCDF4.num2date((np.int64(nc.variables['time'][:])-2400001)*3600*24*1000 + nc.variables['time2'][:], units=millisecond_units)
        nc.renameVariable("time", "old_time")
        nc.sync()

        time = nc.createVariable('time', 'f8', ('time'))
        time.units          = epoch_units
        time.standard_name  = "time"
        time.long_name      = "time of measurement"
        time.calendar       = "gregorian"
        time[:] = netCDF4.date2num(time_data, units=epoch_units)
    except BaseException:
        logger.exception("Error")
        raise
    finally:
        nc_close(nc)


def main(output, do_download):

    # Download files
    download_folder = os.path.abspath(os.path.join(".", "download"))

    if do_download:
        try:
            downloaded_files = download(download_folder)
        except KeyboardInterrupt:
            downloaded_files = []
    else:
        downloaded_files = glob(os.path.join(download_folder, "*"))

    temp_folder = os.path.abspath(os.path.join(".", "temp"))
    shutil.rmtree(temp_folder, ignore_errors=True)
    os.makedirs(temp_folder)

    for down_file in downloaded_files:

        temp_file = os.path.join(temp_folder, os.path.basename(down_file))
        shutil.copy(down_file, temp_file)

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
                station_id = "{0}_{1}".format(project_name, nc.MOORING[0:3]).lower()
                station_name = "{0} ({1})".format(project_name, nc.MOORING[0:3])
            else:
                try:
                    # Mooring ID is the first three numbers of the file
                    station_id = int(os.path.basename(down_file)[0:3])
                    station_id = "{0}_mooring_{0}".format(project_name, station_id)
                    station_name = "{0} Mooring ({0})".format(project_name, station_id)
                except BaseException:
                    logger.error("Could not create a suitable station_id. Skipping {0}.".format(down_file))
                    continue

            latitude  = nc.variables.get("lat")[0]
            longitude = nc.variables.get("lon")[0]
            starting  = datetime.utcfromtimestamp(nc.variables.get("time")[0])
            ending    = datetime.utcfromtimestamp(nc.variables.get("time")[-1])

            station_urn = "urn:ioos:station:{0}:{1}".format('gov.usgs.cmgp', station_id).lower()
            logger.info("STATION: {0}".format(station_urn))

            data_variables  = list()
            other_variables = list()
            coord_vars      = ['time', 'time2', 'old_time', 'depth', 'lat', 'lon']
            for v in nc.variables:
                if v in coord_vars:
                    # Skip coordinate variables
                    continue
                nc_var = nc.variables.get(v)
                if hasattr(nc_var, "cf_role"):
                    other_variables.append(v)
                elif hasattr(nc_var, "standard_name") and nc_var.standard_name.lower() == "time":
                    # Skip time variables
                    pass
                elif hasattr(nc_var, "standard_name") and nc_var.standard_name.lower() in ["latitude", "longitude"]:
                    # Skip lat/lon variables
                    pass
                elif hasattr(nc_var, "standard_name") and nc_var.standard_name.lower() == "surface_altitude":
                    other_variables.append(v)
                elif hasattr(nc_var, "axis"):
                    other_variables.append(v)
                elif hasattr(nc_var, "epic_code") and nc_var.epic_code in epic.metadata_codes:
                    other_variables.append(v)
                elif v == "bindist":
                    other_variables.append(v)
                elif hasattr(nc_var, "standard_name") and nc_var.standard_name.lower() in ['northward_sea_water_velocity', 'eastward_sea_water_velocity']:
                    # We created speed/direction variables, so skip these
                    continue
                elif hasattr(nc_var, "standard_name"):
                    if hasattr(nc_var, "cell_methods"):
                        data_variables.append((v, nc_var.standard_name, nc_var.cell_methods, ))
                    else:
                        data_variables.append((v, nc_var.standard_name, None))
                else:
                    if hasattr(nc_var, 'long_name'):
                        logger.warning("Skipping {0}, no standard_name attribute.  'long_name' is {1}.".format(v, nc_var.long_name))
                    else:
                        logger.warning("Skipping {0}, no standard_name or long_name attribute.".format(v))
                    continue

            for dv, std, cm in data_variables:
                try:
                    logger.info("Exporting: {0}".format(dv))
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

                    file_global_attributes = copy(global_attributes)
                    file_global_attributes['id']               = station_id
                    file_global_attributes['title']            = station_name
                    if project_name in project_metadata:
                        pi, title, summary = project_metadata.get(project_name).split(':')
                        file_global_attributes['contributor_name'] = pi
                        file_global_attributes['project_title']    = title
                        file_global_attributes['project_summary']  = summary

                    ts = nc.variables.get("time")[:]
                    zs = nc.variables.get("depth")[:]

                    times     = np.ma.repeat(ts, zs.size)
                    verticals = np.ma.ravel(np.ma.repeat([zs], ts.size, axis=0))
                    values    = nc.variables.get(dv)[:]

                    assert values.size == verticals.size == times.size
                    create_timeseries_file(output_directory=output_directory, latitude=latitude, longitude=longitude, full_station_urn=station_urn, full_sensor_urn=sensor_urn, global_attributes=file_global_attributes, attributes=nc.variables.get(dv).__dict__, output_filename=file_name, times=times, verticals=verticals, values=values)

                    new_nc = netCDF4.Dataset(os.path.join(output_directory, file_name), 'a')
                    for other in other_variables:
                        old_var = nc.variables.get(other)

                        # Get new variable name
                        variable_name = sensor_urn.split(":")[-1]
                        new_var = new_nc.variables.get(variable_name)

                        """
                        # Switch to using the 'height' dimension
                        dims = [ d if d != 'depth' else 'height' for d in old_var.dimensions ]
                        # No dimensions for
                        dims = filter(None, [ d if d not in ['lat', 'lon'] else None for d in dims ])
                        if 'height' in dims and 'height' not in nc.dimensions:
                            if len(nc.dimensions.get('depth')) == 1:
                                # Remove the depth dimension to match the CF file spec
                                dims = filter(None, [ d if d != 'height' else None for d in dims ])
                            else:
                                logger.info("Skipping: {0}.  It has a Z axis but the core variable '{1}' does not.".format(other, dv))
                                continue
                        """

                        logger.info("Adding: {0}".format(other))
                        other_var = new_nc.createVariable(other, old_var.dtype, new_var.dimensions)
                        for k in old_var.ncattrs():
                            other_var.setncattr(k, old_var.getncattr(k))
                        try:
                            other_var[:] = nc.variables.get(other)[:]
                        except BaseException:
                            logger.info("Skipping: {0}.  It has more dimensions than the core variable '{1}'.".format(other, dv))
                            continue
                        new_nc.sync()
                    nc_close(new_nc)

                except BaseException:
                    logger.exception("Error. Skipping {0} in {1}.".format(down_file, dv))
                    if os.path.exists(os.path.join(output_directory, file_name)):
                        os.remove(os.path.join(output_directory, file_name))
                    continue

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
                        help="Should we download the files or use the temp files?  Useful for debugging.")
    args = parser.parse_args()
    main(args.output, args.download)
