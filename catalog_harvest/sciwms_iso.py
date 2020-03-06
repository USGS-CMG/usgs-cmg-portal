'''
Created on Jul 30, 2015

@author: ayan
'''
import os
import logging
try:
    from urllib import urlretrieve
except ImportError:
    from urllib.request import urlretrieve

from thredds_crawler.crawl import Crawl

from iso_utils import IsoMetadata, construct_dataset_name_from_url
from sciwms_requests import SciWMSApi, determine_grid_type


def get_metadata(thredds_servers, save_dir,
                 skips=Crawl.SKIPS, select=None,
                 debug=True, logger_name=None):
    logger = logging.getLogger(logger_name)
    tsi = thredds_servers.items()
    local_metadata_paths = []
    for subfolder, thredds_url in tsi:
        logger.info("Crawling {0} ({1})".format(subfolder, thredds_url))
        crawler = Crawl(thredds_url, skip=skips, select=select, debug=debug)
        filefolder = os.path.join(save_dir, subfolder)
        if not os.path.exists(filefolder):
            os.makedirs(filefolder)
        isos = [(d.id, s.get("url")) for d in crawler.datasets for s in d.services if s.get("service").lower() == "iso"]
        for iso in isos:
            filename = '{0}{1}'.format(iso[0].replace('/', '_'), '.iso.xml')
            filepath = os.path.join(filefolder, filename)
            try:
                urlretrieve(iso[1], filepath)
            except BaseException:
                logger.exception("Error!")
            else:
                local_metadata_paths.append(filepath)
    return local_metadata_paths


def add_dataset_to_sciwms(rest_url, user,
                          password, metadata_files, 
                          logger_name=None):
    logger = logging.getLogger(logger_name)
    swa = SciWMSApi(rest_url, user, password)
    for metadata_file in metadata_files:
        logger.info("Reading ISO metadata file {0}".format(metadata_file))
        im = IsoMetadata(metadata_file)
        if im.is_cmg_portal_project():
            dataset_url = im.get_dataset_url()
            dataset_name = construct_dataset_name_from_url(dataset_url)
            dataset_abstract = im.get_abstract()
            dataset_type = determine_grid_type(dataset_url)
            logger.info("Posting {0} to {1}".format(dataset_url, rest_url))
            post_dataset_resp = swa.create_dataset(uri=dataset_url,
                                                   dataset_name=dataset_name,
                                                   title=dataset_name,
                                                   abstract=dataset_abstract,
                                                   ds_type=dataset_type,
                                                   )
            logger.info("Posted {0} to {1} with a status code of {2}".format(dataset_url,
                                                                              rest_url,
                                                                              post_dataset_resp.status_code
                                                                              )
                        )
            new_getcaps_url = swa.construct_getcapabilities_url(dataset_name)
            logger.info("Replacing the WMS GetCapabilities endpoint with: {0}".format(new_getcaps_url))
            im.replace_nciso_wms_getcaps_endpoints(new_getcaps_url)
        else:
            logger.info("{0} does not contain the keyword for inclusion in the CMG Portal".format(metadata_file))
