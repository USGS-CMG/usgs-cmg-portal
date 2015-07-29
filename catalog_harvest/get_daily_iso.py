import os
import urllib
import logging
import logging.handlers

from thredds_crawler.crawl import Crawl

from iso_utils import IsoMetadata, construct_dataset_name_from_url
from sciwms_connect import SCIWMS_REST_URL, USER, PASSWORD
from sciwms_requests import SciWMSApi, determine_grid_type


logger = logging.getLogger('thredds_crawler')
fh = logging.handlers.RotatingFileHandler('/usgs/data0/iso/logs/iso_harvest.log', maxBytes=1024*1024*10, backupCount=5)
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
logger.addHandler(fh)
logger.addHandler(ch)
logger.setLevel(logging.DEBUG)

SAVE_DIR="/usgs/data0/iso/iso_records"

THREDDS_SERVERS = {
    "necofs1": "http://www.smast.umassd.edu:8080/thredds/forecasts.html",
    "coawst":   "http://geoport-dev.whoi.edu/thredds/catalog/coawst_4/use/fmrc/catalog.html" 
}

metadata_files = []
for subfolder, thredds_url in THREDDS_SERVERS.items():
    logger.info("Crawling %s (%s)" % (subfolder, thredds_url))
    crawler = Crawl(thredds_url, debug=True)
    isos = [(d.id, s.get("url")) for d in crawler.datasets for s in d.services if s.get("service").lower() == "iso"]
    filefolder = os.path.join(SAVE_DIR, subfolder)
    if not os.path.exists(filefolder):
        os.makedirs(filefolder)
    for iso in isos:
        filename = iso[0].replace("/", "_") + ".iso.xml"
        filepath = os.path.join(filefolder, filename)
        try:
            logger.info("Downloading/Saving %s" % filepath)
            urllib.urlretrieve(iso[1], filepath)
        except BaseException:
            logger.exception("Error!")
        else:
            metadata_files.append(filepath)           
swa = SciWMSApi(SCIWMS_REST_URL, USER, PASSWORD)
for metadata_file in metadata_files:
    logger.info("Reading ISO metadata file %s." % metadata_file)
    im = IsoMetadata(metadata_file)
    if im.is_cmg_portal_project():
        dataset_url = im.get_dataset_url()
        logger.info("Extracted this dataset URL from the metadata: %s." % dataset_url)
        dataset_name = construct_dataset_name_from_url(dataset_url)
        dataset_abstract = im.get_abstract()
        dataset_type = determine_grid_type(dataset_url)
        logger.info("Posting %s to %s." % (dataset_url, SCIWMS_REST_URL))
        post_dataset_resp = swa.create_dataset(uri=dataset_url,
                                               dataset_name=dataset_name,
                                               title=dataset_name,
                                               abstract=dataset_abstract,
                                               ds_type=dataset_type
                                               )
        logger.info("Posted %s to %s with a status code of %s." % (dataset_url,
                                                                   SCIWMS_REST_URL,
                                                                   post_dataset_resp.status_code
                                                                   )
                    )
        new_getcaps_url = swa.construct_getcapabilities_url(dataset_name)
        logger.info("Replacing the WMS GetCapabilities endpoint with: %s." % new_getcaps_url)
        im.replace_nciso_wms_getcaps_endpoints(new_getcaps_url)
    