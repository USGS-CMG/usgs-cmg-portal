import logging
import logging.handlers

from thredds_crawler.crawl import Crawl

from sciwms_connect import SCIWMS_REST_URL, USER, PASSWORD
from sciwms_iso import add_dataset_to_sciwms, get_metadata


logger_name = 'thredds_crawler'
logger = logging.getLogger(logger_name)
fh = logging.handlers.RotatingFileHandler('/opt/docker/harvest/logs/iso_harvest.log', maxBytes=1024*1024*10, backupCount=5)
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
logger.addHandler(fh)
logger.addHandler(ch)
logger.setLevel(logging.DEBUG)

SAVE_DIR="/opt/docker/pycsw/force/iso_records"

THREDDS_SERVERS = {
    "necofs1": "http://www.smast.umassd.edu:8080/thredds/forecasts.html",
    "coawst":   "http://geoport-dev.whoi.edu/thredds/catalog/coawst_4/use/fmrc/catalog.html",
    "estofs_wh": "http://geoport-dev.whoi.edu/thredds/estofs_agg.html",
    "pacioos": "http://oos.soest.hawaii.edu/thredds/idd/ocn_mod.html"
}

metadata_files = get_metadata(thredds_servers=THREDDS_SERVERS,
                              save_dir=SAVE_DIR,
                              logger_name=logger_name
                              )
'''
add_dataset_to_sciwms(rest_url=SCIWMS_REST_URL,
                      user=USER,
                      password=PASSWORD,
                      metadata_files=metadata_files,
                      logger_name=logger_name
                      )
'''
