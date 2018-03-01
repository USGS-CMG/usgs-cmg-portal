import logging
import logging.handlers

from thredds_crawler.crawl import Crawl

from sciwms_connect import SCIWMS_REST_URL, USER, PASSWORD
from sciwms_iso import add_dataset_to_sciwms, get_metadata


SAVE_DIR = "/opt/docker/pycsw/force/iso_records"
 
THREDDS_SERVERS = {
'zdefne-portal': 'https://gamone.whoi.edu/thredds/catalog/sand/usgs/users/zdefne/catalog.html',
'jcwarner-sandy-portal': 'https://gamone.whoi.edu/thredds/catalog/q',
'nganju-portal': 'https://gamone.whoi.edu/thredds/catalog/sand/usgs/users/nganju/portal_runs/catalog.html',
'projects': 'https://gamone.whoi.edu/thredds/catalog/sand/usgs/Projects/catalog.html'
}


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

skips = Crawl.SKIPS + ['.*MATLAB.*', '\..*']
select = ['.*00_dir.*.ncml']

metadata_files = get_metadata(thredds_servers=THREDDS_SERVERS,
                              save_dir=SAVE_DIR,
                              skips=skips, select=select,
                              logger_name=None
                              )
'''
add_dataset_to_sciwms(rest_url=SCIWMS_REST_URL,
                      user=USER,
                      password=PASSWORD,
                      metadata_files=metadata_files,
                      logger_name=None
                      )
'''
