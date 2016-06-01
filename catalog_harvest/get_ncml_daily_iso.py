import logging
import logging.handlers

from thredds_crawler.crawl import Crawl

from sciwms_connect import SCIWMS_REST_URL, USER, PASSWORD
from sciwms_iso import add_dataset_to_sciwms, get_metadata


SAVE_DIR = "/usgs/data0/iso/iso_records"

THREDDS_SERVERS = {
'mvco': 'http://geoport-dev.whoi.edu/thredds/catalog/usgs/data0/mvco_ce/catalog.html',
'bbleh': 'http://geoport-dev.whoi.edu/thredds/catalog/usgs/data0/bbleh/catalog.html',
'bbleh-sandy': 'http://geoport-dev.whoi.edu/thredds/catalog/clay/usgs/users/zdefne/catalog.html',
'sandy-jcwarner': 'http://geoport-dev.whoi.edu/thredds/catalog/clay/usgs/users/jcwarner/Projects/Sandy/catalog.html',
'chinco-ab': 'http://geoport-dev.whoi.edu/thredds/catalog/clay/usgs/users/abeudin/chinco/catalog.html',
'chinco-nkg': 'http://geoport-dev.whoi.edu/thredds/catalog/clay/usgs/users/nganju/chincoteague_bedelevation/catalog.html'
'hudson-sandy': 'http://clancy.whoi.edu:8080/thredds/catalog/data1/dralston/hudson/sandy/catalog.html'
}


logger_name = 'thredds_crawler'
logger = logging.getLogger(logger_name)
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

skips = Crawl.SKIPS + ['.*MATLAB.*', '\..*']
select = ['.*00_dir.*.ncml']

metadata_files = get_metadata(thredds_servers=THREDDS_SERVERS,
                              save_dir=SAVE_DIR,
                              skips=skips, select=select,
                              logger_name=None
                              )
add_dataset_to_sciwms(rest_url=SCIWMS_REST_URL,
                      user=USER,
                      password=PASSWORD,
                      metadata_files=metadata_files,
                      logger_name=None
                      )
