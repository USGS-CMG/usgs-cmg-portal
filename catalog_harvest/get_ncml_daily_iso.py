import os
import urllib
from thredds_crawler.crawl import Crawl

import logging
import logging.handlers
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
'mvco': 'http://geoport-dev.whoi.edu/thredds/catalog/usgs/data0/mvco_ce/catalog.html',
'bbleh': 'http://geoport-dev.whoi.edu/thredds/catalog/usgs/data0/bbleh/catalog.html'
}

for subfolder, thredds_url in THREDDS_SERVERS.items():
    logger.info("Crawling %s (%s)" % (subfolder, thredds_url))
    skips = Crawl.SKIPS + ['.*MATLAB.*', '\..*']
    crawler = Crawl(thredds_url, debug=True, select=['.*00_dir.*.ncml'], skip=skips)
    isos = [(d.id, s.get("url")) for d in crawler.datasets for s in d.services if s.get("service").lower() == "iso"]
    filefolder = os.path.join(SAVE_DIR, subfolder)
    if not os.path.exists(filefolder):
        os.makedirs(filefolder)
    for iso in isos:
        try:
            filename = iso[0].replace("/", "_") + ".iso.xml"
            filepath = os.path.join(filefolder, filename)
            logger.info("Downloading/Saving %s" % filepath)
            urllib.urlretrieve(iso[1], filepath)
        except BaseException:
            logger.exception("Error!")
