import logging
import logging.handlers

#from sciwms_connect import SCIWMS_REST_URL, USER, PASSWORD
from sciwms_iso import add_dataset_to_sciwms, get_metadata


logger_name = 'thredds_crawler'
logger = logging.getLogger('thredds_crawler')
fh = logging.handlers.RotatingFileHandler('/opt/docker/harvest/logs/iso_harvest_ts.log', maxBytes=1024*1024*10, backupCount=5)
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
logger.addHandler(fh)
logger.addHandler(ch)
logger.setLevel(logging.DEBUG)

SAVE_DIR="/opt/docker/pycsw/store/iso_records"

THREDDS_SERVERS = {
   "usgs-ts": "https://geoport.usgs.esipfed.org/thredds/catalog/silt/usgs/Projects/stellwagen/CF-1.6/catalog.html",
   "F7VD6XBF": "https://geoport.usgs.esipfed.org/thredds/catalog/sand/usgs/users/dnowacki/doi-F7VD6XBF/catalog.html", # Ocean Currents and Pressure Time Series at the Upper Florida Keys: Crocker Reef, FL
   "F7RN362H": "https://geoport.usgs.esipfed.org/thredds/catalog/sand/usgs/users/dnowacki/doi-F7RN362H/catalog.html", # Data from coastal circulation and water-column properties in the National Park of American Samoa, February-July 2015
   "F73R0R07": "https://geoport.usgs.esipfed.org/thredds/catalog/sand/usgs/users/dnowacki/doi-F73R0R07/catalog.html", # Wind-wave and suspended-sediment data from Liberty Island and Little Holland Tract, Sacramento-San Joaquin Delta, California (ver. 2.0, September 2019)
   "F7CR5RW8": "https://geoport.usgs.esipfed.org/thredds/catalog/sand/usgs/users/dnowacki/doi-F7CR5RW8/catalog.html", # Oceanographic measurements obtained offshore of the Elwha River delta in coordination with the Elwha River Restoration Project, Washington, USA, 2010-2014
   "F7NG4NS1": "https://geoport.usgs.esipfed.org/thredds/catalog/sand/usgs/users/dnowacki/doi-F7NG4NS1/catalog.html", # Oceanographic measurements and hydrodynamic modeling of the mouth of the Columbia River, Oregon and Washington, 2013
   "F71C1W36": "https://geoport.usgs.esipfed.org/thredds/catalog/sand/usgs/users/dnowacki/doi-F71C1W36/catalog.html", # Bathymetry, acoustic-backscatter, and time-series datasets collected between 2014 and 2016 of a field of crescent-shaped rippled scour depressions in northern Monterey Bay, California
   "P9IXOHID": "https://geoport.usgs.esipfed.org/thredds/catalog/sand/usgs/users/dnowacki/doi-P9IXOHID/catalog.html", # Time series data of oceanographic conditions from La Parguera, Puerto Rico, 2017-2018 Coral Reef Circulation and Sediment Dynamics Experiment
   "P91T185R": "https://geoport.usgs.esipfed.org/thredds/catalog/sand/usgs/users/dnowacki/doi-P91T185R/catalog.html", # Time Series of Autonomous Carbonate System Parameter Measurements in Eastern Gulf of Mexico near Tampa Bay, Florida, USA
   "F7FT8J7Q": "https://geoport.usgs.esipfed.org/thredds/catalog/sand/usgs/users/dnowacki/doi-F7FT8J7Q/catalog.html", # Time-series oceanographic data from the Monterey Canyon, CA October 2015 - March 2017
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
