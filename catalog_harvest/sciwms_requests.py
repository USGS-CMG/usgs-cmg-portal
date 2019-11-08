'''
Created on Jul 28, 2015

@author: ayan
'''
import json
try:
    from urllib import urlencode
    from urlparse import urljoin
except ImportError:
    from urllib.parse import urlencode, urljoin

try:
    from pysgrid import from_ncfile
except ImportError:
    from pysgrid import read_netcdf as from_ncfile
# from pysgrid.custom_exceptions import SGridNonCompliantError
from pyugrid import UGrid
import requests


def determine_grid_type(dataset_url):
    try:
        UGrid.from_ncfile(dataset_url)
        grid_type = 'ugrid'
    except ValueError:
        try:
            from_ncfile(dataset_url)
            grid_type = 'sgrid'
        except SGridNonCompliantError:
            grid_type = None
    return grid_type
            

class SciWMSApi(object):
    
    def __init__(self, api_url, username, password, headers=None):
        self.api_url = api_url
        self.username = username
        self.password = password
        if headers is not None:
            self.headers = headers
        else:
            self.headers = {'Content-type': 'application/json'}
        # create a requests session
        self.rs = requests.Session()
        self.rs.auth = (self.username, self.password)
        self.rs.headers = self.headers
        
    def _construct_url(self, relative_url):
        if self.api_url[-1] != '/':
            safe_url = '{0}/'.format(self.api_url)
        else:
            safe_url = self.api_url
        target_url = urljoin(safe_url, relative_url)
        return target_url
    
    def get_dataset(self, dataset_id=None):
        if dataset_id is None:
            get_url = self._construct_url('datasets')
        else:
            rel_url = 'datasets/{0}/'.format(dataset_id)
            get_url = self._construct_url(rel_url)
        resp = self.rs.get(get_url)
        return resp.json()
    
    def create_dataset(self, uri, dataset_name,
                       title, abstract, ds_type,
                       keep_up_to_date=True):
        raw = {'uri': uri,
               'name': dataset_name,
               'title': title,
               'type': ds_type,
               'abstract': abstract,
               'keep_up_to_date': keep_up_to_date
               }
        if ds_type is not None:
            raw['type'] = ds_type
        payload = json.dumps(raw)
        post_url = self._construct_url('datasets')
        resp = self.rs.post(post_url, data=payload)
        return resp
    
    def sciwms_wms_endpoint(self, dataset_name):
        hostname = self.api_url.replace('rest/', '').replace('rest', '')
        wms_endpoint = '{0}wms/datasets/{1}'.format(hostname, dataset_name)
        return wms_endpoint
    
    def construct_getcapabilities_url(self, dataset_name):
        wms_endpoint = self.sciwms_wms_endpoint(dataset_name)
        params = {'REQUEST': 'GetCapabilities'}
        urlencoded_params = urlencode(params)
        getcaps_url = '{0}?{1}'.format(wms_endpoint, urlencoded_params)
        return getcaps_url
        
    def get_dataset_layer(self, layer_id):
        rel_url = 'layers/{0}/'.format(layer_id)
        get_url = self._construct_url(rel_url)
        resp = self.rs.get(get_url)
        return resp
    
    def modify_dataset_layer(self, layer_id, styles,
                             default_min, default_max,
                             logscale, description,
                             std_name, active, var_name,
                             units
                             ):
        raw = {'styles': styles,
               'default_min': default_min,
               'default_max': default_max,
               'logscale': logscale,
               'description': description,
               'std_name': std_name,
               'active': active,
               'var_name': var_name,
               'units': units
               }
        rel_url = 'layers/{0}/'.format(layer_id)
        put_url = self._construct_url(rel_url)
        payload = json.dumps(raw)
        resp = self.rs.put(put_url, data=payload)
        return resp
    
    def get_dataset_virtual_layer(self, vlayer_id):
        rel_url = 'vlayers/{0}/'.format(vlayer_id)
        get_url = self._construct_url(rel_url)
        resp = self.rs.get(get_url)
        return resp
    
    def modify_dataset_virtual_layer(self, vlayer_id, styles,
                                     default_min, default_max,
                                     logscale, description,
                                     std_name, active, var_name,
                                     units):
        raw = {'styles': styles,
               'default_min': default_min,
               'default_max': default_max,
               'logscale': logscale,
               'description': description,
               'std_name': std_name,
               'active': active,
               'var_name': var_name,
               'units': units
               }
        rel_url = 'vlayers/{0}/'.format(vlayer_id)
        put_url = self._construct_url(rel_url)
        payload = json.dumps(raw)
        resp = self.rs.put(put_url, data=payload)
        return resp
