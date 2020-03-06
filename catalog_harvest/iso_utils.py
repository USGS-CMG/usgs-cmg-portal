'''
Created on Jul 28, 2015

@author: ayan
'''
from urllib.parse import urlparse

from lxml import etree


def construct_dataset_name_from_url(dataset_url, slice_int=-3):
    u = urlparse(dataset_url)
    u_path = u.path
    no_extn = u_path.split('.')[0].split('/')
    dataset_name = '_'.join(no_extn[slice_int:])
    return dataset_name


class IsoMetadata(object):

    namespaces = {'xsi': 'http://www.w3.org/2001/XMLSchema-instance',
                  'gco': 'http://www.isotc211.org/2005/gco',
                  'gmd': 'http://www.isotc211.org/2005/gmd',
                  'gmi': 'http://www.isotc211.org/2005/gmi',
                  'srv': 'http://www.isotc211.org/2005/srv',
                  'gmx': 'http://www.isotc211.org/2005/gmx',
                  'gsr': 'http://www.isotc211.org/2005/gsr',
                  'gss': 'http://www.isotc211.org/2005/gss',
                  'gts': 'http://www.isotc211.org/2005/gts',
                  'gml': 'http://www.opengis.net/gml/3.2',
                  'xlink': 'http://www.w3.org/1999/xlink',
                  'xs': 'http://www.w3.org/2001/XMLSchema'
                  }

    def __init__(self, metadata_path):
        self.md = metadata_path
        self.tree = etree.parse(self.md)

    def is_cmg_portal_project(self):
        xpath_param = (".//gmd:keyword/gco:CharacterString")
        elements = self.tree.xpath(xpath_param, namespaces=self.namespaces)
        check_results = []
        for element in elements:
            e_text = element.text.lower()
            contains_cmg_portal = 'cmg_portal' in e_text
            check_results.append(contains_cmg_portal)
        if sum(check_results) > 0:
            cmg_portal = True
        else:
            cmg_portal = False
        return cmg_portal

    def get_dataset_url(self):
        xpath_param = (".//srv:connectPoint"
                       "/gmd:CI_OnlineResource[gmd:name/gco:CharacterString='OPeNDAP']"
                       "/gmd:linkage/gmd:URL"
                       )
        elements = self.tree.xpath(xpath_param,
                                 namespaces=self.namespaces
                                 )
        dataset_url = elements[0].text
        return dataset_url

    def get_abstract(self):
        element = self.tree.find('.//gmd:abstract/gco:CharacterString', namespaces=self.namespaces)
        try:
            abstract_content = element.text
        except:
            abstract_content = None
        return abstract_content

    def replace_nciso_wms_getcaps_endpoints(self, new_endpoint):
        xpath_param = (".//srv:connectPoint"
                       "/gmd:CI_OnlineResource[gmd:name/gco:CharacterString='OGC-WMS']"
                       "/gmd:linkage/gmd:URL"
                       )
        elements = self.tree.xpath(xpath_param, namespaces=self.namespaces)
        for element in elements:
            element.text = new_endpoint
        modified_content = etree.tostring(self.tree, pretty_print=True)
        with open(self.md, 'w') as f:
            f.write(modified_content)
        return modified_content
