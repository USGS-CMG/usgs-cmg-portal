#!python
# coding=utf-8
import os
import sys
import shutil
import unittest
import tempfile
from os.path import join as pjoin
from os.path import dirname as dname
from os.path import abspath as apath

import logging
logger = logging.getLogger()
logger.level = logging.DEBUG
logger.handlers = [logging.StreamHandler()]

# Hack to get imports working with pytests
sys.path.append(dname(dname(apath(__file__))))
from collect import enhance_ctd, enhance_adcp


class TestToDataframe(unittest.TestCase):

    def test_13N1T_ADCP(self):
        fpath = apath(pjoin(dname(__file__), 'resources', 'RIVET2', 'MCR13N1T_ADCP_AVG60.nc'))
        _, tmpf = tempfile.mkstemp(suffix='.nc', prefix='test')
        shutil.copy2(fpath, tmpf)
        enhance_adcp(tmpf)
        os.remove(tmpf)
        
    def test_13N1T_CTD(self):
        fpath = apath(pjoin(dname(__file__), 'resources', 'RIVET2', 'MCR13N1T_CTD.nc'))
        _, tmpf = tempfile.mkstemp(suffix='.nc', prefix='test')
        shutil.copy2(fpath, tmpf)
        enhance_ctd(tmpf)
        os.remove(tmpf)
