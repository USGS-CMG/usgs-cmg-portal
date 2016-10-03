#!python
# coding=utf-8
import os
import shutil
import unittest
import tempfile
from os.path import join as pjoin
from os.path import dirname as dname
from os.path import abspath as apath


from scobs import ctd, adcp

import logging
import coloredlogs
logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG')


class TestCtd(unittest.TestCase):

    def test_ctd(self):
        fpath = apath(pjoin(
            dname(__file__),
            'resources',
            'tripod',
            'north',
            'ctd',
            'MCR13N1T_CTD.nc'
        ))
        _, tmpf = tempfile.mkstemp(suffix='.nc', prefix='test')
        shutil.copy2(fpath, tmpf)
        ctd.enhance_ctd(tmpf)
        os.remove(tmpf)


class TestAdcp(unittest.TestCase):

    def test_adcp_avg60(self):
        fpath = apath(pjoin(
            dname(__file__),
            'resources',
            'tripod',
            'north',
            'adcp',
            'MCR13N1T_ADCP_AVG60.nc'
        ))
        _, tmpf = tempfile.mkstemp(suffix='.nc', prefix='test')
        shutil.copy2(fpath, tmpf)
        #adcp.enhance_avg(tmpf)
        os.remove(tmpf)

    def test_adcp_wvs(self):
        fpath = apath(pjoin(
            dname(__file__),
            'resources',
            'tripod',
            'north',
            'adcp',
            'MCR13N1T_wvs.nc'
        ))
        _, tmpf = tempfile.mkstemp(suffix='.nc', prefix='test')
        shutil.copy2(fpath, tmpf)
        #adcp.enhance_wvs(tmpf)
        os.remove(tmpf)
