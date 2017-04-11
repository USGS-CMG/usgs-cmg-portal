#!python
# coding=utf-8
import os
import sys
import pytest
import shutil
import unittest
from os.path import join as pjoin
from os.path import dirname as dname
from os.path import abspath as apath

from scobs import ctd, adcp, adv

import logging
import coloredlogs
logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG')


@pytest.mark.skipif(sys.version_info.major < 3,
                    reason="requires python3")
@pytest.mark.skip(reason="These tests are only experimental")
class TestLoad(unittest.TestCase):

    def setup_method(self, method):
        outdir = pjoin(dname(__file__), 'output', type(self).__name__)
        try:
            os.makedirs(outdir)
        except OSError:
            pass
        self.tmpf = pjoin(outdir, '{}.nc'.format(method.__name__))

    def test_ctd(self):
        fpath = apath(pjoin(
            dname(__file__),
            'resources',
            'tripod',
            'north',
            'ctd',
            'MCR13N1T_CTD.nc'
        ))
        shutil.copy2(fpath, self.tmpf)
        for of in ctd.enhance_ctd(self.tmpf):
            logger.info('Wrote output: {}'.format(of))

    def test_adcp_avg60(self):
        fpath = apath(pjoin(
            dname(__file__),
            'resources',
            'tripod',
            'north',
            'adcp',
            'MCR13N1T_ADCP_AVG60.nc'
        ))
        shutil.copy2(fpath, self.tmpf)
        for of in adcp.enhance_avg(self.tmpf):
            logger.info('Wrote output: {}'.format(of))

    def test_adcp_wvs(self):
        fpath = apath(pjoin(
            dname(__file__),
            'resources',
            'tripod',
            'north',
            'adcp',
            'MCR13N1T_wvs.nc'
        ))
        shutil.copy2(fpath, self.tmpf)
        for of in adcp.enhance_wvs(self.tmpf):
            logger.info('Wrote output: {}'.format(of))

    def test_adv_summary(self):
        fpath = apath(pjoin(
            dname(__file__),
            'resources',
            'tripod',
            'north',
            'adv',
            'MCR13N1T05adv1s-cal.nc'
        ))
        shutil.copy2(fpath, self.tmpf)
        for of in adv.enhance_summary(self.tmpf):
            logger.info('Wrote output: {}'.format(of))
