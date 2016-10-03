#!python
# coding=utf-8

from pyaxiom.netcdf.cf import CFDataset
from pyaxiom.netcdf.sensors.dsg import *

from scobs.utils import (
    normalize_time,
    normalize_epic_codes,
    normalize_vectors,
    normalize_units,
    check_compliance
)

import logging
logger = logging.getLogger(__name__)


def enhance_avg(file):
    # Normalize
    # BUG IN NETCDF:
    # https://github.com/Unidata/netcdf-c/commit/586e832f648600d4580f2b9d5e84036096cc3951
    # We need to wait for a release!
    # normalize_time(file)

    normalize_epic_codes(file)
    normalize_vectors(file)
    normalize_units(file)
    normalize_adcp_location(file)
    normalize_adcp_depths(file)

    # Now check for compliance!
    scores = check_compliance(file)
    for s in scores:
        level = 'info'
        if not s.passed:
            level = 'error'
        getattr(logger, level)(s)

    # Load as a CF DSG
    CFDataset.load(file)


def enhance_wvs(file):
    # Normalize
    normalize_time(file)
    normalize_epic_codes(file)
    normalize_vectors(file)
    normalize_units(file)
    #normalize_adcp_location(file)
    #normalize_adcp_depths(file)

    # Now check for compliance!
    scores = check_compliance(file)
    for s in scores:
        level = 'info'
        if not s.passed:
            level = 'error'
        getattr(logger, level)(s)

    # Load as a CF DSG
    CFDataset.load(file)


def normalize_adcp_location(netcdf_file):
    pass


def normalize_adcp_depths(netcdf_file):
    pass
