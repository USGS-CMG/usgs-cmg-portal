#!python
# coding=utf-8

from utils import (
    normalize_time,
    normalize_epic_codes,
    normalize_vectors,
    normalize_units,
    normalize_ctd_location,
    normalize_ctd_depths,
    normalize_adcp_location,
    normalize_adcp_depths
)

import logging
logger = logging.getLogger()
logger.addHandler(logging.NullHandler())


def enhance_ctd(file):
    # Normalize
    normalize_time(file)
    normalize_epic_codes(file)
    normalize_vectors(file)
    normalize_units(file)
    normalize_ctd_location(file)
    normalize_ctd_depths(file)


def enhance_adcp(file):
    # Normalize
    # BUG IN NETCDF: https://github.com/Unidata/netcdf-c/commit/586e832f648600d4580f2b9d5e84036096cc3951
    # normalize_time(file)

    normalize_epic_codes(file)
    normalize_vectors(file)
    normalize_units(file)
    normalize_adcp_location(file)
    normalize_adcp_depths(file)
