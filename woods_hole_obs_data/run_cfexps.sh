#!/bin/bash
# script to control running collect.py on experiments, 1 by 1
#  make sure (ioos) prepends the command line label by:
# export PATH="/home/usgs/miniconda/bin:$PATH"
# source activate ioos
# cd to /usgs/data2/emontgomery/stellwagen/usgs-cmg-portal/woods_hole_obs_data
#
# comment out any that went with the runn all so you re-run just the ones you need to
#
python collect.py -p ARGO_MERCHANT --folder download -o ../../CF-1.6new
python collect.py -p BARNEGAT --folder download -o ../../CF-1.6new
python collect.py -p BARNEGAT_2014 --folder download -o ../../CF-1.6new
python collect.py -p BUZZ_BAY --folder download -o ../../CF-1.6new
python collect.py -p BW2011 --folder download -o ../../CF-1.6new
python collect.py -p CALCSTWET --folder download -o ../../CF-1.6new
python collect.py -p CAMP --folder download -o ../../CF-1.6new
python collect.py -p CAPE_COD_BAY --folder download -o ../../CF-1.6new
python collect.py -p CC_MISC --folder download -o ../../CF-1.6new
python collect.py -p CHANDELEUR_13 --folder download -o ../../CF-1.6new
python collect.py -p CHINCOTEAGUE --folder download -o ../../CF-1.6new
python collect.py -p DAUPHIN --folder download -o ../../CF-1.6new
python collect.py -p DEEP_REEF --folder download -o ../../CF-1.6new
python collect.py -p DIAMONDSHOALS --folder download -o ../../CF-1.6new
python collect.py -p DWDS_106 --folder download -o ../../CF-1.6new
python collect.py -p ECOHAB_I --folder download -o ../../CF-1.6new
python collect.py -p ECOHAB_II --folder download -o ../../CF-1.6new
python collect.py -p EUROSTRATAFORM --folder download -o ../../CF-1.6new
python collect.py -p FARALLONES --folder download -o ../../CF-1.6new
python collect.py -p FI12 --folder download -o ../../CF-1.6new
python collect.py -p FI14 --folder download -o ../../CF-1.6new
python collect.py -p GB_SED --folder download -o ../../CF-1.6new
python collect.py -p GLOBEC_GB --folder download -o ../../CF-1.6new
python collect.py -p GLOBEC_GSC --folder download -o ../../CF-1.6new
python collect.py -p GOMEX_CI --folder download -o ../../CF-1.6new
python collect.py -p GULF_MAINE --folder download -o ../../CF-1.6new
python collect.py -p HUDSON_SVALLEY --folder download -o ../../CF-1.6new
python collect.py -p HURRIRENE_bb --folder download -o ../../CF-1.6new
python collect.py -p KARIN_RIDGE --folder download -o ../../CF-1.6new
python collect.py -p LYDONIA_C --folder download -o ../../CF-1.6new
python collect.py -p MAB_SED --folder download -o ../../CF-1.6new
python collect.py -p MAMALA_BAY --folder download -o ../../CF-1.6new
python collect.py -p MBAY_CIRC --folder download -o ../../CF-1.6new
python collect.py -p MBAY_IWAVE --folder download -o ../../CF-1.6new
python collect.py -p MBAY_LT --folder download -o ../../CF-1.6new
python collect.py -p MBAY_LTB --folder download -o ../../CF-1.6new
python collect.py -p MBAY_STELL --folder download -o ../../CF-1.6new
python collect.py -p MBAY_WEST --folder download -o ../../CF-1.6new
python collect.py -p MOBILE_BAY --folder download -o ../../CF-1.6new
python collect.py -p MONTEREY_BAY --folder download -o ../../CF-1.6new
python collect.py -p MONTEREY_CAN --folder download -o ../../CF-1.6new
python collect.py -p MVCO_11 --folder download -o ../../CF-1.6new
python collect.py -p MVCO_14 --folder download -o ../../CF-1.6new
python collect.py -p MYRTLEBEACH --folder download -o ../../CF-1.6new
python collect.py -p NE_SLOPE --folder download -o ../../CF-1.6new
python collect.py -p NEARSHORE --folder download -o ../../CF-1.6new
python collect.py -p OCEANOG_C --folder download -o ../../CF-1.6new
python collect.py -p ORANGE_COUNTY --folder download -o ../../CF-1.6new
python collect.py -p PONCHARTRAIN --folder download -o ../../CF-1.6new
python collect.py -p PV_SHELF --folder download -o ../../CF-1.6new
python collect.py -p PV_SHELF04 --folder download -o ../../CF-1.6new
python collect.py -p PV_SHELF07 --folder download -o ../../CF-1.6new
python collect.py -p RCNWR --folder download -o ../../CF-1.6new
python collect.py -p SAB_SED --folder download -o ../../CF-1.6new
python collect.py -p SANDWICH2016 --folder download -o ../../CF-1.6new
python collect.py -p SOUTHERN_CAL --folder download -o ../../CF-1.6new
python collect.py -p STRESS --folder download -o ../../CF-1.6new
python collect.py -p WFAL --folder download -o ../../CF-1.6new
python collect.py -p WRIGHTSVILLE --folder download -o ../../CF-1.6new
