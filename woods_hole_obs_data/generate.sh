#!/bin/bash

RAW_FOLDER=${1:-./download}
CF_FOLDER=${2:-./cf}
DV_FOLDER=${3:-./dv}
DOWNLOAD=${4:-false}

DO_DOWNLOAD=""
if [ "$DOWNLOAD" = "true" ]; then
    DO_DOWNLOAD="-d"
fi
echo "#/bin/bash" > ./run_cfexps.sh
echo "#/bin/bash" > ./run_devices.sh

for PROJECT in $(awk -F '","' '{print $1}' project_metadata.csv| sed '1d;/#/d;s/"//' | sort); do
    echo "python collect.py $DO_DOWNLOAD -p $PROJECT -f $RAW_FOLDER -o $CF_FOLDER" >> ./run_cfexps.sh
    echo "python devices.py -p $PROJECT -f $CF_FOLDER -o $DV_FOLDER" >> ./run_devices.sh
done

echo "Raw Folder (first argument):    $RAW_FOLDER"
echo "CF Folder (second argument):    $CF_FOLDER"
echo "Device Folder (third argument): $DV_FOLDER"
echo ""
echo "Files ./run_cfexps.sh and ./run_devices.sh have been generated!"
echo ""

exit 0

# To make parallel files run these:"
# cat run_cfexps.sh | sed 1d | split -l $(echo "$(cat ./run_cfexps.sh | sed 1d | wc -l) / 4" | bc) - run_cf_
# cat run_devices.sh | sed 1d | split -l $(echo "$(cat ./run_devices.sh | sed 1d | wc -l) / 4" | bc) - run_dv_
