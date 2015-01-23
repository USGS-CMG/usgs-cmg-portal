### Woods Hold CMG Observation Data

This is a python script to download and process NetCDF data from `http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/catalog.html` and turn it into CF-1.6 compliant NetCDF DSG files.

### Installation

```bash
$ git clone https://github.com/axiom-data-science/usgs-cmg-portal.git
$ cd usgs-cmg-portal/woods_hole_obs_data
$ pip install -r requirements.txt
```


### CSV Metadata file

The script uses a CSV metadata file to get information about each project. There is a CSV file called `project_metadata.csv` that is included with the project.  The `project_name` and `catalog_xml` columns are **REQUIRED**.  Any other column names are added as global attributes to each file produced from the specifed `catalog_xml` link.

##### Example

Below, the global attributes `contributer_name`, `project_title`, and `project_summary` are added to the resulting NetCDF files.

```csv
project_name,contributor_name,project_title,project_summary,catalog_xml
"ARGO_MERCHANT","B. Butman","Argo Merchant Experiment","A moored array deployed after the ARGO MERCHANT ran aground onNantucket Shoals designed to help understand the fate of the spilled oil.","http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/ARGO_MERCHANT/catalog.xml"
"BARNEGAT","N. Ganju","Light attenuation and sediment resuspension in Barnegat Bay New Jersey"," Light attenuation is a critical parameter governing the ecological function of shallow estuaries.  Near-bottom and mid-water observations of currents, pressure, chlorophyll, and fDOM were collected at three pairs of sites sequentially at different locations in the estuary to characterize the conditions.","http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/BARNEGAT/catalog.xml"
"BUZZ_BAY","B. Butman","Currents and Sediment Transport in Buzzards Bay","Investigation of the near-bottom circulation in Buzzards Bay and consequent transport of fine-grained sediments that may be contaminated with PCBs from inner New Bedford Harbor.","http://geoport.whoi.edu/thredds/catalog/usgs/data2/emontgomery/stellwagen/Data/BUZZ_BAY/catalog.xml"
...
```


### Running

##### See what the program supports
```bash
$ python collect.py --help
usage: collect.py [-h] -o [OUTPUT] [-d] [-f FOLDER]
                  [-p [PROJECTS [PROJECTS ...]]] [-c [CSV_METADATA_FILE]]
                  {axiom,cf16}

positional arguments:
  {axiom,cf16}          Which type of file to produce. You most likely want
                        'cf16'.

optional arguments:
  -h, --help            show this help message and exit
  -o [OUTPUT], --output [OUTPUT]
                        Directory to output NetCDF files to
  -d, --download        Should we download a new set of files or use the files
                        that have already been downloaded? Useful for
                        debugging. Downloaded files are never altered in any
                        way so you can rerun the processing over and over
                        without having to redownload any files.
  -f FOLDER, --folder FOLDER
                        Specify the folder location of NetCDF files you wish
                        to translate. If this is used along with '--download',
                        the files will be downloaded into this folder and then
                        processed. If used without the '--download' option,
                        this is the location of the root folder you wish to
                        translate into NetCDF files.
  -p [PROJECTS [PROJECTS ...]], --projects [PROJECTS [PROJECTS ...]]
                        Specific projects to process (optional).
  -c [CSV_METADATA_FILE], --csv_metadata_file [CSV_METADATA_FILE]
                        CSV file to load metadata about each project from.
                        Defaults to 'project_metadata.csv'.
```


### Examples

##### Download and create CF-1.6 NetCDF files for a single project
```bash
python collect.py --download \
                  --projects SOUTHERN_CAL \
                  --output=./output/ \
                  cf16
```

##### Download and create NetCDF files for multiple projects
```bash
python collect.py --download \
                  --projects SOUTHERN_CAL DIAMONDSHOALS MYRTLEBEACH \
                  --output=./output/ \
                  cf16
```


##### Specify a different CSV metadata file
```bash
python collect.py --download \
                  --projects SOUTHERN_CAL DIAMONDSHOALS MYRTLEBEACH \
                  --output=./output/ \
                  --csv_metadata_file /some/path/to/your/file.csv \
                  cf16
```

##### Reprocess already downloaded files with new CSV metadata
```bash
python collect.py --projects SOUTHERN_CAL \
                  --output=./output/ \
                  cf16
```
