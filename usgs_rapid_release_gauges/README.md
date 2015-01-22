### USGS - Hurrican Sandy Rapid Response Stations

This is a python script to download and process NetCDF data from `http://ga.water.usgs.gov/flood/hurricane/sandy/datafiles/` and turn it into CF-1.6 compliant NetCDF DSG files.

### Installation

```bash
$ git clone https://github.com/axiomalaska/usgs-cmg-portal.git
$ cd usgs-cmg-portal/usgs_rapid_release_gauges
$ pip install -r requirements.txt
```

### Running

##### See what the program supports
```bash
$ python collect.py --help
usage: collect.py [-h] -o [OUTPUT] {axiom,cf16}

positional arguments:
  {axiom,cf16}          Which type of file to produce. You most likely want
                        'cf16'.

optional arguments:
  -h, --help            show this help message and exit
  -o [OUTPUT], --output [OUTPUT]
                        Directory to output NetCDF files to
```

### Examples

##### Download and create CF-1.6 NetCDF
```bash
python collect.py --output=./output cf16
```
