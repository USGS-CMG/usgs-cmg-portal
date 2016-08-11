### USGS - Hurrican Sandy Rapid Response Stations

This is a python script to download and process NetCDF data from `http://ga.water.usgs.gov/flood/hurricane/sandy/datafiles/` and turn it into CF-1.6 compliant NetCDF DSG files.

### Installation

```bash
$ git clone https://github.com/axiom-data-science/usgs-cmg-portal.git
$ cd usgs-cmg-portal/usgs_rapid_release_gauges
```

##### Using `conda`

```
$ conda crate -n usgs_rapid python=3.5
$ source activate usgs_rapid
$ conda install -c conda-forge --file requirements.txt
```

### Running

##### See what the program supports
```bash
$ python collect.py --help
usage: collect.py [-h] [-d] [-o OUTPUT] [-f FOLDER] [-l [FILES [FILES ...]]]
                  {axiom,cf16}

positional arguments:
  {axiom,cf16}          Which type of file to produce. You most likely want
                        'cf16'.

optional arguments:
  -h, --help            show this help message and exit
  -d, --download        Should we download a new set of files or use the files
                        that have already been downloaded? Useful for
                        debugging. Downloaded files are never altered in any
                        way so you can rerun the processing over and over
                        without having to redownload any files.
  -o OUTPUT, --output OUTPUT
                        Directory to output NetCDF files to
  -f FOLDER, --folder FOLDER
                        Specify the folder location of ASCII files you wish to
                        translate. If this is used along with '--download',
                        the files will be downloaded into this folder and then
                        processed. If used without the '--download' option,
                        this is the location of the root folder you wish to
                        translate into NetCDF files.
  -l [FILES [FILES ...]], --files [FILES [FILES ...]]
                        Specific files to process (optional).
```

### Examples

##### Download and create CF-1.6 NetCDF
```bash
python collect.py cf16 --output=./output
```
