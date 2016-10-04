# Santa Cruz Obs -> NetCDF

This is a proof of concept repository for discovering a better way to convert
USGS CMG observational data into CF compliant netCDF files.

This is not fully functional. It contains methods to convert (4) specific types
of files.

* Tripod - ADCP AVG
* Tripod - ADCP WVS
* Tripod - ADV
* Tripod - CTD


### Getting Started

Create a new conda environemnt

```python
conda create -n usgscmg35 python=3.5
```

Install dependencies

```python
conda install -c conda-forge -c axiom-data-science --file requirements.txt
```

Run the tests

```
py.test -s -rxs -v
```
