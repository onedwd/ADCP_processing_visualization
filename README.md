# ADCP_processing_visualization
This repository contains files associated with using the `ADCP` package (https://github.com/hhourston/mooredDataProcessing_adcp). Usage of the `ADCP` package is directed here.

Sample files for using the `ADCP` package: \
ADCP metadata template: \
R processing script: \

Sample file for creating N/E and along- and cross-shore current velocity plots from ADCP files in netCDF format: \
Python plotting script: \

Sample files for using `ADCP` from Python: \
Python script for calling R: \
Modified version of *ADCP_lvl1_process.R* that is callable from Python: \

Sample ADCP files: \
Raw ADCP file (not named according to conventions): \
Filled-out csv metadata file for the above raw ADCP file: \
Output netCDF file from the above raw ADCP file, produced using *ADCP_lvl1_process.R*: \

Other versions of adcpToolbox.R besides adcpToolbox_P01.R that is currently in the forked repository mooredDataProcessing_adcp.

## Installation
It is not necessary that this whole repository be cloned by a user. It is recommended that files in this repository be downloaded individually as needed.

Steps to download an individual file from a GitHub repository:
1. Open a file in a GitHub repository and click on the button "raw" in the top right corner to view the raw file in a new browser tab. Copy the url of the raw file (e.g. https://raw.githubusercontent.com/username/reponame/path/to/file).
2. To save this file to your computer, open a terminal window and type
    
        wget https://raw.githubusercontent.com/username/reponame/path/to/file
    
### Dependencies of ADCP_lvl1_process.R and ADCP_lvl1_process_pycalled.R
To install dependencies:

    install.packages("tools")
    install.packages("gsw")
    install.packages("english")
