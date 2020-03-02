# ADCP_processing_visualization
This repository contains files associated with using the `ADCP` package (https://github.com/hhourston/mooredDataProcessing_adcp). Usage of the `ADCP` package is directed here.

#### Sample files for using the `ADCP` package
ADCP metadata template: *ADCP_metadata_template/ADCP_metadata_template_L1.csv*

R processing script: *ADCP_process_lvl1.R*

#### Plotting from netCDF files
Python script for creating N/E and along- and cross-shore current velocity plots from a netCDF-format ADCP file: *plot_westcoast_nc_P01.py*

#### Sample files for using `ADCP` from Python: in the subfolder *callR_fromPython*
Python script for calling R: *callR_ADCPprocessing.py*

Modified version of *ADCP_process_lvl1.R* that is callable from Python: *ADCP_process_lvl1_pycalled.R*

#### Sample ADCP files: in the subfolder *sample_files*
Raw ADCP file (not named according to conventions): *A1_UU_8745.000*

Filled-out csv metadata file for the above raw ADCP file: *a1_20080430_20080918_0489m_meta_L1.csv*

Output netCDF file from the above raw ADCP file, produced using *ADCP_lvl1_process.R*: *a1_20080430_20080918_0489m.adcp.L1.csv*

#### adcpToolbox.R versions
Other versions of adcpToolbox.R besides adcpToolbox_P01.R that is currently in the forked repository mooredDataProcessing_adcp can be found in **adcpToolbox_versions**. These versions are deprecated and should not be used.

#### Opening netCDF files
Scripts for opening and closing netCDF files in R and Python can be found in **open_nc**.

## Installation
It is not necessary that this whole repository be cloned by a user. It is recommended that files in this repository be downloaded individually as needed.

Steps to download an individual file from a GitHub repository:
1. Open a file in a GitHub repository and click on the button "raw" in the top right corner to view the raw file in a new browser tab. Copy the url of the raw file (e.g. https://raw.githubusercontent.com/username/reponame/path/to/file).
2. To save this file to your computer, open a terminal window and type
    
        wget https://raw.githubusercontent.com/username/reponame/path/to/file
    
### Dependencies of ADCP_lvl1_process.R and ADCP_lvl1_process_pycalled.R
To install dependencies in R:

    install.packages("tools")
    install.packages("gsw")
    install.packages("english")
        
