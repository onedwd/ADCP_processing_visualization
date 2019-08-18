Aug. 18, 2019:
This folder contains several versions of the adcpToolbox.R function code file from E.Chisholm's ADCP processing package.
The version currently a part of my fork of the package is adcpToolbox_P01.R, which follows the BODC SeaDataNet P01 variable naming conventions and uses variable names reflecting that the data have been corrected for magnetic declination.

adcpToolbox_GF3.R uses the GF3 variable naming scheme from https://www.nodc.noaa.gov/woce/woce_v3/wocedata_1/sss/documents/liste_param.htm. The variable names used reflect that the data have been corrected for magnetic declination.
The GF3 codes are no longer governed.

adcpToolbox_P01_L0.R uses the BODC SeaDataNet P01 variable naming scheme from http://seadatanet.maris2.nl/v_bodc_vocab_v2/search.asp?lib=P01. The variable names used reflect that the data were not corrected for magnetic declination. For low-level (level 0) processing only.
