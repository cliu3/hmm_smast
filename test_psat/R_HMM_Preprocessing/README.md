# R_HMM_Preprocessing
R code for preprocessing DST raw data for the Matlab HMM Geolocation toolbox [hmm_smast](https://github.com/cliu3/hmm_smast). 
Required imput data include meta tagging information from an Excel spreadsheet and the raw temperature and depth/pressure data downloaded from Star-ODDI or Lotek digital storage tags.

# Usage
1. Fill in meta tagging information in the Excel file, including the filenames of the raw data files (see `Inventory of yellowtail DST Recaptures.xlsx` for an example). Put raw data files in `raw_tags` folder.

2. Edit project-wide information in `metadata.R`.

3. In `process_tags.R`, specify the tags you desire to process by assigning corresponding list of tag IDs to the `ptags` variable.

4. Run `process_tags.R`. The processed tags will be stored in .mat format in `processed_tags` folder.

# Dependencies
* [gdata](https://cran.r-project.org/web/packages/gdata/index.html)
* [R.matlab](https://cran.r-project.org/web/packages/R.matlab/index.html)
