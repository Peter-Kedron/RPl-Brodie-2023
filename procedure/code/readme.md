# Code
Store computational code-based research procedures here.
Document an index of files stored here in [procedures_index.csv](../procedure_index.csv) the root [procedure](../) folder.

Before running the code, please download the shared data from figshare and Asia and Pacific protected area database from WDPA website, and put them in `data/public`.

# Structure

- kick_off.R: The function to load libraries and functions for any usage.
- clean_pa.R: The function to query and clean Protected Areas (PA), then cast into single polygons, and finally cluster the polygons into different groups for mammals.
- pre_dm.R: The function to prepare the data model to calculate connectivity metrics of flux and area weighted flux based on four assumptions.
- calc_conn.R: The wrap function to load inputs and calculate connectivity metrics.
- clean_data.R: The function to prepare training data for the following models.
- identify_outliers.R: The function to detect outliers in the prepared data.
- model_pa_efficacy.R: The function to run models for PA efficacy.
- model_pa_spillover: The function to run models for testing PA spillover.
- compile_models.R: The wrap function to run all needed models.
- archive: holds old scripts.
- raw: holds the original scripts shared with the paper.
