## -------------------------------------------------------------------
## Script name: clean_data
## Purpose of script: Clean the training data: join original data with 
## landscape meatures, scale the relevant variables, and subset the data
## for relevant variables only.
## Author: Lei Song, Peter Kedron
## Date Created: 2024-08-13
## Email: lsong@ucsb.edu

## Import from package:
## dplyr

## Inputs:
## taxon (character): The taxon name. Either mammal or bird.
## conn_metrics (character): The connectivity measure to use. The value is 
## one of ["flux_ptg", "awf_ptg", "flux_gtg", "awf_gtg", "flux_rst_ptp", 
## "awf_rst_ptp", "flux_rst_ptp2", "awf_rst_ptp2"]. The default is "awf_ptg".
## See details in function prep_dm.
## src_dir (character): The directory of data source to read original CSVs.
## conn_dir (character): The directory to read CSVs of connectivity measures.
## Check function calc_conn for correct setting.
## dst_dir (character): The directory to save files to.

## Outputs:
## The data.frame of the cleaned data. A CSV file named "dat_analysis_*" is 
## saved out as well.
## Variables are: station, study_area (only mammal), country, PA, utm_east, 
## utm_north, dist_to_PA, PA_size_km2, asymptPD, maxFRic, SR.mean, med_dist, 
## utm_east.z, utm_north.z, HDI.z, access_log10.z, PA_size_km2.z, dist_to_PA.z,
## forest_structure, connectivity.z, BigPA, CloseToPA

## Usage example:
# taxon <- "bird"
# conn_metrics <- 'awf_ptg'
# src_dir <- "data/raw/public"
# conn_dir <- "data/derived/public"
# dst_dir <- "data/derived/public"
# dat_clean <- clean_data(taxon, conn_metrics, src_dir, conn_dir, dst_dir)
## -------------------------------------------------------------------

clean_data <- function(taxon, 
                       conn_metrics = "awf_ptg", 
                       src_dir, 
                       conn_dir, 
                       dst_dir){
  # Check inputs
  if (!taxon %in% c("bird", "mammal")){
    stop("Taxon must be either bird or mammal.")
  }
  if (!conn_metrics %in% c("flux_ptg", "awf_ptg", "flux_gtg", "awf_gtg", 
                           "flux_rst_ptp", "awf_rst_ptp", "flux_rst_ptp2", 
                           "awf_rst_ptp2")){
    stop("Wrong conn_metrics setting. Check function prep_dm for correct one.")
  }
  
  # Load the data provided by Brodie et al. (2023)
  dat_brodie <- data.frame(read.csv(
    file.path(src_dir, sprintf("training/%s_data_corrected_240122.csv", taxon)), 
    header = T))
  
  # Create data.frame containing the subset of variable used in the analysis
  ## Difference between bird and mammal is that using country and study_area 
  ## to match for mammals, but only use country for bird.
  nms <- c('station', 'country', 'PA', 'utm_east', 'utm_north', 
          'Hansen_recentloss', 'access_log10', 'HDI', 'dist_to_PA', 
          'PA_size_km2', 'rh_95_a0.pred', 'asymptPD', 'maxFRic', 'SR.mean')
  if(taxon == "mammal"){
      nms <- append(nms, "study_area", after = 1)} else {nms <- nms}
  
  dat <- dat_brodie %>% select(all_of(nms))
  
  # Add the connectivity variables for each station calculated with 
  # calc_conn_metrics.R
  dat_conn_metrics <- data.frame(read.csv(
    file.path(conn_dir, sprintf("conn_flux_%s_10_150.csv", taxon)), 
    header = T)) %>% 
    select(all_of(c("station", "med_dist", conn_metrics)))
  dat <- left_join(dat, dat_conn_metrics, by = "station")
  
  # Scale subset of continuous variables in dat
  dat_scale <- subset(
    dat, select = c("utm_east", "utm_north", 
                    "HDI", "access_log10", 
                    "PA_size_km2", 
                    "dist_to_PA", 
                    "rh_95_a0.pred", 
                    conn_metrics)) %>% 
    scale(. , center = TRUE, scale = TRUE) %>% data.frame()
  
  # Append scaled variables to data with .z suffixes
  dat[paste0(names(dat_scale), '.z')] <- dat_scale
  
  # Rename scaled, predicted relative canopy height at 95% (rh_95_a0.pred.z) 
  # as forest_structure to match DAG 
  names(dat)[names(dat) == "rh_95_a0.pred.z"] <- "forest_structure"
  
  # Exclude stations that underwent recent forest loss as defined by 
  # Hansen et al. (2013)
  dat_clean <- subset(dat, Hansen_recentloss == 0)
  
  # Subset data.frame to only contain relevant variables
  dat_clean <- dat_clean %>% 
    select(-all_of(c("Hansen_recentloss", "access_log10", "HDI", 
              "rh_95_a0.pred", conn_metrics)))
  
  # Rename the selected connectivity measure to a common name
  conn_nm <- sprintf("%s.z", conn_metrics)
  names(dat_clean)[names(dat_clean) == conn_nm] <- "connectivity.z" # match others
  
  # Add spillover variables
  ## PA size
  dat_clean$BigPA <- NA
  dat_clean[dat_clean$PA == 0, "BigPA"] <- 
      ifelse(dat_clean[dat_clean$PA == 0, "PA_size_km2"] < 500, 0, 1)
  
  ## Distance to PA
  dat_clean$CloseToPA <- NA
  dat_clean[dat_clean$PA == 0, "CloseToPA"] <- 
      ifelse(dat_clean[dat_clean$PA == 0, "dist_to_PA"] > 2, 0, 1)
  
  # Save out
  fname <- file.path(dst_dir, sprintf("dat_analysis_%s.csv", taxon))
  write.csv(dat_clean, fname, row.names = FALSE)
  
  # Return the data.frame
  return(dat_clean)
}
