## -------------------------------------------------------------------
## Script name: identify_outliers
## Purpose of script: A tool function to identify and optionally remove the 
## outliers in the cleaned data from function clean_data.

## Inputs:
## dat (data.frame): The data.frame of the processed data to use. dat must have
## two columns: y for independent variable and station for station index. It
## also should include multiple feature columns to define outliers.
## keep (logical): Keep the outliers in the input data or not. TRUE to keep,
## and FALSE to remove. The default is FALSE.

## Outputs:

## -------------------------------------------------------------------

# Because outlier detection depends on the modeling objective, it is easier
# to deal with the data.frame directly, rather than loading csv files.
identify_outliers <- function(dat){
    # Check structure
    if (!all(c("y", "station") %in% names(dat))){
        stop("dat must have column y and station to process.")
    }
    
  # TODO: Detect outliers here
  outliers <- NULL
  
  # Return
  outliers
}
