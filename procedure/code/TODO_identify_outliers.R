## -------------------------------------------------------------------
## Script name: identify_outliers
## Purpose of script: A tool function to identify and optionally remove the 
## outliers in the cleaned data from function clean_data.

## Inputs:
## dat (data.frame): The data.frame to process.
## keep (logical): Keep the outliers in the input data or not. TRUE to keep,
## and FALSE to remove. The default is FALSE.

## Outputs:

## -------------------------------------------------------------------

# Because outlier detection depends on the modeling objective, it is easier
# to deal with the data.frame directly, rather than loading csv files.
identify_outliers <- function(dat,
                              keep = FALSE){
  # TODO: detect outliers here
  outliers <- NULL
  
  if (!keep){
    dat <- dat %>% filter(!station %in% outliers)
  }
  
  # Return: if keep the outliers in dat, then return outliers, otherwise return
  # the cleaned dat.
  if (keep){
    outliers
  } else{
    dat
  }
}
