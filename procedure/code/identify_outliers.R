## -----------------------------------------------------------------------------
## identify_outliers
## -----------------------------------------------------------------------------
#
# Purpose of script: -----------------------------------------------------------
# Use hatvalues to detect outliers in the data.

## Author: Wenxin Yang, Lei Song
# Date Created: 2024-08-22
# Last Update:  2024-09-01

## Import from package: nlme

# Inputs: ---------------------------------------------------------------------
## dat (data.frame): The data.frame of the processed data to use. dat must be
##                   the same data.frame to create mod.
## mod (lme model object): The model object to diagnose.
## leverage_threshold (numeric): The threshold in hatvalues to detect outliers.
##                               The default is NULL. Then the function will 
##                               use quantile(hatvalues, 0.99) as the threshold.

## Outputs:
## outliers (vector): A vector of station ids for outliers.

## WARNING:
## 1. Only using hatvalues is not a reliable way to detect outliers.
## The influence of each data point is quantified by both leverage and discrepancy.
## For instance, car package extend the generic function `influence` to diagnose
## linear mixed models. But it is extremely slow which I have no idea why.
## 2. Using hatvalues typically needs manual checking by plotting the points.
## Even though it is recommended to use 3 * mean(hatvalues) as the cut-line.
## -------------------------------------------------------------------

identify_outliers <- function(dat, mod,
                              leverage_threshold = NULL){
    # Calculate hatvalues
    lev_values <- hat(model.matrix(mod))
    
    # Select outliers
    dat$lev_value <- lev_values
    
    # # Calculate DFFITS
    # ei <- residuals(mod)
    # s_e <- summary(mod)$sigma
    # estar <- ei / (s_e * sqrt(1 - lev_values))
    # dffits <- estar * sqrt(lev_values / (1 - lev_values))
    # cutoff <- 2 * sqrt(length(mod$coefficients$fixed) / nrow(dat))
    
    # Calculate leverage_threshold
    leverage_threshold <- ifelse(
        is.null(leverage_threshold), 
        quantile(lev_values, 0.99), 
        leverage_threshold)
    
    # Get the station ids for outliers
    outliers <- dat[dat$lev_value >= leverage_threshold, "station"]
    
    # Return
    return(outliers)
}
