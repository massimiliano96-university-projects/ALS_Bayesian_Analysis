# Data Preprocessing

After we get our Alsfrs-R in longitudinal way (see Make ALSFRS longitudinal folder), we need to build a dataset with other covariates,
fixed or longitudinal.

To do, we use left_join function with on the left ALSFRS-R, and we join by (ID) for fixed covariates, and by (ID, Delta) for longitudinal
covariates. 

In this way, if we got for example alsfrs with columns (ID, Delta, ALSFRS-R) and a longitudinal variable v: (ID, Delta, V)
If we call left_join(alsfrs,v) we get all the rows from alsfrs, and we match it by (ID, Delta) with v, and the result will be
a dataset with (ID, Delta, ALSFRS-R, V) and when a row from alsfrs doesn't match with v, we get a NA in column V.

## Interpolation for covariates

With this strategy, we can loose informations from the covariates v when the Delta of a patient does not match. For this reason 
we use a linear interpolation only for the longitudinal covariates (because of course the column of ALSFRS will be full) through two 
functions "join_alsfrs_interpolation" and "interpol".

## Threshold for NA

Even if we use an interpolation, some patients simply dont have that covariate and they have anyway a NA. We need to fill it with 
other techniques, but first of all we discard the columns with more than 35% of NA.

## Filling Method

For what concern fixed covariates, we use some tecniques to fill the NA, for example fill with the median.
For longitudinal covariates we use the package MICE (Multivariate Imputation by Chained Equations) with n=5 imputation.

The dataset resulting from this preprocessing is "alsfrs_riempito.csv"
