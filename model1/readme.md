# Model Delta
This is the first model of this project, the aim is investigating the trend of ALSFRS with respect ot the time (Delta). Of course we expect the ALSFRS as a decreasing function of time.

## Model
We adopted a Mixed effect framework for ore models:

ALSFRS_(patient_i, time_j) = \beta_0 + \beta_1Delta(j) + \theta_1i + \theta_2i*Delta(j)

here the betas are the fixed coefficients and the thetas the random ones.

## Prior

For the betas, we assume a Normal prior with fix variances, instead for theta we assume a normal prior of N(0,tau^2), where tau comes from a inv-gamma prior. 

## Results

![alt text](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model1/images/plot_model1_beta2.%20jpeg)

the posterior credible interval is totally negative, hence we deduce that our covariate is significative, and ALSFRS decrease with time.

## Random Effect Outliers

We discovered the existence of some patients for which the thetas are very large:

![alt text](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model1/images/random_effect_outliers.%20jpeg)

we plotted these patients:

![alt text](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model1/images/outliers_score.%20jpeg)

As you can see, there is a very big variability in the intercept and the slope, this fact justify the addition of new covariates, let's proceed with Model2.
