# Model Delta
This is the first model of this project, the aim is investigating the trend of ALSFRS with respect ot the time (Delta). Of course we expect the ALSFRS as a decreasing function of time.

## Model
We adopted a Mixed effect fraimwork for ore models:

ALSFRS(t) = beta_0 + beta_1 * t + theta_0 + theta_1* t + epsilon

here the betas are the fixed coefficients and the thetas the random ones.

## Prior

we chose a normal - inv-gamma models

## Results

![alt text](https://github.com/massimiliano96/ALS_Bayesian_Analysis/tree/master/model1/images/plot_model1_beta2.jpeg)
