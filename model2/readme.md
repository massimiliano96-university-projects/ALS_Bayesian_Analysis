# Model Delta & Onset Type

With this model we want to investigate the difference between the two different types of ALS onset location: Bulbar or Limb.
In the medicine literature, we can see that the Bulbar type is more aggressive than the Limb one; we want to prove this.

Here is a plot of 100 patient for both groups:
![images bulbar vs limb](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model2/images/bulbar_vs_limb.png)

Bulbar group seems to have a higher intercept than Limb, but also a faster decrease: we want to investigate if being in Bulbar
group means to have a fast progression of the disease.

# Bayesian Analysis

We use a Bayesian Mixed Effect Model, where the groups are made by each patient

We have introduced a dummy variable for Bulbar type: I_bulbar(patient_i) = 1 if patient_i has a Bulbar ALS, 0 otherwise

## formula

ALSFRS_(patient_i, time_j) = \beta_0 + \beta_1*Delta(j) + \beta_2* I_bulbar(i) + \beta_3*Delta(j)*I_bulbar(i) + \theta_1i + \theta_2i*Delta(j)

For the betas, we assume a Normal prior with fix variances, instead for theta we assume a normal prior of N(0,tau^2), where 
tau comes from a inv-gamma (or a cauchy) prior. We've tried also to model a correlation between the two theta parameters (while the betas are assumed to be independent). A LKJ correlation distribution was used in one of the three different models, using Cholesky
factorization.

## Prior choice 

Model with:
Normal prior for the betas
inv-gamma(2,1) prior for thetas and no correlation among them

# Results 


![images traceplot](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model2/images/traceplot_betas.png)

These are the traceplots of our betas, we can see that beta2 is bigger than zero and beta3 less than zero, as we can expect.

More in detail, if we analyze the posterior density of beta3 we can see that its 95% Credible Interval Bounds are both less than zero,
and this can confirm our hypotesis that Bulbar patient has in general a faster decrease of the disease.
![images traceplot](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model2/images/beta3.png)
The difference is more or less -0.01 * day in mean, that is -0.3 ALSFRS-R point/month, quite significant.

# Error
