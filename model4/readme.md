# Covariates selection: HORSESHOE PRIOR

And finally we add all the covariates that we select in "Data Preprocessing". With more than 20 betas, we have to face the problem
of overfitting due to the number of degree of freedom of our model.
For this reason we use a regularization in the prior for the beta, 
in particular we use the Horseshoe Prior (Carvalho, Polson, and Scott 2009).

