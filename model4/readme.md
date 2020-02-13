# Covariates selection: HORSESHOE PRIOR

And finally we add all the covariates that we select in "Data Preprocessing". With more than 20 betas, we have to face the problem
of overfitting due to the number of degree of freedom of our model.
For this reason we use a regularization in the prior for the beta, 
in particular we use the Horseshoe Prior (Carvalho, Polson, and Scott 2009).

If we want to understand how this regularization works, 
we can see the comparison of the unit ball from classical regul (Ridge, Lassu, Cauchy) and the Horseshoe.
![images regularization](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/images/horseshoe.jpeg)

The difference here is that in the Horseshoe reg. the unit ball contain all the axes; for this reason, we can delete the covariates
that will have alll the density around zero.

# Covariates
We can see the result for the first betas
![betas model4](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/images/betas_model4.jpg)

and for the betas of interaction between Delta and Onset_site, Medication and Riluzole
![beta2_model4](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/images/betas2_model4.jpg)

There are some variables that have a density all in zero, like BMI_0, BMI_diff, Sex, Glucose etc.

In the next model (model 5) we don't use these variables, in order to have a simpler model and to save computational time.

# Errors
| Error 1  |  |
| ------------- | ------------- |
| Training set  | 1.358297  |
| Train. with new obs | 3.639182  |
| Test set | 6.036255  |

| Error 2  |  |
| ------------- | ------------- |
| Training set  | 10.56538 |
| Train. with new obs | 19.77135 |
| Test set | 37.31195  |
