# Final model: Covariates and Cosine in Time


After the covariate selection, we choose to discard BMI_0, BMI_{diff}, Bilirubin, Glucose, Hematocrit, Hemoglobin, Urin PH, Riluzole and Sex.

![betas_model_final](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model5/images/betas_model_final.jpg)

For what concerne the coefficients estimate by the model (now with Normal priors as in the first models), we can see that
they are all statistically significative except from some lab data (like SGPT, SGOT and RBC) and Study Treatment, 
in accordance with what we found in the "Treatment Model".
While the coefficient of "Study Treatment" is perfectly centered in zero, the coefficient of Interaction between Study Treatment
and Delta contains zero in its 95 IC but it is close to the boundary lines; for this reason, we can left it in our model.

If we analyze the two errors that we have defined to evaluate these model, we can see that we can 

Errors:

| Error 1  |  |
| ------------- | ------------- |
| Training set  |  1.048792 |
| Train. with new obs |  3.582663 |
| Test set | 5.977326  |

| Error 2  |  |
| ------------- | ------------- |
| Training set  | 9.024221  |
| Train. with new obs | 20.45746   |
| Test set |  35.94814  |

If we focus on the errors for the second case, where we have the patients in the training set but we evaluate the score on some new temporal observation 
(that is what a doctor should do in a real case: monitoring a patient for a few months and try to predict how he can behave in the future), we have an error1 of 3.5, that is quite good
on a score that ranges from 0 to 48, and an error2 of 20, that means an uncertainty of plus or minus 8 points more or less (in the case where we consider a credible interval of 95%).

The error on the training set is very small (error 1 is 1.04) because we can predict well the behaviour of these patients using all the covariates in the fixed effect part but also
using the random effect part, that can be used to model the variability among the patients.

Finally, the error on some new patients, i.e. patients such that their ID is not present in the training dataset, is bigger than the previous two cases, and this is due to the fact
that now we can't sample directly from the random effect specifically of these patients, because we don't have it in our model, but we have to sample it from the posterior
distribution of the random effects that is basically a normal distribution with zero mean and a big variance (thanks to the eteronegeneity of our training patients). 

We will see some examples of prediction of train patient but on new temporal observation (that is, for the reason explained above, the most intresting prediction we can make) 
in the last section of this report.



