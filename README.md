# ALS Bayesian Analysis

(Authors: Blerta Begu, Daniele Ceccarelli and Massimiliano Riva)

The aim of this project is to study ALS (Amyotrophic Lateral Sclerosis) Progression with a Bayesian approach. We are using the PRO-ACT
database (https://nctu.partners.org/ProACT) 

## ALSFRS
We want to study the progression of the disease in terms of ALSFRS-R (ALS Functional Rating Score - Rev), that is a score used by doctors to evaluate ALS of a patient with 12 questions, rated from 0 to 4, for a maximum of 48 points.
The questions want to evaluate the functionality of different part of the body: more info at (https://www.outcomes-umassmed.org/ALS/sf12.aspx#scale).

## Data preparation

First, we need to work on alsfrs dataset in order to have it all in ALSFRS-Revised way. We end up with a dataset with columns:
- ID of patient
- Delta, time from time 0 (when patient enters in the dataset) when we evaluate ALSFRS-R
- ALSFRS-R: score.

Every patient has more than one temporal evaluation (and so more than one rows in our dataset), for this reason we can use a Longitudinal approach.
See [read me](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/make%20ALSFRS%20longitudinal/readme.md) in the folder "Make ALSFRS longitudinal" for details.

After that, we need to build a dataset with other covariates, both fixed and longitudinal. We include some data of patients like:
- Personal data: Age, Sex, BMI
- Disease data: Onset syte of ALS and Onset Delta (time from Time 0 when the patient start to see the symptoms)
- Lab data: some lab test like FVC, SGPT, SGOT etc.
- Medication data: Riluzole and Treatment.

See [read me](https://github.com/massimiliano96/ALS_Bayesian_Analysis/tree/master/Data%20Preprocessing) in the folder "Data Preprocessing" for details.

# Bayesian Approach: Longitudinal Data and First Mixed Effect Model
We use a Mixed Effect model with i = 1,...,m (num of patients), and j = 1,...,n_i (with n_i number of obs for patient i).
![images approach](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model1/images/bayesian%20approach.jpg)

We can see the fixed effect part (X_ij * beta) and the random effect part (Z_ij * theta_i), where theta_i is the vector of random
coeff for the patient i.
We assume as priors a Normal prior with fixed variance for the betas, while for theta_i a normal with a normal-invgamma model.

We use as covariate matrix X_ij both the time (for example Delta and later also Cosine_Delta) and the other covariates that we add 
in Data Preprocessing (both fixed and longitudinal); while for the random effect part we are going to use only time variables.

# Models
## How can we evaluate our models?

We have split our dataset in 3 parts: a training set and two different test sets.
Firstly, we select a group of patient ID to use them in our model ("alsfrs_train"). Another group of ID is used as test set ("alsfrs_test).
Then, we take our training set and we split it in two part:
For some patient (50%) in training set, we take the second half temporal evaluation and we put it in a new test set, 
and we call it : "alsfrs_train_patient_test_observation".

We have done this second split for two reasons:
1) we can test our model on two different test set, with different characteristics;
2) mixed effect model predict in a different way if the patient ID is present or not in the training set.

To evaluate our model, in addiction to standard Bayesian methods like WAIC and LOO, we compute two different prediction error,
and we evaluate on the three different dataset that we have (alsfrs_train, alsfrs_train_patient_test_observation, alsfrs_test).

The first error is just a mean absolute error: 
        mean( | mean(Y_predict) - Y_true | )

The second error takes into account also the variability of prediction: 
        mean( | mean(Y_predict) - Y_true | + |lenght(IC_95)| )

## Model 1: Delta

In our first model, we use as covariate only Delta. 
In particular: 
mean(ALSFRS_ij) = beta_0 + beta_1 * Delta_ij + theta_0i + theta_1i* Delta_ij

See [read me](https://github.com/massimiliano96/ALS_Bayesian_Analysis/tree/master/model1) for details and results.

## Model 2: Delta + Onset_site

In second model, we add our first covariate, in this case a factor varible: Onset_Site, that can be equal to Bulbar or Limb.
In the medicine literature, they state that the Bulbar type of ALS is more aggressive than the Limb one.

![images bulbar vs limb](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model2/images/bulbar_vs_limb.png)

mean(ALSFRS_ij) = beta_0 + beta_1 * Delta_ij + beta_2 * Onset_Site_LIMB + beta_3 * Delta_ij * Onset_Site_LIMB +
                  theta_0i + theta_1i * Delta_ij

See [read me](https://github.com/massimiliano96/ALS_Bayesian_Analysis/edit/master/model2/readme.md) for results.

| Error 1  |  |
| ------------- | ------------- |
| Training set  | 1.353407  |
| Train. with new obs | 3.808556  |
| Test set | 6.555238  |

| Error 2  |  |
| ------------- | ------------- |
| Training set  | 10.62259 |
| Train. with new obs | 20.34279  |
| Test set | 40.21576  |

## Model 3: Delta + Onset_site + Medication (Riluzole and Treatment)
At this step we wanted to investigate if some clinical behaviour have a statistical evidence effect on the progression of the disease. With this aim, we added the following covariates:

### Treatment:
The exact medications used in the trials are not specified, as part of our effort to avoid identification of the patients involved. Information is available as to whether any individual patient received medication or placebo.

### Riluzole :
Riluzole is the only drug approved for treating Amyotrophic Lateral Sclerosis (ALS). Information is included about its use.

See [read me](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model3/readme.md) for results.

| Error 1  |  |
| ------------- | ------------- |
| Training set  | 1.353605  |
| Train. with new obs | 3.7856  |
| Test set | 20.29393  |

| Error 2  |  |
| ------------- | ------------- |
| Training set  |  10.61728 |
| Train. with new obs | 20.29393   |
| Test set | 40.2572   |

## Model 4: Covariates selection

And finally we add all the covariates that we select in "Data Preprocessing". With more than 20 betas, we have to face the problem
of overfitting due to the number of degree of freedom of our model.
For this reason we use a regularization in the prior for the beta, 
in particular we use the Horseshoe Prior (Carvalho, Polson, and Scott 2009).

If we want to understand how this regularization works, 
we can see the comparison of the unit ball from classical regul (Ridge, Lassu, Cauchy) and the Horseshoe.
![images regularization](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model4/images/horseshoe.jpeg)

The difference here is that in the Horseshoe reg. the unit ball contain all the axes; for this reason, we can delete the covariates
that will have alll the density around zero.

### Covariates
We can see the result for the betas
![betas model4](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model4/images/betas_model4.jpg)


There are some variables that have a density all in zero, like BMI_0, BMI_diff, Sex, Glucose etc.

In the next model (model 5) we don't use these variables, in order to have a simpler model and to save computational time.

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

## Model 5: Covariates + Cosine

After the selection made in model 4, we add some cosine splines in time, both in fixed and random part, 
and we discard BMI_0, BMI_{diff}, Bilirubin, Glucose, Hematocrit, Hemoglobin, Urin PH, Riluzole and Sex.


For a patient i, at time j :
- cos_delta_1_ij = cos(pi * (t_ij - t_min,i) / (t_max,i - t_min,i))
- cos_delta_2_ij = cos(2 * pi * (t_ij - t_min,i) / (t_max,i - t_min,i))
- cos_delta_3_ij = cos(3 * pi * (t_ij - t_min,i) / (t_max,i - t_min,i))

These new terms can capture the sinusoidal behaviour of some patients, that in particular have a first period of slow decrease,
then a fast decrease, and finally a "stabilization".

![plot patient cos](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model5/images/plot_pazienti_cos.jpeg)

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



# Prediction on patient from training set on NEW observation

Probably the most useful prediction from this method could be a prediction for a patient in the training set but on 
new observation (later in time).
With this prediction, we can sample from the posterior of thetas of that specific patient (in addiction to the fixed effect betas)
and we can take into account the variability among the different patients.

Let's see two example:
![paziente1](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model5/images/prediction_2_patients.jpeg)

We can see in blue the evaluation in the training set, while in green the new observations. The prediction are quite good, the estimate
are just a bit bigger than the real value.

# Patient far from "the mean"

We can use this type of model (Mixed Effect Model) to found patients that behave very differently from the rest.
If for example we analyze their random effect part (in mean, since we have a posterior sample), we can see these type of plot:

![images_pair](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model5/images/random_effect_pairs.jpeg)

Here we can detect, for example patients with anomalous random intercept:
![random_intercept_outlier](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model5/images/outliers_intercept.jpeg)

for others, see:
![images](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model5/images)


