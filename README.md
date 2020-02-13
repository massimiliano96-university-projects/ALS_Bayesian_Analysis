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
![images approach](
https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/images/bayesian%20approach.jpg)

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

![images bulbar vs limb](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/Model%20Delta%20%26%20Onset%20Type/images/bulbar_vs_limb.png)

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

## Model 4: Covariates selection

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
