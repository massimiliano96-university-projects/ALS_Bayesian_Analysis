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

...
