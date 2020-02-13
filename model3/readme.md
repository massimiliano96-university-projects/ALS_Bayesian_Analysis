# Model with Onset Site + Medication (Riluzole and New Treatment)

At this step we wanted to investigate if some clinical behaviour have a statistical evidence effect on the progression of the disease. With this aim, we added the following covariates:

## Treatment:
The exact medications used in the trials are not specified, as part of our effort to avoid identification of the patients involved. Information is available as to whether any individual patient received medication or placebo.

## Riluzole :
Riluzole is the only drug approved for treating Amyotrophic Lateral Sclerosis (ALS). Information is included about its use.

## Model:

$\beta$

ALSFRS_(patient_i, time_j) = \beta_0 + \beta_1Delta(j) + \beta_2 I_bulbar(i) + \beta_3*Delta(j)I_bulbar(i) + \theta_1i + \theta_2iDelta(j)

![alt text](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model3/images/plot_b_Treatment.%20jpeg)

![alt text](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model3/images/plot_b_Treatment_interaction.%20jpeg)

![alt text](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model3/images/plot_b_Riluzole.%20jpeg)

![alt text](https://github.com/massimiliano96/ALS_Bayesian_Analysis/blob/master/model3/images/plot_b_Riluzole_interaction.%20jpeg)
