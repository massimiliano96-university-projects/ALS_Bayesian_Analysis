# Model Delta & Onset Type

With this model we want to investigate the difference between the two different types of ALS onset location: Bulbar or Limb.
In the medicine literature, we can see that the Bulbar type is more aggressive than the Limb one; we want to prove this.

We use a Mixed Effect Model, where the groups are made by each patient:

We have introduced a dummy variable for Bulbar type: I_bulbar(patient_i) = 1 if patient_i has a Bulbar ALS, 0 otherwise

ALSFRS_(patient_i, time_j) = \beta_0 + \beta_1 * I_bulbar(i) + \beta_2*Delta(j) + \beta_3*Delta(j)*I_bulbar(i) + \theta_1i + \theta_2i*Delta(j)
