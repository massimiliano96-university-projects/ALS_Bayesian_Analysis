#Model_Delta_Onset_Type

With this model we want to investigate the difference between the two different types of ALS onset location: Bulbar or Limb.
In the medicine literature, we can see that the Bulbar type is more aggressive than the Limb one; we want to prove this.

We use a Mixed Effect Model, where the groups are made by each patient:

We have introduced a dummy variable for Bulbar type: I_bulbar(patient_i) = 1 if patient_i has a Bulbar ALS, 0 otherwise

ALSFRS_(patient_i, time_j) = beta0 + beta1*I_bulbar(i) + beta2*Delta(j) + beta3*Delta(j)*I_bulbar(i) + theta1_i + theta2_i*Delta(j)
