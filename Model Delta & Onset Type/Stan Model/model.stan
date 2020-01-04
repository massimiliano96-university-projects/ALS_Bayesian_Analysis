data
{
	int<lower = 0> M;       // number of patients
	int<lower = 0> N[M];    // numerosities in each patient: array 
	vector<lower = 0, upper = 48>[sum(N)] Y;     // ALSFRS
	int<lower = 1, upper = M> idx[sum(N)];   // which patient the datum belongs to!
	vector<lower = 0, upper = 1>[sum(N)] dummy_bulbar; //dummy for bulbar or limb
	vector[sum(N)] days; //days
	vector[sum(N)] interaction; //interaction
}

parameters
{
	real theta1[M];		// array of the mean into each genus 
	real theta2[M];
	real<lower = 0> sigma2;
	real<lower = 0> sd_1;
	real<lower = 0> sd_2;
	vector[4] beta;
}

transformed parameters 
{
	vector[sum(N)] mu;

  	mu = beta[1] + beta[2]*dummy_bulbar + beta[3]*days + beta[4]*interaction;	
}

model
{
	// Prior:
	sigma2 ~ inv_gamma(2, 1.);
	sd_1 ~ inv_gamma(2, 1.);
	sd_2 ~ inv_gamma(2, 1.);
	beta[1] ~ normal(0., 5.); 
	beta[2] ~ normal(0., 1.);
	beta[3] ~ normal(0., 1.);
	beta[4] ~ normal(0., 1.);

	theta1 ~ normal(0, sd_1);  // iid sample
	theta2 ~ normal(0, sd_2);
	

	// Likelihood:
	for(i in 1:(sum(N)))
	{
		Y[i] ~ normal(mu[i] + theta1[idx[i]] + theta2[idx[i]]*days[i], 
							sigma2);
	}
}