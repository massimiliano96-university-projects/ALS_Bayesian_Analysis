data
{
	int<lower = 0> M;       // number of patients
	int<lower = 0> N[M];    // numerosities in each patient: array 
	vector<lower = 0>[sum(N)] Y;     // ALSFRS
	int<lower = 1, upper = M> idx[sum(N)];   // which patient the datum belongs to!
	vector[sum(N)] time; 
	vector[sum(N)] Delta_FVC; 
	vector[sum(N)] Delta_URIC_ACID; 
	vector[sum(N)] Delta_CREATININE;
	vector[sum(N)] AGE;
	//vector[sum(N)] SEX;
}

parameters
{
	real theta1[M];		// array of the mean into each genus 
	real theta2[M];
	real<lower = 0> sigma2;
	//real<lower = 0> tau2;
	real beta0;
	real beta1;
	real beta2;
	real beta3;
	real beta4;
	real beta5;
	//real beta6;
}

transformed parameters 
{
	vector[sum(N)] mu;
  	mu = beta0 + beta1*AGE + beta2*time + beta3*Delta_FVC + beta4*Delta_URIC_ACID + beta5*Delta_CREATININE ; //+ beta2*SEX
}

model
{
	// Prior:
	sigma2 ~ inv_gamma(2., 1.);
	beta0 ~ normal(0., 1.); 
	beta1 ~ normal(0., 1.);
	beta2 ~ normal(0., 1.);
	beta3 ~ normal(0., 1.);
	beta4 ~ normal(0., 1.);
	beta5 ~ normal(0., 1.);
	//beta6 ~ normal(0., 1.);


	//tau2 ~ inv_gamma(2., 1.);
	theta1 ~ normal(0, 1);  // iid sample
	theta2 ~ normal(0, 0.1);

	// Likelihood:
	for(i in 1:(sum(N)))
	{
		Y[i] ~ normal(mu[i] + theta1[idx[i]] + theta2[idx[i]]*time[i], 
							pow(sigma2, 0.5));
	}	

}
