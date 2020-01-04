library(rstan)
library(coda)
library(dplyr)
library(loo)



alsfrs <- read.csv("alsfrs_R_complete_longitudinal.csv", header = TRUE) ###we're gonna use ALSFRS-R
alsfrs <- alsfrs[,c(-4)] ### delete the previous counter
onset_type <- read.csv("alstype_final.csv") ### We need dummy_bulbar 

alsfrs <- inner_join(alsfrs,onset_type) #add dummy_bulbar

names(alsfrs)[1] <- "ID"
names(alsfrs)[2] <- "Delta"

alsfrs <- na.omit(alsfrs)

###alsfrs = alsfrs[which(alsfrs$Delta>=0),] ###not necessary

ID = unique(alsfrs$ID)
n = dim(alsfrs)[1]
m = length(ID)
numerosity <- as.vector(table(alsfrs[,1])) ### number of obs for each patient

alsfrs$counter = 0
i=1
j=1
while(i<=m){
  alsfrs$counter[j:(j+numerosity[i]-1)] = numerosity[i]
  j = j+numerosity[i]
  i = i+1
}
alsfrs = alsfrs[which(alsfrs$counter > 5),]



#####################################################################
######### up to now, we use a small subset of patients -> computational time problem
######### we will use the complete dataset soon



numerosity <- as.vector(table(alsfrs[,1])) ### number of obs for each patient
n1 = sum(numerosity[1:500])
########## subset -> training set
alsfrs_train = alsfrs[1:n1,] ### 500 patients

ID_train = unique(alsfrs_train$ID)
n_train = dim(alsfrs_train)[1]
m_train = length(ID_train)
numerosity_train <- as.vector(table(alsfrs_train[,1])) ### number of obs for each patient

########## subset -> test set
n2 = sum(numerosity[501:600])
alsfrs_test = alsfrs[(n1+1):(n1+n2),] ### 100 patients

ID_test = unique(alsfrs_test$ID)
n_test = dim(alsfrs_test)[1]
m_test = length(ID_test)
numerosity_test <- as.vector(table(alsfrs_test[,1])) ### number of obs for each patient



### new variable with ID a number from 1 to m_train
new_ID_train = rep(0,m_train)
i=1
j=1
while(i<=m_train){
  new_ID_train[j:(j+numerosity[i]-1)] = i
  j = j+numerosity[i]
  i = i+1
}


stan_data <- list(M = m_train, N = numerosity_train, 
                  Y = alsfrs_train$ALSFRS_R_Total,
                  dummy_bulbar = alsfrs_train$dummy_bulbar,
                  days = alsfrs_train$Delta,
                  interaction = alsfrs_train$dummy_bulbar*alsfrs_train$Delta,
                  idx = new_ID_train
)


fit1 <- stan(file = "model.stan", 
             data = stan_data, 
             iter = 20000, thin = 5,
             chains = 1, warmup = 1000,
             cores = 4,
             algorithm = 'NUTS',
             control = list(adapt_delta = .975, max_treedepth = 20),
             verbose = TRUE,
             seed = 42)


save(fit1, file = 'model_stan.Rdata')

load("model_stan.Rdata")
summary(fit1)
plot(fit1, pars = c("beta[4]"))

data.out=as.matrix(fit1) # trasform the mcmc.list into a matrix,
data.out=data.frame(data.out)     # or, better, into a dataframe (easiest to handle in R)

chain = data.out[, c("sigma2", "sd_1", "sd_2", "beta.1.", "beta.2.", "beta.3.", "beta.4.")]
#plot
acf(chain[,7],lwd=3,col="red3",main="autocorrelation of beta4")


#### posterior prediction for test dataset!
beta_post = chain[,c(4:7)]
sigma_post = chain[,c(1)]
sd_post = chain[,c(2:3)]

test = alsfrs_test
n = dim(test)[1]


# Function for simulating y based on new x
gen_y_test <- function(test,n,iteration) {
  
  lin_comb = matrix(rep(0,n*iteration), iteration, n)   ### y for all new patients
  
  m = length(unique(test[,1])) ### number of new patient
  
  numerosity = as.vector(table(test[,1]))
  n_prec = 1
  
  for(i in 1:m){   ### for every new patient separately
    n_i = numerosity[i]
    lin_comb_i = matrix(rep(0,n_i*iteration), iteration, n_i)   ### y of new patient
    x = test[n_prec:(n_prec+n_i-1),] ### covariates of new patient
    
    for(k in 1:iteration){
      sigma = sample(sigma_post, size = n_i)
      sd1 = sample(sd_post[,1], size = n_i)
      sd2 = sample(sd_post[,2], size = n_i)
      theta1 = rep(0,n_i)
      theta2 = rep(0,n_i)
      
      for(j in 1:n_i){
        theta1[j] = rnorm(1,0,sd1[j])
        theta2[j] = rnorm(1,0,sd2[j])
      }
      
      lin_comb_i[k,1:n_i] <- sample(beta_post[,1], size = n_i) + ###beta1
      sample(beta_post[,2], size = n_i)*x[,4] +   ###beta2*Dummy_bulb
      sample(beta_post[,3], size = n_i)*x[,2] +   ###beta3*Days
      sample(beta_post[,4], size = n_i)*(x[,2]*x[,4]) +   ###beta4*Dummy_bulb*Days
      theta1 + theta2*x[,2] ###theta1 + theta2*Days
      
      error = rep(0, n_i) ###epsilon
      for(j in 1:n_i){
        error[j] = rnorm(1, 0, sigma[j])
      }
      
      lin_comb_i[k,1:n_i] = lin_comb_i[k,1:n_i] + error
    }
    lin_comb[1:iteration,n_prec:(n_prec+n_i-1)] = lin_comb_i[1:iteration, 1:n_i ] ###put the obs
                                                          ### of a single patient
                                                          ### in the final container
    n_prec = n_prec + n_i
  }  
  return(lin_comb)
}

# Run the function on x_test
set.seed(56)
y_pred_r <- gen_quant_r(test,n,1000)
mean_y = colMeans(y_pred_r)
true_y = test$ALSFRS_R_Total
quantile_1_y = apply(y_pred_r,2,quantile,c(.25))   ### lower bound quantile 
quantile_2_y = apply(y_pred_r,2,quantile,c(.75))  ### upper bound quantile

error = abs(mean_y - true_y)
hist(error)


###### plot result for a test patient
nn = 45 #number of test patient, between 1 and 100!

n1 = ifelse(nn == 1, 1, sum(numerosity[1:(nn-1)])+1)
n2 = sum(numerosity[1:nn])

plot(test$Delta[n1:n2], true_y[n1:n2], pch = 20, col = "blue", ylim = c(0,48))
points(test$Delta[n1:n2], mean_y[n1:n2], type = "b", pch = 20, col = "red", ylim = c(0,48))
points(test$Delta[n1:n2], quantile_1_y[n1:n2], type = "b", pch = 20, col = "yellow")
points(test$Delta[n1:n2], quantile_2_y[n1:n2], type = "b", pch = 20, col = "yellow")
