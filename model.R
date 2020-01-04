longitudinal_data <- read.csv("alsfrs_new.csv", header = TRUE)
longitudinal_data <- longitudinal_data[,c(-1)]

longitudinal_data$Delta_FVC = longitudinal_data$Delta * longitudinal_data$FVC
longitudinal_data$Delta_URIC_ACID = longitudinal_data$Delta * longitudinal_data$URIC_ACID
longitudinal_data$Delta_CREATININE = longitudinal_data$Delta * longitudinal_data$CREATININE

id <- unique(longitudinal_data$subject_id)

n = 36830
m = 4142

longitudinal_data$counter = 0 ### number of obs for each ID
j=1
i = 1

while(i < dim(longitudinal_data)[1]){
  if(longitudinal_data$subject_id[i] == longitudinal_data$subject_id[j])
  {i = i+1}
  else
  {
    longitudinal_data$counter[j:(i-1)] = i-j
    j = i
    i = i+1
  }
}

# R & STAN are friends!
library(rstan)
library(coda)

# for plots
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(ggsci)


new_ID = longitudinal_data$ID
id <- longitudinal_data$subject_id[1]
cont <- 1
for(i in 1:n){
  if( longitudinal_data$subject_id[i]==id ) 
  {
    new_ID[i] = cont
  }
  else
  {
    id <- longitudinal_data$subject_id[i]
    cont <- cont + 1
    new_ID[i] = cont
  }
}

numerosity <- as.vector(table(longitudinal_data[,1])) ### number of obs for each patient


stan_data <- list(M = m, N = numerosity, 
                  Y = longitudinal_data$ALSFRS,
                  time = longitudinal_data$Delta,
                  AGE = longitudinal_data$Age,
                  Delta_FVC = longitudinal_data$Delta_FVC,
                  Delta_URIC_ACID = longitudinal_data$URIC_ACID,
                  Delta_CREATININE = longitudinal_data$Delta_CREATININE,
                  idx = as.numeric(new_ID))


fit1 <- stan(file = "stan_model.stan", 
             data = stan_data, 
             iter = 5000, thin = 10,
             chains = 1, warmup = 1000, 
             algorithm = 'NUTS', 
             diagnostic_file = 'diagnostic.txt', 
             verbose = TRUE,
             seed = 42)

save(fit1, file = 'stan_model.Rdata')

x11()
grid()

rstan::traceplot(fit1, pars = c('beta0', 'beta1', 'beta2', 'beta3', 'beta4'), inc_warmup = FALSE)

save(fit1, file = 'stan_model.Rdata')


# Diagnostic ------------------------------------------------------------

coda_chain <- As.mcmc.list(fit1, pars = c("beta0", "beta1", "beta2", "beta3", "beta4"))
summary(coda_chain)

# autocorrelation and plot
plot(fit1, plotfun = "stan_ac", pars = c( 'beta0', 'beta1', 'beta2', 'beta3', 'beta4'))

data.out=as.matrix(fit1) # trasform the mcmc.list into a matrix,
data.out=data.frame(data.out)     # or, better, into a dataframe (easiest to handle in R)
attach(data.out)
n.chain=dim(data.out)[1] 
n.chain # this is the final sample size (5000-1000)/5

beta.post <- data.out[,c('beta0','beta1','beta2', 'beta3', 'beta4')]
#posterior mean of the beta parameters
beta.bayes  <- apply(beta.post,2,"mean")
beta.bayes

chain <- beta.post[,c('beta2')]
#Divide the plot device in three sub-graph regions
#two square on the upper and a rectangle on the bottom
x11()
layout(matrix(c(1,2,3,3),2,2,byrow=T))
#trace-plot of the posterior chain
plot(chain,type="l",main="Trace plot of beta3")
# autocorrelation plot
acf(chain,lwd=3,col="red3",main="autocorrelation of FVC")
#Histogram
hist(chain,nclass="fd",freq=F,main="Posterior of FVC",col="gray") 
## Overlap the kernel-density 
lines(density(chain),col="blue",lwd=2)
## Posterior credible interval of beta3
quantile(chain,prob=c(0.025,0.5,0.975))


## Display the posterior credible interval on the graph
abline(v=quantile(chain,prob=c(0.025)),col="red",lty=2,lwd=2)
abline(v=quantile(chain,prob=c(0.5)),col="red",lty=1,lwd=2)
abline(v=quantile(chain,prob=c(0.975)),col="red",lty=2,lwd=2)
## Add the legend to the plot
legend("topright",legend=c("posterior median", "95% Credible bounds","kernel density smoother"),lwd=c(2,2,2), col=c("red","red","blue"),lty=c(1,2,1))


# PLOTTO UN PAZIENTE CASUALE
patient_1 <- longitudinal_data[1:12, c(6,4)]
x11()
plot(patient_1$Months, patient_1$ALSFRS_Total, type= "l", col= "red")
data_1 <-  longitudinal_data[1:12, c(6,7)]
y <- beta.post$beta0 + beta.post$beta1* data_1$Months + beta.post$beta2* data_1$Delta_FVC
data_pred <- cbind(y, data_1$Months)
data_pred <- as.data.frame(data_pred)
points(data_pred$V2 , data_pred$y, type = "l", col= "blue")
