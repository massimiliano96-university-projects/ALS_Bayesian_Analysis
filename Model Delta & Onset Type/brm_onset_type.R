library(brms)

library(rstan)
library(coda)

library(dplyr)



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



######### plot of the two different groups: are they different? ###################

set.seed(1)
n = 100

id_0 = unique(alsfrs[which(alsfrs$dummy_bulbar == 0),]$ID) ###ID which Als type = LIMB
id_1 = unique(alsfrs[which(alsfrs$dummy_bulbar == 1),]$ID) ###ID which Als type = BULBAR

sample_0 = sample(id_0, n) ### LIMB sample of n patients
sample_1 = sample(id_1, n) ### BULBAR sample of n patients


#plot original data
x11()
par(mfrow = c(1,2))

#### DUMMY = 0
plot(alsfrs[which(alsfrs$ID == sample_0[1]),]$Delta, 
     alsfrs[which(alsfrs$ID == sample_0[1]),]$ALSFRS_R_Total, type = "b", col = "red", ylim = c(0,48),xlim = c(0,700), main = "Limb, dummy=0", xlab = "days", ylab = "ALSFRS-R")

for(i in 2:n){
  points(alsfrs[which(alsfrs$ID == sample_0[i]),]$Delta, 
         alsfrs[which(alsfrs$ID == sample_0[i]),]$ALSFRS_R_Total, type = "b", col = "red")
}

#### DUMMY = 1
plot(alsfrs[which(alsfrs$ID == sample_1[1]),]$Delta, 
     alsfrs[which(alsfrs$ID == sample_1[1]),]$ALSFRS_R_Total, type = "b", col = "forestgreen", ylim = c(0,48), xlim = c(0,700), main = "Bulbar, dummy = 1", xlab = "days", ylab = "ALSFRS-R")

for(i in 2:n){
  points(alsfrs[which(alsfrs$ID == sample_1[i]),]$Delta, 
         alsfrs[which(alsfrs$ID == sample_1[i]),]$ALSFRS_R_Total, type = "b", col = "forestgreen")
}
########## Comments on the plot: ##############
### we can see a small difference in the distribution of 
### the two groups, Bulbar disease
### seems to be more aggressive than Limb 
##################################################



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

### new variable with ID a number from 1 to m_test
new_ID_test = rep(0,m_test)
i=1
j=1
while(i<=m_test){
  new_ID_test[j:(j+numerosity[i]-1)] = i
  j = j+numerosity[i]
  i = i+1
}


############### we can use this package -> brms
########## it use STAN and it's useful to model also the correlation prior!
fit1 <- brm(data = alsfrs_train,
              family = gaussian,
              ALSFRS_R_Total ~ 1 + Delta + dummy_bulbar + 
                Delta*dummy_bulbar + (1 + Delta | ID),
              
              prior = c(prior(normal(0, 5), class = Intercept),
                        prior(normal(0, 1), class = b),
                        prior(cauchy(0, 1), class = sd),                ### cauchy prior
                        prior(cauchy(0, 1), class = sigma),             ###
                        prior(lkj_corr_cholesky(1.5), class = cor)),    ### prior to the correlation
                                                                        ### (lkj)
              iter = 20000, warmup = 1000, chains = 1, cores = 4,
              control = list(adapt_delta = .975, max_treedepth = 20),
              seed = 190831)

save(fit1, file = "brm_model_onset_type.Rdata")

#fit2 <- brm(data = alsfrs_train,
#              family = gaussian,
#              ALSFRS_R_Total ~ 1 + Delta + dummy_bulbar + 
#                Delta*dummy_bulbar + (1 + Delta | ID),
#              
#              prior = c(prior(normal(0, 5), class = Intercept),
#                        prior(normal(0, 1), class = b),
#                        prior(inv_gamma(2, 1), class = sd),             ### inv_gamma prior
#                        prior(inv_gamma(2, 1), class = sigma),          ###
#                        prior(lkj_corr_cholesky(1.5), class = cor)),    ### prior to the correlation
#                                                                        ### (lkj)
#              iter = 20000, warmup = 1000, chains = 1, cores = 4,
#              control = list(adapt_delta = .975, max_treedepth = 20),
#              seed = 190831)
#
#save(fit2, file = "brm_model_onset_type_inv_gamma.Rdata")


#fit3 <- brm(data = alsfrs_train,
#              family = gaussian,
#              ALSFRS_R_Total ~ 1 + Delta + dummy_bulbar + 
#                Delta*dummy_bulbar + (1 + Delta || ID),          ##### NO CORRELATION ####
#              
#              prior = c(prior(normal(0, 5), class = Intercept),
#                        prior(normal(0, 1), class = b),
#                        prior(inv_gamma(2, 1), class = sd),
#                        prior(inv_gamma(2, 1), class = sigma),
#                       
#              iter = 20000, warmup = 1000, chains = 1, cores = 4,
#              control = list(adapt_delta = .975, max_treedepth = 20),
#              seed = 190831)
#
#save(fit3, file = "brm_model_onset_type_inv_gamma_no_corr.Rdata")


load("brm_model_onset_type.Rdata")


summary(fit1)

x11()
plot(fit1)

marginal_effects(fit1)

### to recover the stancode -> for complex model
code = stancode(fit1)
code
standata = standata(fit1)

loo1 = loo(fit1) ### leave-one-out cross-validation 
