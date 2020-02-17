library(brms)
library(rstan)
library(coda)
library(dplyr)



alsfrs <- read.csv("alsfrs_riempito.csv", header = TRUE) ###we're gonna use ALSFRS-R

numerosity <- as.vector(table(alsfrs[,1])) ### number of obs for each patient

m = length(unique(alsfrs$ID)) 

alsfrs$counter = 0
i=1
j=1
while(i<=m){
  alsfrs$counter[j:(j+numerosity[i]-1)] = numerosity[i]
  j = j+numerosity[i]
  i = i+1
}

alsfrs = alsfrs[which(alsfrs$counter >= 5),]
#alsfrs = alsfrs[which(alsfrs$counter > 1),]

alsfrs = alsfrs[,c(-22)] ### delete counter, change number

numerosity <- as.vector(table(alsfrs[,1])) ### number of obs for each patient
n1 = sum(numerosity[1:500]) ###
########## subset -> training set
alsfrs_train = alsfrs[1:n1,] ### 500 patients

ID_train = unique(alsfrs_train$ID)
n_train = dim(alsfrs_train)[1]
m_train = length(ID_train)
numerosity_train <- as.vector(table(alsfrs_train[,1])) ### number of obs for each patient

########## subset -> test set
n2 = sum(numerosity[501:800])
alsfrs_test = alsfrs[(n1+1):(n1+n2),] ### 300 patients

ID_test = unique(alsfrs_test$ID)
n_test = dim(alsfrs_test)[1]
m_test = length(ID_test)
numerosity_test <- as.vector(table(alsfrs_test[,1])) ### number of obs for each patient



##################### split some of the train patients in train and set observations
alsfrs_train$train_obs = 1
num_train <- as.vector(table(alsfrs_train[,1])) 
m_train = length(unique(alsfrs_train$ID)) 

i=1
j=1
while(i <= floor(m_train/2)){
  alsfrs_train$train_obs[(j + floor(num_train[i]/2)):(j+num_train[i]-1)] = 0
  j = j+num_train[i]
  i = i+1
}

alsfrs_train_patient_test_observation = alsfrs_train[which(alsfrs_train$train_obs == 0), c(-22)]
alsfrs_train = alsfrs_train[which(alsfrs_train$train_obs == 1), c(-22)]



########## fit model with brm #####################
fit1 <- brm(data = alsfrs_train,
            family = gaussian,
            ALSFRS_R_Total ~ 1 + Delta + Onset_Site + Onset_Delta + 
              Age + FVC + BMI_0 + BMI_diff + SGPT + SGOT + BILIRUBIN + 
              GLUCOSE + HEMATOCRIT + HEMOGLOBIN + RBC + WBC + URINE_PH +
              RILUZOLE + Study_Treatment + Sex + 
              Delta*Onset_Site + Delta*Study_Treatment + Delta*RILUZOLE +
              (1 + Delta || ID),
            
            prior = c(prior(normal(0, 10), class = Intercept),
                      prior(horseshoe(1), class = b),
                      prior(inv_gamma(2, 1), class = sd),
                      prior(inv_gamma(2, 1), class = sigma)),    
            iter = 20000, warmup = 5000, chains = 1, cores = 4,
            control = list(adapt_delta = .9, max_treedepth = 15),
            seed = 190831)
save(fit1, file = "brm_model_covariates_500_horseshoe3.Rdata")


load("brm_model_covariates_500_horseshoe3.Rdata")


summary(fit1)

x11()
plot(fit1)

marginal_effects(fit1)

### to recover the stancode -> for complex model
code = stancode(fit1)
code
standata = standata(fit1)

loo1 = loo(fit1) ### leave-one-out cross-validation 



################## model checking on training dataset ###################
x11()
par(mfrow = c(2,3))
for(n in 11*(1:6)){
  new_train = alsfrs_train[which(alsfrs_train$ID == ID_train[n]), ]
  data = predict(fit1, newdata = new_train, allow_new_levels = TRUE)
  
  plot(new_train$Delta, data[,1], ylim = c(0,48), type = "b", col = "red", pch = 20, xlab = "days", ylab = "ALSFRS")
  points(new_train$Delta, alsfrs_train[which(alsfrs_train$ID == ID_train[n]), c(3)], type = "b", col = "blue", pch = 17)
  points(new_train$Delta, data[,3], ylim = c(0,48), type = "b", col = "orange", pch = 20)
  points(new_train$Delta, data[,4], ylim = c(0,48), type = "b", col = "orange", pch = 20)
}


################## ###################
x11()
par(mfrow = c(2,3))
for(n in 11*(1:6)){
  new_train = alsfrs_train_patient_test_observation[which(alsfrs_train_patient_test_observation$ID == ID_train[n]), ]
  data = predict(fit1, newdata = new_train, allow_new_levels = TRUE)
  
  plot(new_train$Delta, data[,1], ylim = c(0,48), type = "b", col = "red", pch = 20, xlab = "days", ylab = "ALSFRS")
  points(new_train$Delta, alsfrs_train_patient_test_observation[which(alsfrs_train_patient_test_observation$ID == ID_train[n]), c(3)], type = "b", col = "blue", pch = 17)
  points(new_train$Delta, data[,3], ylim = c(0,48), type = "b", col = "orange", pch = 20)
  points(new_train$Delta, data[,4], ylim = c(0,48), type = "b", col = "orange", pch = 20)
}



####################### prediction on new patients from test set ########################
x11()
par(mfrow = c(2,3))
for(n in 12*(1:6)){
  new_test = alsfrs_test[which(alsfrs_test$ID == ID_test[n]), ]
  data = predict(fit1, newdata = new_test, allow_new_levels = TRUE)
  
  plot(new_test$Delta, data[,1], ylim = c(0,48), type = "b", col = "red", pch = 20, xlab = "days", ylab = "ALSFRS")
  points(new_test$Delta, alsfrs_test[which(alsfrs_test$ID == ID_test[n]), c(3)], type = "b", col = "blue", pch = 17)
  points(new_test$Delta, data[,3], ylim = c(0,48), type = "b", col = "orange", pch = 20)
  points(new_test$Delta, data[,4], ylim = c(0,48), type = "b", col = "orange", pch = 20)
}
###################################################################################



######## errors ##########
data = predict(fit1, newdata = alsfrs_train)
mean(abs(data[,1] - alsfrs_train$ALSFRS_R_Total))

data = predict(fit1, newdata = alsfrs_train_patient_test_observation)
mean(abs(data[,1] - alsfrs_train_patient_test_observation$ALSFRS_R_Total))

data = predict(fit1, newdata = alsfrs_test, allow_new_levels = TRUE)
mean(abs(data[,1] - alsfrs_test$ALSFRS_R_Total))

# 1.358297
# 3.639182
# 6.036255

hist(abs(data[,1] - alsfrs_test$ALSFRS_R_Total))
