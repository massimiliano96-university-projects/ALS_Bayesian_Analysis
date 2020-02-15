library(brms)
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

alsfrs = alsfrs[,c(1,2,3,4,5,6,7,10,11,16,17,20)]

numerosity <- as.vector(table(alsfrs[,1])) ### number of obs for each patient

m = length(unique(alsfrs$ID)) 

alsfrs$Cos_Delta_1 = 0
alsfrs$Cos_Delta_2 = 0
alsfrs$Cos_Delta_3 = 0
i=1
j=1
while(i<=m){
  
  t_min = alsfrs$Delta[j]
  t_max = alsfrs$Delta[j+numerosity[i]-1]
  
  alsfrs$Cos_Delta_1[j:(j+numerosity[i]-1)] = 
    cos(pi*(alsfrs$Delta[j:(j+numerosity[i]-1)] - t_min)/(t_max-t_min))
  
  alsfrs$Cos_Delta_2[j:(j+numerosity[i]-1)] = 
    cos(2*pi*(alsfrs$Delta[j:(j+numerosity[i]-1)] - t_min)/(t_max-t_min))
  
  alsfrs$Cos_Delta_3[j:(j+numerosity[i]-1)] = 
    cos(3*pi*(alsfrs$Delta[j:(j+numerosity[i]-1)] - t_min)/(t_max-t_min))
  
  j = j+numerosity[i]
  i = i+1
}

numerosity <- as.vector(table(alsfrs[,1])) ### number of obs for each patient
n1 = sum(numerosity[1:500])
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
while(i <= floor(m_train/2)){ ### 50% of training patient
  
  ### 50% of obs in train and 50% on test
  alsfrs_train$train_obs[(j + floor(num_train[i]/2)):(j+num_train[i]-1)] = 0
  j = j+num_train[i]
  i = i+1
}

alsfrs_train_patient_test_observation = alsfrs_train[which(alsfrs_train$train_obs == 0), c(-22)]
alsfrs_train = alsfrs_train[which(alsfrs_train$train_obs == 1), c(-22)]



########## fit model with brm #####################
fit1 <- brm(data = alsfrs_train,
            family = gaussian,
            ALSFRS_R_Total ~ 1 + Delta + Cos_Delta_1 + Cos_Delta_2 + Cos_Delta_3 +
              Onset_Site + Onset_Delta + 
              Age + FVC + SGPT + SGOT + RBC + WBC +
              Study_Treatment +  
              Delta*Onset_Site + Delta*Study_Treatment +
              (1 + Delta + Cos_Delta_1 + Cos_Delta_2 + Cos_Delta_3 || ID),
            
            prior = c(prior(normal(0, 10), class = Intercept),
                      prior(normal(0,10), class = b),
                      prior(inv_gamma(2, 1), class = sd),
                      prior(inv_gamma(2, 1), class = sigma)),    
            iter = 20000, warmup = 5000, chains = 1, cores = 4,
            #control = list(adapt_delta = .9, max_treedepth = 15),
            seed = 190831)
save(fit1, file = "brm_model_final.Rdata")


load("brm_model_final.Rdata")


summary(fit1)

x11()
plot(fit1)

marginal_effects(fit1)

### to recover the stancode
code = stancode(fit1)
code
standata = standata(fit1)

loo1 = loo(fit1) ### leave-one-out cross-validation 



################## model checking on training dataset ###################
x11()
par(mfrow = c(2,3))
for(n in (11*(1:6)+250)){
  new_train = alsfrs_train[which(alsfrs_train$ID == ID_train[n]), ]
  data = predict(fit1, newdata = new_train)
  
  plot(new_train$Delta, data[,1], ylim = c(0,48), type = "b", col = "red", pch = 20, xlab = "days", ylab = "ALSFRS")
  points(new_train$Delta, alsfrs_train[which(alsfrs_train$ID == ID_train[n]), c(3)], type = "b", col = "blue", pch = 17)
  points(new_train$Delta, data[,3], ylim = c(0,48), type = "b", col = "orange", pch = 20)
  points(new_train$Delta, data[,4], ylim = c(0,48), type = "b", col = "orange", pch = 20)
}



################# prediction on patients from train set but on new obs ########################
x11()
par(mfrow = c(2,3))
for(n in 11*(1:6)){
  new_test = alsfrs_train_patient_test_observation[which(alsfrs_train_patient_test_observation$ID == ID_train[n]), ]
  data = predict(fit1, newdata = new_test, allow_new_levels = TRUE)
  
  plot(new_test$Delta, data[,1], ylim = c(0,48), type = "b", col = "red", pch = 20, xlab = "days", ylab = "ALSFRS")
  points(new_test$Delta, alsfrs_train_patient_test_observation[which(alsfrs_train_patient_test_observation$ID == ID_train[n]), c(3)], type = "b", col = "blue", pch = 17)
  points(new_test$Delta, data[,3], ylim = c(0,48), type = "b", col = "orange", pch = 20)
  points(new_test$Delta, data[,4], ylim = c(0,48), type = "b", col = "orange", pch = 20)
}
###################################################################################



####################### prediction on new patients from test set ###################
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


### prediction for two patients in the training set 
### in blue observation in train, in green obs. in the train_patient_test_observation
x11()
par(mfrow = c(1,2))
new_test = rbind(alsfrs_train[10:15,], alsfrs_train_patient_test_observation[10:15,])
data = predict(fit1, newdata = new_test)
plot(new_test$Delta, data[,1], ylim = c(0,48), type = "b", col = "red", pch = 20, xlab = "days", ylab = "ALSFRS")
points(new_test$Delta[1:6], new_test$ALSFRS_R_Total[1:6], type = "b", col = "blue", pch = 17)
points(new_test$Delta[7:12], new_test$ALSFRS_R_Total[7:12], type = "b", col = "green", pch = 17)
points(new_test$Delta, data[,3], ylim = c(0,48), type = "b", col = "orange", pch = 20)
points(new_test$Delta, data[,4], ylim = c(0,48), type = "b", col = "orange", pch = 20)


new_test = rbind(alsfrs_train[24:31,], alsfrs_train_patient_test_observation[26:33,])
data = predict(fit1, newdata = new_test)
plot(new_test$Delta, data[,1], ylim = c(0,48), type = "b", col = "red", pch = 20, xlab = "days", ylab = "ALSFRS")
points(new_test$Delta[1:8], new_test$ALSFRS_R_Total[1:8], type = "b", col = "blue", pch = 17)
points(new_test$Delta[9:16], new_test$ALSFRS_R_Total[9:16], type = "b", col = "green", pch = 17)
points(new_test$Delta, data[,3], ylim = c(0,48), type = "b", col = "orange", pch = 20)
points(new_test$Delta, data[,4], ylim = c(0,48), type = "b", col = "orange", pch = 20)



######## errors ##########
data = predict(fit1, newdata = alsfrs_train)
err1_train = mean(abs(data[,1] - alsfrs_train$ALSFRS_R_Total))
err2_train = err1_train + mean(abs(data[,4] - data[,3]))
# 1.048792
# 9.024221

data = predict(fit1, newdata = alsfrs_train_patient_test_observation)
err1_train_patient_test_observation = mean(abs(data[,1] - alsfrs_train_patient_test_observation$ALSFRS_R_Total))
err2_train_patient_test_observation = err1_train_patient_test_observation + mean(abs(data[,4] - data[,3]))
# 3.582663
# 20.45746

data = predict(fit1, newdata = alsfrs_test, allow_new_levels = TRUE)
err1_test = mean(abs(data[,1] - alsfrs_test$ALSFRS_R_Total))
err2_test = err1_test + mean(abs(data[,4] - data[,3]))
# 5.977326
# 35.94814




data.out=as.matrix(fit1) 
data.out=data.frame(data.out)

chain = data.out[,1:16]
### density of the coefficients
x11()
par(mfrow = c(4,4))
for(i in 1:16){
  plot(density(chain[,i]),col="blue",lwd=2, xlab ="", ylab = "", main = names(chain)[i])
  abline(v=0, col = "black")
  abline(v=quantile(chain[,i],prob=c(0.025)),col="red",lty=2,lwd=2)
  abline(v=quantile(chain[,i],prob=c(0.5)),col="red",lty=1,lwd=2)
  abline(v=quantile(chain[,i],prob=c(0.975)),col="red",lty=2,lwd=2)
}


################## Anomaly detection using random effect #####################

mean_intercept = colMeans(data.out[,c(24:523)])
mean_slopes = colMeans(data.out[,c(524:1023)])
mean_cos1 = colMeans(data.out[,c(1024:1523)])
mean_cos2 = colMeans(data.out[,c(1524:2023)])
mean_cos3 = colMeans(data.out[,c(2024:2523)])

mean_random_effect = cbind(mean_intercept, mean_slopes, mean_cos1, mean_cos2, mean_cos3)

x11()
pairs(mean_random_effect)
summary(mean_random_effect)

### we can see which patients use more than the others the random effect parameter:
### this is a way to find some anomaly patient that cause this big eterogeneity in the dataset

### Let's plot the patient that achieve the min and the max for each random effect

### anomaly patient in mean_intercept
# for example
nn = c(which(mean_intercept == min(mean_intercept)), which(mean_intercept == max(mean_intercept)))
nn

x11()
par(mfrow = c(1,2))
for(n in nn){
  ID_patient = ID_train[n]
  plot(alsfrs_train[which(alsfrs_train == ID_patient),2], 
       alsfrs_train[which(alsfrs_train == ID_patient),3], ylim = c(0,48), type = "b",xlab = "days", ylab = "ALSFRS")
}

### anomaly patient in mean_slopes
# for example
nn = c(which(mean_slopes == min(mean_slopes)), which(mean_slopes == max(mean_slopes)))
nn

x11()
par(mfrow = c(1,2))
for(n in nn){
  ID_patient = ID_train[n]
  plot(alsfrs_train[which(alsfrs_train == ID_patient),2], 
       alsfrs_train[which(alsfrs_train == ID_patient),3], ylim = c(0,48), type = "b",xlab = "days", ylab = "ALSFRS")
}


### anomaly patient in mean_cos1
# for example
nn = c(which(mean_cos1 == min(mean_cos1)), which(mean_cos1 == max(mean_cos1)))
nn

x11()
par(mfrow = c(1,2))
for(n in nn){
  ID_patient = ID_train[n]
  plot(alsfrs_train[which(alsfrs_train == ID_patient),2], 
       alsfrs_train[which(alsfrs_train == ID_patient),3], ylim = c(0,48), type = "b",xlab = "days", ylab = "ALSFRS")
}

### anomaly patient in mean_cos2
# for example
nn = c(which(mean_cos2 == min(mean_cos2)), which(mean_cos2 == max(mean_cos2)))
nn

x11()
par(mfrow = c(1,2))
for(n in nn){
  ID_patient = ID_train[n]
  plot(alsfrs_train[which(alsfrs_train == ID_patient),2], 
       alsfrs_train[which(alsfrs_train == ID_patient),3], ylim = c(0,48), type = "b",xlab = "days", ylab = "ALSFRS")
}


### anomaly patient in mean_cos3
# for example
nn = c(which(mean_cos3 == min(mean_cos3)), which(mean_cos3 == max(mean_cos3)))
nn

x11()
par(mfrow = c(1,2))
for(n in nn){
  ID_patient = ID_train[n]
  plot(alsfrs_train[which(alsfrs_train == ID_patient),2], 
       alsfrs_train[which(alsfrs_train == ID_patient),3], ylim = c(0,48), type = "b",xlab = "days", ylab = "ALSFRS")
}

