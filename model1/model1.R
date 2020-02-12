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

alsfrs = alsfrs[which(alsfrs$counter >= 5),] ## we use only patients with more then 5 observations

alsfrs = alsfrs[,c(-22)] ### delete counter, change number

numerosity <- as.vector(table(alsfrs[,1])) ### number of obs for each patient
n1 = sum(numerosity[1:500]) ## number of rows in order to take 500 patients
alsfrs_train = alsfrs[1:n1,] ### 500 patients

ID_train = unique(alsfrs_train$ID) ## ID patients in the training set
n_train = dim(alsfrs_train)[1] ## number of rows
m_train = length(ID_train) ## number of patients
numerosity_train <- as.vector(table(alsfrs_train[,1])) ### number of obs for each patient

##### Test set
n2 = sum(numerosity[501:800])
alsfrs_test = alsfrs[(n1+1):(n1+n2),] ### 300 patients

ID_test = unique(alsfrs_test$ID)
n_test = dim(alsfrs_test)[1]
m_test = length(ID_test)
numerosity_test <- as.vector(table(alsfrs_test[,1])) ### number of obs for each patient


############# split some of the train patients in train and set observations #######
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
            ALSFRS_R_Total ~ 1 + Delta + (1 + Delta || ID),
            
            prior = c(prior(normal(0, 10), class = Intercept),
                      prior(normal(0, 10), class = b),
                      prior(inv_gamma(2, 1), class = sd),
                      prior(inv_gamma(2, 1), class = sigma)),    
            iter = 20000, warmup = 5000, chains = 1, cores = 4,
            #control = list(adapt_delta = .975, max_treedepth = 20),
            seed = 190831)
save(fit1, file = "brm_model_Delta.Rdata")


load("brm_model_Delta.Rdata")

####### Summary
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

savePlot(filename = paste("plot_model_checking_train.","jpeg"), device = dev.cur())


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
savePlot(filename = paste("plot_model_checking_test.","jpeg"), device = dev.cur())

###################################################################################



######## errors ##########
data = predict(fit1, newdata = alsfrs_train)
mean(abs(data[,1] - alsfrs_train$ALSFRS_R_Total))

data = predict(fit1, newdata = alsfrs_train_patient_test_observation)
mean(abs(data[,1] - alsfrs_train_patient_test_observation$ALSFRS_R_Total))

data = predict(fit1, newdata = alsfrs_test, allow_new_levels = TRUE)
mean(abs(data[,1] - alsfrs_test$ALSFRS_R_Total))

x11()
hist(abs(data[,1] - alsfrs_test$ALSFRS_R_Total), main = "Histogram of error on test set" ,col = "blue", xlab = "abs(ALSFRS_fitted - ALSFRS_real)", border = "white", xlim = c(0,35))
savePlot(filename = paste("plot_model1_hist_err.","jpeg"), device = dev.cur())



#### Analysis of betas

data.out=as.matrix(fit1) # trasform the mcmc.list into a matrix,
data.out=data.frame(data.out)     # or, better, into a dataframe (easiest to handle in R)
attach(data.out)
n.chain=dim(data.out)[1] 
n.chain # this is the final sample size 

chain <- data.out[,2]
x11()
layout(matrix(c(1,2,3,3),2,2,byrow=T))
plot(chain,type="l",main="Trace plot of beta2")
acf(chain,lwd=3,col="red3",main="autocorrelation of beta2")
hist(chain,nclass="fd",freq=F,main="Posterior of beta2",col="gray") 
lines(density(chain),col="blue",lwd=2)
quantile(chain,prob=c(0.025,0.5,0.975))


abline(v=quantile(chain,prob=c(0.025)),col="red",lty=2,lwd=2)
abline(v=quantile(chain,prob=c(0.5)),col="red",lty=1,lwd=2)
abline(v=quantile(chain,prob=c(0.975)),col="red",lty=2,lwd=2)
legend("topright",legend=c("posterior median", "95% Credible bounds","kernel density smoother"),lwd=c(2,2,2), col=c("red","red","blue"),lty=c(1,2,1))

savePlot(filename = paste("plot_model1_beta2.","jpeg"), device = dev.cur())




################## Outlier detenction in random effects #####################
mean_intercept = colMeans(data.out[,c(7:506)])
mean_slopes = colMeans(data.out[,c(507:1006)])

x11()
plot(mean_intercept, mean_slopes, col = "black", pch = 20)
points(mean(mean_intercept), mean(mean_slopes), col = "red", lwd = 5, pch = 10)
savePlot(filename = paste("random_effect_outliers.","jpeg"), device = dev.cur())

outliers <- data.frame(unique(alsfrs_train$ID))
outliers[,2] <- mean_intercept
outliers[,3] <- mean_slopes
colnames(outliers)[1] <- c("ID")
colnames(outliers)[2] <- c("mean_intercept")
colnames(outliers)[3] <- c("mean_slope")
outliers1 <- outliers[outliers$mean_intercept < -17,]
outliers2 <- outliers[abs(outliers$mean_slope) > 0.07,]
outliers3 <- outliers[ outliers$mean_intercept >11,]

outliers <- rbind(outliers1,outliers2, outliers3)

outliers_data <- alsfrs_train[alsfrs_train$ID %in% outliers$ID,]

x11()
par(mfrow = c(3,5))
for(n in c(1:length(outliers$ID))){
  new_test = alsfrs_train[which(alsfrs_train$ID == outliers$ID[n]), ]
  data = predict(fit1, newdata = new_test, allow_new_levels = TRUE)
  
  plot(new_test$Delta, data[,1], ylim = c(0,48), type = "b", col = "red", pch = 20, xlab = "days", ylab = "ALSFRS")
  points(new_test$Delta, alsfrs_train[which(alsfrs_train$ID == outliers$ID[n]), c(3)], type = "b", col = "blue", pch = 17)
  points(new_test$Delta, data[,3], ylim = c(0,48), type = "b", col = "orange", pch = 20)
  points(new_test$Delta, data[,4], ylim = c(0,48), type = "b", col = "orange", pch = 20)

#### fixed part
  new_test = alsfrs_train[alsfrs_train$ID == outliers$ID[n],]
  new_test$ID = 1 ## new ID, so we treat it as a new patient!
  data = predict( fit1, newdata = new_test, allow_new_levels = TRUE)
  points(new_test$Delta, data[,1], ylim = c(0,48), type = "b", col = "grey", pch = 17)
  points(new_test$Delta, data[,3], ylim = c(0,48), type = "b", col = "grey", pch = 20)
  points(new_test$Delta, data[,4], ylim = c(0,48), type = "b", col = "grey", pch = 20)
  
}
savePlot(filename = paste("outliers_score.","jpeg"), device = dev.cur())



