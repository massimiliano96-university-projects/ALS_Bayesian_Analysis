library(dplyr)

################ FUNCTIONS ##################
interpol <- function(alsfrs, data) ### Linear Interpolation of variable v
{
  na_n = which(is.na(alsfrs$v)) # we need interpolation only when left_join(alsfrs, ...)
                                # gives us a NA for the variable at that specific delta
  
  min_delta = min(data$Delta)
  max_delta = max(data$Delta)
  
  for(i in na_n){
    if(alsfrs$Delta[i]<min_delta) # if we need to fill a NA with time before the first
    {                             # time evaluation of v, we just copy the first v we have
      alsfrs$v[i] = data$v[1]
    }
    else{
      if(alsfrs$Delta[i]<max_delta)
      { # linear interpolation between an evaluation before and after the time we need
        t = alsfrs$Delta[i]
        
        index_time_low = max(which(data$Delta <= t))
        index_time_up = min(which(data$Delta >= t))
        
        incremental = (data$v[index_time_up] - data$v[index_time_low])/
          (data$Delta[index_time_up] - data$Delta[index_time_low])
        
        alsfrs$v[i] = data$v[index_time_low] + incremental*(t-data$Delta[index_time_low])
      }
      else # if we need to fill a NA with time after the last
      {    # time evaluation of v, we just copy the last v we have
        alsfrs$v[i] = data$v[length(data$v)]
      }
    }
  }
  return(alsfrs$v)
}


join_alsfrs_interpolation <- function( alsfrs, data )
{ # Join alsfrs with a new longitudinal variable v after an interpolation to match the Deltas
  ###
  ### data must be of the form: data[,1] = ID, data[,2] = delta, data[,3] = variable
  ###
  old_name = names(data)[3]   # replace variable name with a generic "v"
  names(data)[3] = "v"
  
  new_alsfrs = alsfrs
  new_alsfrs$ID_present = 0
  
  # to see for which patients we can make an interpolation (only if we have an evaluation in data)
  new_alsfrs[which(new_alsfrs$ID %in% unique(data$ID)),]$ID_present = 1
  
  new_alsfrs = left_join(new_alsfrs, data) # left join with on the left alsfrs
  # we have only the rows of alsfrs, and if we dont have v, we get a NA
  
  n = dim(new_alsfrs)[2]
  m = length(unique(new_alsfrs$ID))
  nn = as.vector(table(new_alsfrs[,1]))
  
  i = 1
  j = 1
  while(i<=m) ### loop on all the m patients
  {
    if(new_alsfrs$ID_present[j] == 1) # if the ID is present in data
    {
      if(sum(is.na(new_alsfrs[j:(j+nn[i]-1),]$v)) > 0) # if we have at least one NA for that patient
      {
        new_ID = new_alsfrs$ID[j]
        new_alsfrs$v[j:(j+nn[i]-1)] = 
          interpol(new_alsfrs[j:(j+nn[i]-1),], data[which(data$ID == new_ID),] )
        # call to function interpol with input the rows and the data of that patient  
      }
    }
    
    j = j+nn[i]
    i = i+1  
  }
  
  names(new_alsfrs)[n] = old_name # put the right name of variable v
  return(new_alsfrs[,-c(n-1)]) # cancel the column ID_present
}






###we're gonna use ALSFRS-R
alsfrs <- read.csv("alsfrs_R_complete_longitudinal.csv", header = TRUE)
names(alsfrs)[1] = "ID"
names(alsfrs)[2] = "Delta"
alsfrs = alsfrs[,c(-4)]  ### delete counter


### Bulbar or Limb
onset_type = read.csv("alstype_char_final.csv", header = TRUE) 
names(onset_type)[1] = "ID"

alsfrs = inner_join(alsfrs, onset_type)


### Onset_Delta
onset_delta = read.csv("ALSHistory_RP_Site.csv", header = TRUE)
names(onset_delta)[1] = "ID"
onset_delta = onset_delta[,c(-2,-5)]

alsfrs = left_join(alsfrs, onset_delta)


### Age
age = read.csv("Demographics_RP.csv", header = TRUE)
age = age[,c(1,3,4)]
names(age)[1] = "ID"

# if we don't have age we can compute it thanks to "Date_of_Birth" column
age[which(is.na(age$Age)),2] = age[which(is.na(age$Age)),3]/(-365.25)
age = age[,-3]

alsfrs = left_join(alsfrs, age)


### FVC
fvc = read.csv("FVC_RP.csv", header = TRUE)
fvc = fvc[which(fvc$Repeated_Record == 1),c(1,8,2)]
names(fvc)[1] = "ID"
names(fvc)[2] = "Delta"
names(fvc)[3] = "FVC"

alsfrs = join_alsfrs_interpolation(alsfrs, fvc) ### join with interpolation! FVC is longitudinal


### BMI
bmi = read.csv("WeightHeight_RP.csv", header = TRUE)

n_bmi = dim(bmi)[1] ### Compute the BMI
for(i in 1:n_bmi){
  if(is.na(bmi$BMI[i])){
    if(is.na(bmi$BMI_Kilo[i])){
      if(is.na(bmi$BMI_Pounds[i])){
        bmi$BMI[i] = bmi$Final_Weight[i]/((bmi$Final_Height[i]/100)^2)
      }
      else
        bmi$BMI[i] = bmi$BMI_Pounds[i]*2.205
    }
    else
      bmi$BMI[i] = bmi$BMI_Kilo[i]
  }
}
bmi = bmi[,c(1,2,19)]
names(bmi)[1] = "ID"
names(bmi)[2] = "Delta"
bmi = na.omit(bmi)

bmi$BMI_0 = 0 ### we want to use the BMI at time zero and the % difference in time from that
bmi$BMI_diff = 0

m = length(unique(bmi$ID))
numerosity <- as.vector(table(bmi[,1])) ### number of obs for each patient

i=1
j=1
while(i<=m){
  bmi$BMI_0[j:(j+numerosity[i]-1)] = bmi$BMI[j]
  bmi$BMI_diff[j:(j+numerosity[i]-1)] = 
    (bmi$BMI[j:(j+numerosity[i]-1)] - bmi$BMI_0[j:(j+numerosity[i]-1)])/bmi$BMI_0[j:(j+numerosity[i]-1)]
  j = j+numerosity[i]
  i = i+1
}

bmi = bmi[,c(-3)]


m = length(unique(bmi$ID))
numerosity <- as.vector(table(bmi[,1])) ### number of obs for each patient

n_bmi = dim(bmi)[1]
bmi$rep = 1
for(i in 1:(n_bmi-1)){
  if(bmi$ID[i] == bmi$ID[i+1])
    if(bmi$Delta[i] == bmi$Delta[i+1])
      bmi$rep[i+1] = 2
}
bmi = bmi[which(bmi$rep == 1),c(-5)] ### take only one repetition

alsfrs = join_alsfrs_interpolation(alsfrs, bmi[,c(1,2,3)]) #joi with interpolation
alsfrs = join_alsfrs_interpolation(alsfrs, bmi[,c(1,2,4)])


### lab categorical variables
labcat = read.csv("LabsCategorical_RP.csv", header = TRUE)
labcat = labcat[which(labcat$Duplicated_Date == 1),]
levels(labcat$Test_name)
numerosity <- as.vector(table(labcat[,1])) ### number of obs for each patient
which(numerosity>19900) ### select the variables more used
#34 49 58
#"URINE APPEARANCE"
#"URINE GLUCOSE"
#"URINE PROTEIN"
labcat1 = labcat[which(labcat$Test_name == "URINE APPEARANCE"),c(2,3,11)]
numerosity <- as.vector(table(labcat1[,3]))
labcat1$RevisedTest_result = as.character.factor(labcat1$RevisedTest_result)
labcat1$RevisedTest_result = as.factor(labcat1$RevisedTest_result)
levels(labcat1$RevisedTest_result)

names(labcat1)[1] = "ID"
names(labcat1)[2] = "Delta"
names(labcat1)[3] = "URINE APPEARANCE"


labcat2 = labcat[which(labcat$Test_name == "URINE GLUCOSE"),c(2,3,11)]
numerosity <- as.vector(table(labcat2[,3]))
labcat2$RevisedTest_result = as.character.factor(labcat2$RevisedTest_result)
labcat2$RevisedTest_result = as.factor(labcat2$RevisedTest_result)
levels(labcat2$RevisedTest_result)

names(labcat2)[1] = "ID"
names(labcat2)[2] = "Delta"
names(labcat2)[3] = "URINE GLUCOSE"


labcat3 = labcat[which(labcat$Test_name == "URINE PROTEIN"),c(2,3,11)]
numerosity <- as.vector(table(labcat3[,3]))
labcat3$RevisedTest_result = as.character.factor(labcat3$RevisedTest_result)
labcat3$RevisedTest_result = as.factor(labcat3$RevisedTest_result)
levels(labcat3$RevisedTest_result)

names(labcat3)[1] = "ID"
names(labcat3)[2] = "Delta"
names(labcat3)[3] = "URINE PROTEIN"

labcat3$rep = 1
n_labcat3 = dim(labcat3)[1]
for(i in 1:(n_labcat3-1)){
  if(labcat3$ID[i] == labcat3$ID[i+1])
    if(labcat3$Delta[i] == labcat3$Delta[i+1])
      labcat3$rep[i+1] = 2
}
labcat3 = labcat3[which(labcat3$rep == 1),c(-4)]


### lab numerical variables
labnum = read.csv("LabsNumerical_RP.csv", header = TRUE)
labnum = labnum[which(labnum$Duplicated_Date == 1),]
levels(labnum$Test_name)
numerosity <- as.vector(table(labnum[,1])) ### number of obs for each patient
which(numerosity>60000) #cambiare!
#18  22  31  56  59  60 105 143
#"ALT(SGPT)"
#"AST(SGOT)"
#"BILIRUBIN (TOTAL)"
#"GLUCOSE"
#"HEMATOCRIT"
#"HEMOGLOBIN"
#"RED BLOOD CELLS (RBC)"
#"WHITE BLOOD CELL (WBC)" 


### + 135 "URINE PH"  +   124 "URIC ACID"
labnum1 = labnum[which(labnum$Test_name == "ALT(SGPT)"), c(2,3,11)]
labnum2 = labnum[which(labnum$Test_name == "AST(SGOT)"), c(2,3,11)]
labnum3 = labnum[which(labnum$Test_name == "BILIRUBIN (TOTAL)"), c(2,3,11)]
labnum4 = labnum[which(labnum$Test_name == "GLUCOSE"), c(2,3,11)]
labnum5 = labnum[which(labnum$Test_name == "HEMATOCRIT"), c(2,3,11)]
labnum6 = labnum[which(labnum$Test_name == "HEMOGLOBIN"), c(2,3,11)]
labnum7 = labnum[which(labnum$Test_name == "RED BLOOD CELLS (RBC)"), c(2,3,11)]
labnum8 = labnum[which(labnum$Test_name == "WHITE BLOOD CELL (WBC)"), c(2,3,11)]
labnum9 = labnum[which(labnum$Test_name == "URINE PH"), c(2,3,11)]
labnum10 = labnum[which(labnum$Test_name == "URIC ACID"), c(2,3,11)]

names(labnum1)[1] = "ID"
names(labnum1)[2] = "Delta"
names(labnum1)[3] = "SGPT"

names(labnum2)[1] = "ID"
names(labnum2)[2] = "Delta"
names(labnum2)[3] = "SGOT"

names(labnum3)[1] = "ID"
names(labnum3)[2] = "Delta"
names(labnum3)[3] = "BILIRUBIN"

names(labnum4)[1] = "ID"
names(labnum4)[2] = "Delta"
names(labnum4)[3] = "GLUCOSE"

names(labnum5)[1] = "ID"
names(labnum5)[2] = "Delta"
names(labnum5)[3] = "HEMATOCRIT"

names(labnum6)[1] = "ID"
names(labnum6)[2] = "Delta"
names(labnum6)[3] = "HEMOGLOBIN"

names(labnum7)[1] = "ID"
names(labnum7)[2] = "Delta"
names(labnum7)[3] = "RBC"

names(labnum8)[1] = "ID"
names(labnum8)[2] = "Delta"
names(labnum8)[3] = "WBC"

names(labnum9)[1] = "ID"
names(labnum9)[2] = "Delta"
names(labnum9)[3] = "URINE PH"

names(labnum10)[1] = "ID"
names(labnum10)[2] = "Delta"
names(labnum10)[3] = "URIC ACID"

alsfrs = join_alsfrs_interpolation(alsfrs, labnum1)
alsfrs = join_alsfrs_interpolation(alsfrs, labnum2)
alsfrs = join_alsfrs_interpolation(alsfrs, labnum3)
alsfrs = join_alsfrs_interpolation(alsfrs, labnum4)
alsfrs = join_alsfrs_interpolation(alsfrs, labnum5)
alsfrs = join_alsfrs_interpolation(alsfrs, labnum6)
alsfrs = join_alsfrs_interpolation(alsfrs, labnum7)
alsfrs = join_alsfrs_interpolation(alsfrs, labnum8)
alsfrs = join_alsfrs_interpolation(alsfrs, labnum9)
alsfrs = join_alsfrs_interpolation(alsfrs, labnum10)

alsfrs = left_join(alsfrs, labcat1)
alsfrs = left_join(alsfrs, labcat2)
alsfrs = left_join(alsfrs, labcat3)


### RILUZOLE, SVC, Treatment (Active or Placebo), Vitals
riluzole = read.csv("Riluzole_RP_Added.csv", header = TRUE)
svc = read.csv("SVC_RP.csv", header = TRUE)
treatment = read.csv("Treatment.csv", header = TRUE)
vitals = read.csv("Vitals_RP.csv", header = TRUE)

riluzole_id = as.data.frame(unique(riluzole$subject_id))

riluzole_id$RILUZOLE = 1
names(riluzole_id)[1] = "ID"

alsfrs = left_join(alsfrs, riluzole_id)
alsfrs[which(is.na(alsfrs$RILUZOLE)),]$RILUZOLE = 0


svc = svc[,c(1,9,2)]
names(svc)[1] = "ID"
names(svc)[2] = "Delta"
names(svc)[3] = "SVC"
svc = na.omit(svc)


alsfrs = join_alsfrs_interpolation(alsfrs, svc)


treatment = treatment[,c(1,2)]
names(treatment)[1] = "ID"
names(treatment)[2] = "Study Treatment"

alsfrs = left_join(alsfrs, treatment)


### adjust the names of columns with " "
names(alsfrs)[19] = "URINE_PH"
names(alsfrs)[20] = "URIC_ACID"
names(alsfrs)[21] = "URINE_APPEARANCE"
names(alsfrs)[22] = "URINE_GLUCOSE"
names(alsfrs)[23] = "URINE_PROTEIN"
names(alsfrs)[26] = "Study_Treatment"


### Sex
sex = read.csv("Demographics_RP.csv", header = TRUE)
sex = sex[,c(1,5)]
sex = sex[which(sex$Sex != ""),]
names(sex)[1] = "ID"
alsfrs = left_join(alsfrs,sex)



##### choose the variables with less than 35% of NA -> Threshold
###
perc = rep(0,27)
n_dim = dim(alsfrs)[1]
for(i in 1:27){
  perc[i] = sum(is.na(alsfrs[,c(i)]))/n_dim
}

which(perc<=0.35)


### select only these columns #####
alsfrs = alsfrs[, which(perc<=0.35)]
alsfrs_riempito = alsfrs




####### how can we fill these NA? -> MICE + Other tecniques

######### How can we simulate fixed covariates?
#### 1) BMI_0
boxplot(alsfrs[which(alsfrs$Sex == "Male"),]$BMI_0)
mean(na.omit(alsfrs[which(alsfrs$Sex == "Male"),]$BMI_0))
median(na.omit(alsfrs[which(alsfrs$Sex == "Male"),]$BMI_0))
### for the male, we can set 25.8 as BMI_0

boxplot(alsfrs[which(alsfrs$Sex == "Female"),]$BMI_0)
mean(na.omit(alsfrs[which(alsfrs$Sex == "Female"),]$BMI_0))
median(na.omit(alsfrs[which(alsfrs$Sex == "Female"),]$BMI_0))
### for the female, we can set 24.4 as BMI_0

alsfrs_riempito[which(is.na(alsfrs_riempito$BMI_0) & alsfrs_riempito$Sex == "Male"),]$BMI_0 = 25.8
alsfrs_riempito[which(is.na(alsfrs_riempito$BMI_0) & alsfrs_riempito$Sex == "Female"),]$BMI_0 = 24.4


#### 2) Onset_Delta
summary(alsfrs$Onset_Delta) ### fill with the median
alsfrs_riempito[which(is.na(alsfrs_riempito$Onset_Delta)),]$Onset_Delta = -567.0


### now, for the time-variant covariates -> MICE
library(mice)

alsfrs_riempito <- mice(alsfrs_riempito, m = 10)

alsfrs_finale = alsfrs_riempito
alsfrs_finale = complete(alsfrs_finale)

write.csv(alsfrs, "alsfrs_non_riempito.csv", row.names = FALSE)
write.csv(alsfrs_finale, "alsfrs_riempito.csv", row.names = FALSE)

