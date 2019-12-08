
fvc = read.csv("FVC_RP.csv", header = TRUE)
fvc = fvc[,c(1,2,8,14)] #uso subject_liter_trial_1
newdata = fvc
fvc = na.omit(fvc)

fvc = fvc[which(fvc$Forced_Vital_Capacity_Delta >= 0),]
fvc = fvc[which(fvc$Forced_Vital_Capacity_Delta < 400),]
fvc = fvc[which(fvc$Repeated_Record == 1),]


fvc$maxtime = 0
j=1
i = 1
while(i < dim(fvc)[1]){
  if(fvc$subject_id[i] == fvc$subject_id[j])
  {i = i+1}
  else
  {
    fvc$maxtime[j:(i-1)] = fvc$Forced_Vital_Capacity_Delta[i-1]
    j = i
    i = i+1
  }
}
#prendo gli ID con le valutazioni ad almeno 300 giorni
fvc = fvc[which(fvc$maxtime>300),]


fvc$counter = 0
j=1
i = 1
while(i < dim(fvc)[1]){
  if(fvc$subject_id[i] == fvc$subject_id[j])
  {i = i+1}
  else
  {
    fvc$counter[j:(i-1)] = i-j
    j = i
    i = i+1
  }
}
#prendo gli ID con almeno 6 osservazioni
fvc = fvc[which(fvc$counter>5),]


ID = unique(fvc$subject_id)

ID_0 = unique(fvc[which(fvc$Forced_Vital_Capacity_Delta ==0),]$subject_id) #tot = 2293 pazienti
ID_0 = sort(ID_0)

fvc = fvc[which(fvc$subject_id %in% ID_0),]

t = seq(0,360,30)
n = length(t)
n_id = length(ID_0)
data = matrix(1:(n*n_id), nrow = n_id*n, ncol = 1)
data = as.data.frame(data)

data$ID = 0
data$Delta = 0
data$fvc = 0

for(j in c(1:n_id)){
  for(i in c(1:n)){
    data$ID[(j-1)*n + i] = ID_0[j]
  }
}

for(j in c(1:n_id)){
  for(i in c(1:n)){
    data$Delta[(j-1)*n + i] = t[i]
  }
}

##################################

ID_0 = sort(ID_0)
for(i in 1:n_id){
  data$fvc[(i-1)*n +1] = fvc[which(fvc$subject_id == ID_0[i] & fvc$Forced_Vital_Capacity_Delta == 0),]$Subject_Liters_Trial_1
}

fvc$incremental = 0

i = 1

while(i < dim(fvc)[1]){
  if(fvc$subject_id[i] == fvc$subject_id[i+1])
  {
    fvc$incremental[i] = (fvc$Subject_Liters_Trial_1[i+1]-fvc$Subject_Liters_Trial_1[i])/(fvc$Forced_Vital_Capacity_Delta[i+1]-fvc$Forced_Vital_Capacity_Delta[i])
    i = i+1
  }
  else
  {
    i = i+1
  }
}


for(i in 1:n_id){
  obs = fvc[which(fvc$subject_id == ID_0[i]),c(2,3,7)] #we need score,delta and incremental
  j = 2
  m = dim(obs)[1]
  while(j <= n){
    t = 30*(j-1)
    l = 1
    for(k in 1:m){
      if(obs[k,]$Forced_Vital_Capacity_Delta <= t)
        l = k
    }
    data$fvc[(i-1)*n + j] = obs[l,]$Subject_Liters_Trial_1 + obs[l,]$incremental*(t - obs[l,]$Forced_Vital_Capacity_Delta)
    j = j+1
  }
}


data = data[,c(-1)]

write.csv(data, file = "fvc_longitudinal_1year.csv", row.names = FALSE)
