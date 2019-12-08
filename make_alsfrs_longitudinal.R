
alsfrs = read.csv("Alsfrs_RP.csv", header = TRUE)
alsfrs = alsfrs[,c(1,13,14,19)]
newdata = alsfrs
alsfrs = na.omit(alsfrs)

alsfrs = alsfrs[which(alsfrs$ALSFRS_Delta >= 0),]
alsfrs = alsfrs[which(alsfrs$ALSFRS_Delta < 400),]
alsfrs = alsfrs[which(alsfrs$Repeated_Record == 1),]


alsfrs$maxtime = 0
j=1
i = 1
while(i < dim(alsfrs)[1]){
  if(alsfrs$subject_id[i] == alsfrs$subject_id[j])
  {i = i+1}
  else
  {
    alsfrs$maxtime[j:(i-1)] = alsfrs$ALSFRS_Delta[i-1]
    j = i
    i = i+1
  }
}
#prendo gli ID con le valutazioni ad almeno 300 giorni
alsfrs = alsfrs[which(alsfrs$maxtime>300),]


alsfrs$counter = 0
j=1
i = 1
while(i < dim(alsfrs)[1]){
  if(alsfrs$subject_id[i] == alsfrs$subject_id[j])
  {i = i+1}
  else
  {
    alsfrs$counter[j:(i-1)] = i-j
    j = i
    i = i+1
  }
}
#prendo gli ID con almeno 6 osservazioni
alsfrs = alsfrs[which(alsfrs$counter>5),]


ID = unique(alsfrs$subject_id)

ID_0 = unique(alsfrs[which(alsfrs$ALSFRS_Delta ==0),]$subject_id) #tot = 1656 pazienti
ID_0 = sort(ID_0)

alsfrs = alsfrs[which(alsfrs$subject_id %in% ID_0),]

t = seq(0,360,30)
n = length(t)
n_id = length(ID_0)
data = matrix(1:(n*n_id), nrow = n_id*n, ncol = 1)
data = as.data.frame(data)

data$ID = 0
data$Delta = 0
data$ALSFRS = 0

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
  data$ALSFRS[(i-1)*n +1] = alsfrs[which(alsfrs$subject_id == ID_0[i] & alsfrs$ALSFRS_Delta == 0),]$ALSFRS_Total
}

alsfrs$incremental = 0

i = 1

while(i < dim(alsfrs)[1]){
  if(alsfrs$subject_id[i] == alsfrs$subject_id[i+1])
  {
    alsfrs$incremental[i] = (alsfrs$ALSFRS_Total[i+1]-alsfrs$ALSFRS_Total[i])/(alsfrs$ALSFRS_Delta[i+1]-alsfrs$ALSFRS_Delta[i])
    i = i+1
  }
  else
  {
    i = i+1
  }
}


for(i in 1:n_id){
  obs = alsfrs[which(alsfrs$subject_id == ID_0[i]),c(2,3,7)] #we need score, delta and incremental
  j = 2
  m = dim(obs)[1]
  while(j <= n){
    t = 30*(j-1)
    l = 1
    for(k in 1:m){
      if(obs[k,]$ALSFRS_Delta <= t)
        l = k
    }
    data$ALSFRS[(i-1)*n + j] = obs[l,]$ALSFRS_Total + obs[l,]$incremental*(t - obs[l,]$ALSFRS_Delta)
    j = j+1
  }
}


data = data[,c(-1)]

write.csv(data, file = "alsfrs_longitudinal_1year.csv", row.names = FALSE, col.names = TRUE)
