alsfrs = read.csv("Alsfrs_RP.csv", header = TRUE)
alsfrs = alsfrs[,c(1,12,13,14,15,19)]

alsfrs = alsfrs[which(alsfrs$Repeated_Record == 1),]


alsfrs$counter = 0 ### number of obs for each ID
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
#select patients with at least 6 obs.
alsfrs = alsfrs[which(alsfrs$counter>5),]

alsfrs = alsfrs[order(alsfrs$subject_id, alsfrs$ALSFRS_Delta),]
ID = unique(alsfrs$subject_id)

### up to now, we have two different score: ALSFRS and ALSFRS_R
### most of the patients have only one of them, but we can compute ALSFRS_R from
### ALSFRS by adding two times the question 10 
alsfrs_complete = alsfrs

alsfrs_complete$ALSFRS_Total = alsfrs_complete$ALSFRS_Total + 2*alsfrs_complete$Q10_Respiratory
n = dim(alsfrs_complete)[1]

for(i in 1:n){
  if(is.na(alsfrs_complete$ALSFRS_R_Total[i])){
    alsfrs_complete$ALSFRS_R_Total[i] = alsfrs_complete$ALSFRS_Total[i]
  }
}
alsfrs_complete = alsfrs_complete[,c(1,3,5,7)]

alsfrs_complete = na.omit(alsfrs_complete)
alsfrs_normal = na.omit(alsfrs[,c(1,3,4,6)])
alsfrs_R = na.omit(alsfrs[,c(1,3,5,6)])

write.csv(alsfrs_normal, file = "alsfrs_normal_longitudinal.csv", row.names = FALSE)
write.csv(alsfrs_R, file = "alsfrs_R_longitudinal.csv", row.names = FALSE)
write.csv(alsfrs_complete, file = "alsfrs_R_complete_longitudinal.csv", row.names = FALSE)
