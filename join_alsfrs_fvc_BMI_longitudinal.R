
alsfrs = read.csv("alsfrs_longitudinal_1year.csv", header = TRUE)
fvc = read.csv("fvc_longitudinal_1year.csv", header = TRUE)
BMI = read.csv("BMI_longitudinal_1year.csv", header = TRUE)

library(dplyr)

longitudinal_data = inner_join(alsfrs,fvc)
longitudinal_data = inner_join(longitudinal_data, BMI)
# total of 1160 patients

write.csv(longitudinal_data, "alsfrs_fvc_BMI_longitudinal.csv", row.names = FALSE)
