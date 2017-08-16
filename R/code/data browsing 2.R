# This code is meant for data browsing and orientation only. The data concerned regards the Survey of Health, Ageing and 
# Retirement in Europe.

# convert SPSS file to R format
#install.packages("Hmisc")
library(foreign)
library(Hmisc)
#file.choose()
health_w1 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew1_rel6-0-0_ALL_datasets_spss\\sharew1_rel6-0-0_gv_health.sav",
                      to.data.frame=TRUE)
health_w2 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew2_rel6-0-0_ALL_datasets_spss\\sharew2_rel6-0-0_gv_health.sav",
                      to.data.frame=TRUE)
health_w4 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew4_rel6-0-0_ALL_datasets_spss\\sharew4_rel6-0-0_gv_health.sav",
                      to.data.frame=TRUE)
health_w5 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew5_rel6-0-0_ALL_datasets_spss\\sharew5_rel6-0-0_gv_health.sav",
                      to.data.frame=TRUE)
health_w6 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew6_rel6-0-0_ALL_datasets_spss\\sharew6_rel6-0-0_gv_health.sav",
                      to.data.frame=TRUE)
ph_w1 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew1_rel6-0-0_ALL_datasets_spss\\sharew1_rel6-0-0_ph.sav",
                  to.data.frame=TRUE)
ph_w2 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew2_rel6-0-0_ALL_datasets_spss\\sharew2_rel6-0-0_ph.sav",
                  to.data.frame=TRUE)
ph_w4 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew4_rel6-0-0_ALL_datasets_spss\\sharew4_rel6-0-0_ph.sav",
                  to.data.frame=TRUE)
ph_w5 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew5_rel6-0-0_ALL_datasets_spss\\sharew5_rel6-0-0_ph.sav",
                  to.data.frame=TRUE)
ph_w6 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew6_rel6-0-0_ALL_datasets_spss\\sharew6_rel6-0-0_ph.sav",
                  to.data.frame=TRUE)
cv_w1 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew1_rel6-0-0_ALL_datasets_spss\\sharew1_rel6-0-0_cv_r.sav",
                  to.data.frame=TRUE)
cv_w2 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew2_rel6-0-0_ALL_datasets_spss\\sharew2_rel6-0-0_cv_r.sav",
                  to.data.frame=TRUE)
cv_w3 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew4_rel6-0-0_ALL_datasets_spss\\sharew4_rel6-0-0_cv_r.sav",
                  to.data.frame=TRUE)
cv_w4 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew5_rel6-0-0_ALL_datasets_spss\\sharew5_rel6-0-0_cv_r.sav",
                  to.data.frame=TRUE)
cv_w5 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew6_rel6-0-0_ALL_datasets_spss\\sharew6_rel6-0-0_cv_r.sav",
                  to.data.frame=TRUE)
dn_w1 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew1_rel6-0-0_ALL_datasets_spss\\sharew1_rel6-0-0_dn.sav",
                  to.data.frame=TRUE)
read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew6_rel6-0-0_ALL_datasets_spss\\sharew6_rel6-0-0_cv_r.sav",
          to.data.frame=TRUE)
read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew6_rel6-0-0_ALL_datasets_spss\\sharew6_rel6-0-0_cv_r.sav",
          to.data.frame=TRUE)
read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew6_rel6-0-0_ALL_datasets_spss\\sharew6_rel6-0-0_cv_r.sav",
          to.data.frame=TRUE)

## try to merge the four waves on basis of ph, with attrition from first wave onwards
ph_mergeid <- list(mergeid1 = ph_w1[,"mergeid"], mergeid2 = ph_w2[,"mergeid"],
                   mergeid4 = ph_w4[,"mergeid"], mergeid5 = ph_w5[,"mergeid"], mergeid6 = ph_w6[,"mergeid"])
cv_mergeid <- list(mergeid1 = cv_w1[,"mergeid"], mergeid2 = cv_w2[,"mergeid"],
                   mergeid3 = cv_w3[,"mergeid"], mergeid4 = cv_w4[,"mergeid"], mergeid5 = cv_w5[,"mergeid"])
health_mergeid <- list(mergeid1 = health_w1[,"mergeid"], mergeid2 = health_w2[,"mergeid"],
                       mergeid4 = health_w4[,"mergeid"], mergeid5 = health_w5[,"mergeid"], mergeid6 = health_w6[,"mergeid"])
phvec <- unlist(ph_mergeid)
ph_unique <- sort(unique(phvec))
ph_waves <- sort(rep(ph_unique,5))
ph_retention <- matrix(0,length(ph_waves),1)
cv_retention <- matrix(0,length(ph_waves),1)
health_retention <- matrix(0,length(ph_waves),1)
index1 <- seq(1,(-4 + length(ph_waves)),by=5)
index2 <- seq(2,(-3 + length(ph_waves)),by=5)
index3 <- seq(3,(-2 + length(ph_waves)),by=5)
index4 <- seq(4,(-1 + length(ph_waves)),by=5)
index5 <- seq(5,(     length(ph_waves)),by=5)
for(i in 1:5){ 
  index <- eval(parse(text=paste("index",i,sep="")))
  ph_retention[index] <- match(ph_unique,ph_mergeid[[i]]) 
  cv_retention[index] <- match(ph_unique,cv_mergeid[[i]])
  health_retention[index] <- match(ph_unique,health_mergeid[[i]])
}
