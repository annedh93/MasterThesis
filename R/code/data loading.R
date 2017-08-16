# This code is meant for data browsing and orientation only. The data concerned regards the Survey of Health, Ageing and 
# Retirement in Europe.

# convert SPSS file to R format
install.packages("Hmisc")
install.packages("memisc")
library(foreign)
library(Hmisc)
library(memisc)
#file.choose()

# generalized variable imputations
imp_1 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew1_rel6-0-0_ALL_datasets_spss\\sharew1_rel6-0-0_gv_imputations.sav",
                   to.data.frame=TRUE)
imp_w1 <- imp_1[,c("mergeid","chronic","thinc","yedu","sphus","mstat","nchild","gali","adl","iadl","maxgrip")]
rm(imp_1)
imp_2 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew2_rel6-0-0_ALL_datasets_spss\\sharew2_rel6-0-0_gv_imputations.sav",
                   to.data.frame=TRUE)
imp_w2 <- imp_2[,c("mergeid","chronic","thinc","yedu","sphus","mstat","nchild","gali","adl","iadl","maxgrip")]
rm(imp_2)
imp_3 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew4_rel6-0-0_ALL_datasets_spss\\sharew4_rel6-0-0_gv_imputations.sav",
                   to.data.frame=TRUE)
imp_w3 <- imp_3[,c("mergeid","chronic","thinc","yedu","sphus","mstat","nchild","gali","adl","iadl","maxgrip")]
rm(imp_3)
imp_4 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew5_rel6-0-0_ALL_datasets_spss\\sharew5_rel6-0-0_gv_imputations.sav",
                   to.data.frame=TRUE)
imp_w4 <- imp_4[,c("mergeid","chronic","thinc","yedu","sphus","mstat","nchild","gali","adl","iadl","maxgrip")]
rm(imp_4)
imp_5 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew6_rel6-0-0_ALL_datasets_spss\\sharew6_rel6-0-0_gv_imputations.sav",
                   to.data.frame=TRUE)
imp_w5 <- imp_5[,c("mergeid","chronic","thinc","yedu","sphus","mstat","nchild","gali","adl","iadl","maxgrip")]
rm(imp_5)
# coverscreen on individual level
cv_1 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew1_rel6-0-0_ALL_datasets_spss\\sharew1_rel6-0-0_cv_r.sav",
                  to.data.frame=TRUE)
cv_2 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew2_rel6-0-0_ALL_datasets_spss\\sharew2_rel6-0-0_cv_r.sav",
                  to.data.frame=TRUE)
cv_3 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew4_rel6-0-0_ALL_datasets_spss\\sharew4_rel6-0-0_cv_r.sav",
                  to.data.frame=TRUE)
cv_4 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew5_rel6-0-0_ALL_datasets_spss\\sharew5_rel6-0-0_cv_r.sav",
                  to.data.frame=TRUE)
cv_5 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew6_rel6-0-0_ALL_datasets_spss\\sharew6_rel6-0-0_cv_r.sav",
                  to.data.frame=TRUE)
# demographics and networks
dn_w1 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew1_rel6-0-0_ALL_datasets_spss\\sharew1_rel6-0-0_dn.sav",
                  to.data.frame=TRUE)
dn_w2 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew2_rel6-0-0_ALL_datasets_spss\\sharew2_rel6-0-0_dn.sav",
                  to.data.frame=TRUE)
dn_w3 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew4_rel6-0-0_ALL_datasets_spss\\sharew4_rel6-0-0_dn.sav",
                  to.data.frame=TRUE)
dn_w4 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew5_rel6-0-0_ALL_datasets_spss\\sharew5_rel6-0-0_dn.sav",
                  to.data.frame=TRUE)
dn_w5 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew6_rel6-0-0_ALL_datasets_spss\\sharew6_rel6-0-0_dn.sav",
                  to.data.frame=TRUE)
# children 
ch_w1 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew1_rel6-0-0_ALL_datasets_spss\\sharew1_rel6-0-0_ch.sav",
                  to.data.frame=TRUE)
ch_w2 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew2_rel6-0-0_ALL_datasets_spss\\sharew2_rel6-0-0_ch.sav",
                  to.data.frame=TRUE)
ch_w3 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew4_rel6-0-0_ALL_datasets_spss\\sharew4_rel6-0-0_ch.sav",
                  to.data.frame=TRUE)
ch_w4 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew5_rel6-0-0_ALL_datasets_spss\\sharew5_rel6-0-0_ch.sav",
                  to.data.frame=TRUE)
ch_w5 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew6_rel6-0-0_ALL_datasets_spss\\sharew6_rel6-0-0_ch.sav",
                  to.data.frame=TRUE)
# physicial health
ph_w1 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew1_rel6-0-0_ALL_datasets_spss\\sharew1_rel6-0-0_ph.sav",
                  to.data.frame=TRUE)
ph_w2 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew2_rel6-0-0_ALL_datasets_spss\\sharew2_rel6-0-0_ph.sav",
                  to.data.frame=TRUE)
ph_w3 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew4_rel6-0-0_ALL_datasets_spss\\sharew4_rel6-0-0_ph.sav",
                  to.data.frame=TRUE)
ph_w4 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew5_rel6-0-0_ALL_datasets_spss\\sharew5_rel6-0-0_ph.sav",
                  to.data.frame=TRUE)
ph_w5 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew6_rel6-0-0_ALL_datasets_spss\\sharew6_rel6-0-0_ph.sav",
                  to.data.frame=TRUE)
# employment and pensions
ep_w1 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew1_rel6-0-0_ALL_datasets_spss\\sharew1_rel6-0-0_ep.sav",
                  to.data.frame=TRUE)
ep_w2 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew2_rel6-0-0_ALL_datasets_spss\\sharew2_rel6-0-0_ep.sav",
                  to.data.frame=TRUE)
ep_w3 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew4_rel6-0-0_ALL_datasets_spss\\sharew4_rel6-0-0_ep.sav",
                  to.data.frame=TRUE)
ep_w4 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew5_rel6-0-0_ALL_datasets_spss\\sharew5_rel6-0-0_ep.sav",
                  to.data.frame=TRUE)
ep_w5 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew6_rel6-0-0_ALL_datasets_spss\\sharew6_rel6-0-0_ep.sav",
                  to.data.frame=TRUE)
# activities
ac_w1 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew1_rel6-0-0_ALL_datasets_spss\\sharew1_rel6-0-0_ac.sav",
                  to.data.frame=TRUE)
ac_w2 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew2_rel6-0-0_ALL_datasets_spss\\sharew2_rel6-0-0_ac.sav",
                  to.data.frame=TRUE)
ac_w3 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew4_rel6-0-0_ALL_datasets_spss\\sharew4_rel6-0-0_ac.sav",
                  to.data.frame=TRUE)
ac_w4 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew5_rel6-0-0_ALL_datasets_spss\\sharew5_rel6-0-0_ac.sav",
                  to.data.frame=TRUE)
ac_w5 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew6_rel6-0-0_ALL_datasets_spss\\sharew6_rel6-0-0_ac.sav",
                  to.data.frame=TRUE)
# housing
ho_w1 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew1_rel6-0-0_ALL_datasets_spss\\sharew1_rel6-0-0_ho.sav",
                  to.data.frame=TRUE)
ho_w2 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew2_rel6-0-0_ALL_datasets_spss\\sharew2_rel6-0-0_ho.sav",
                  to.data.frame=TRUE)
ho_w3 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew4_rel6-0-0_ALL_datasets_spss\\sharew4_rel6-0-0_ho.sav",
                  to.data.frame=TRUE)
ho_w4 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew5_rel6-0-0_ALL_datasets_spss\\sharew5_rel6-0-0_ho.sav",
                  to.data.frame=TRUE)
ho_w5 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew6_rel6-0-0_ALL_datasets_spss\\sharew6_rel6-0-0_ho.sav",
                  to.data.frame=TRUE)
# generalized variable health
health_w1 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew1_rel6-0-0_ALL_datasets_spss\\sharew1_rel6-0-0_gv_health.sav",
                      to.data.frame=TRUE)
health_w2 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew2_rel6-0-0_ALL_datasets_spss\\sharew2_rel6-0-0_gv_health.sav",
                      to.data.frame=TRUE)
health_w3 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew4_rel6-0-0_ALL_datasets_spss\\sharew4_rel6-0-0_gv_health.sav",
                      to.data.frame=TRUE)
health_w4 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew5_rel6-0-0_ALL_datasets_spss\\sharew5_rel6-0-0_gv_health.sav",
                      to.data.frame=TRUE)
health_w5 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew6_rel6-0-0_ALL_datasets_spss\\sharew6_rel6-0-0_gv_health.sav",
                      to.data.frame=TRUE)
# generalized variable isced
isced_w1 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew1_rel6-0-0_ALL_datasets_spss\\sharew1_rel6-0-0_gv_isced.sav",
                      to.data.frame=TRUE)
isced_w2 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew2_rel6-0-0_ALL_datasets_spss\\sharew2_rel6-0-0_gv_isced.sav",
                      to.data.frame=TRUE)
isced_w3 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew4_rel6-0-0_ALL_datasets_spss\\sharew4_rel6-0-0_gv_isced.sav",
                      to.data.frame=TRUE)
isced_w4 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew5_rel6-0-0_ALL_datasets_spss\\sharew5_rel6-0-0_gv_isced.sav",
                      to.data.frame=TRUE)
isced_w5 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew6_rel6-0-0_ALL_datasets_spss\\sharew6_rel6-0-0_gv_isced.sav",
                      to.data.frame=TRUE)
# technical variables
tech_w1 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew1_rel6-0-0_ALL_datasets_spss\\sharew1_rel6-0-0_technical_variables.sav",
                    to.data.frame=TRUE)
tech_w2 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew2_rel6-0-0_ALL_datasets_spss\\sharew2_rel6-0-0_technical_variables.sav",
                    to.data.frame=TRUE)
tech_w3 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew4_rel6-0-0_ALL_datasets_spss\\sharew4_rel6-0-0_technical_variables.sav",
                    to.data.frame=TRUE)
tech_w4 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew5_rel6-0-0_ALL_datasets_spss\\sharew5_rel6-0-0_technical_variables.sav",
                    to.data.frame=TRUE)
tech_w5 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew6_rel6-0-0_ALL_datasets_spss\\sharew6_rel6-0-0_technical_variables.sav",
                    to.data.frame=TRUE)

# edit cv so that only observations are selected that have information on physical welfare
cv_w1 <- cv_1[cv_1[,"mergeid"]%in%ph_w1[,"mergeid"],]
cv_w2 <- cv_2[cv_2[,"mergeid"]%in%ph_w2[,"mergeid"],]
cv_w3 <- cv_3[cv_3[,"mergeid"]%in%ph_w3[,"mergeid"],]
cv_w4 <- cv_4[cv_4[,"mergeid"]%in%ph_w4[,"mergeid"],]
cv_w5 <- cv_5[cv_5[,"mergeid"]%in%ph_w5[,"mergeid"],]

# merge the four waves on basis of ph, with attrition from first wave onwards
cv_mergeid <- list(mergeid1 = cv_w1[,"mergeid"], mergeid2 = cv_w2[,"mergeid"],
                   mergeid3 = cv_w3[,"mergeid"], mergeid4 = cv_w4[,"mergeid"], mergeid5 = cv_w5[,"mergeid"])
dn_mergeid <- list(mergeid1 = dn_w1[,"mergeid"], mergeid2 = dn_w2[,"mergeid"],
                   mergeid3 = dn_w3[,"mergeid"], mergeid4 = dn_w4[,"mergeid"], mergeid5 = dn_w5[,"mergeid"])
ch_mergeid <- list(mergeid1 = ch_w1[,"mergeid"], mergeid2 = ch_w2[,"mergeid"],
                   mergeid3 = ch_w3[,"mergeid"], mergeid4 = ch_w4[,"mergeid"], mergeid5 = ch_w5[,"mergeid"])
ph_mergeid <- list(mergeid1 = ph_w1[,"mergeid"], mergeid2 = ph_w2[,"mergeid"],
                   mergeid3 = ph_w3[,"mergeid"], mergeid4 = ph_w4[,"mergeid"], mergeid5 = ph_w5[,"mergeid"])
ep_mergeid <- list(mergeid1 = ep_w1[,"mergeid"], mergeid2 = ep_w2[,"mergeid"],
                   mergeid3 = ep_w3[,"mergeid"], mergeid4 = ep_w4[,"mergeid"], mergeid5 = ep_w5[,"mergeid"])
ac_mergeid <- list(mergeid1 = ac_w1[,"mergeid"], mergeid2 = ac_w2[,"mergeid"],
                   mergeid3 = ac_w3[,"mergeid"], mergeid4 = ac_w4[,"mergeid"], mergeid5 = ac_w5[,"mergeid"])
ho_mergeid <- list(mergeid1 = ho_w1[,"mergeid"], mergeid2 = ho_w2[,"mergeid"],
                   mergeid3 = ho_w3[,"mergeid"], mergeid4 = ho_w4[,"mergeid"], mergeid5 = ho_w5[,"mergeid"])
health_mergeid <- list(mergeid1 = health_w1[,"mergeid"], mergeid2 = health_w2[,"mergeid"],
                       mergeid3 = health_w3[,"mergeid"], mergeid4 = health_w4[,"mergeid"], mergeid5 = health_w5[,"mergeid"])
isced_mergeid <- list(mergeid1 = isced_w1[,"mergeid"], mergeid2 = isced_w2[,"mergeid"],
                       mergeid3 = isced_w3[,"mergeid"], mergeid4 = isced_w4[,"mergeid"], mergeid5 = isced_w5[,"mergeid"])
tech_mergeid <- list(mergeid1 = tech_w1[,"mergeid"], mergeid2 = tech_w2[,"mergeid"],
                     mergeid3 = tech_w3[,"mergeid"], mergeid4 = tech_w4[,"mergeid"], mergeid5 = tech_w5[,"mergeid"])

# edit imputations so that only one dataset is selected
impind <- seq(1,nrow(imp_w1),5) # only select one of the imputations provided
imp_w1 <- imp_w1[impind,]
impind <- seq(1,nrow(imp_w2),5) # only select one of the imputations provided
imp_w2 <- imp_w2[impind,]
impind <- seq(1,nrow(imp_w3),5) # only select one of the imputations provided
imp_w3 <- imp_w3[impind,]
impind <- seq(1,nrow(imp_w4),5) # only select one of the imputations provided
imp_w4 <- imp_w4[impind,]
impind <- seq(1,nrow(imp_w5),5) # only select one of the imputations provided
imp_w5 <- imp_w5[impind,]
imp_mergeid <- list(mergeid1 = imp_w1[,"mergeid"], mergeid2 = imp_w2[,"mergeid"],
                    mergeid3 = imp_w3[,"mergeid"], mergeid4 = imp_w4[,"mergeid"], mergeid5 = imp_w5[,"mergeid"])

# find indices of observations in different variables in different waves
vec <- unlist(ph_mergeid)
#vec <- unlist(cv_mergeid)
unique <- sort(unique(vec))
id <- sort(rep(unique,5))

cv_retention <- matrix(0,length(id),1)
dn_retention <- matrix(0,length(id),1)
ch_retention <- matrix(0,length(id),1)
ph_retention <- matrix(0,length(id),1)
ep_retention <- matrix(0,length(id),1)
ac_retention <- matrix(0,length(id),1)
ho_retention <- matrix(0,length(id),1)
health_retention <- matrix(0,length(id),1)
imp_retention <- matrix(0,length(id),1)
isced_retention <- matrix(0,length(id),1)
tech_retention <- matrix(0,length(id),1)

# assign indices to find the right wave
index1 <- seq(1,(-4 + length(id)),by=5)
index2 <- seq(2,(-3 + length(id)),by=5)
index3 <- seq(3,(-2 + length(id)),by=5)
index4 <- seq(4,(-1 + length(id)),by=5)
index5 <- seq(5,(     length(id)),by=5)
for(i in 1:5){ 
  index <- eval(parse(text=paste("index",i,sep="")))
  cv_retention[index] <- match(unique,cv_mergeid[[i]])
  dn_retention[index] <- match(unique,dn_mergeid[[i]])
  ch_retention[index] <- match(unique,ch_mergeid[[i]])
  ph_retention[index] <- match(unique,ph_mergeid[[i]]) 
  ep_retention[index] <- match(unique,ep_mergeid[[i]])
  ac_retention[index] <- match(unique,ac_mergeid[[i]])
  ho_retention[index] <- match(unique,ho_mergeid[[i]])
  health_retention[index] <- match(unique,health_mergeid[[i]])
  #imp_retention[index] <- match(unique,imp_mergeid[[i]])
  isced_retention[index] <- match(unique,isced_mergeid[[i]])
  tech_retention[index] <- match(unique,tech_mergeid[[i]])
}

# cbind.data.frame(ph_unique[1:40],(imp_retention[index2])[1:40],imp_mergeid[[2]][1:40])

# Pieter SHARE
#install.packages("readstata13")
#library(readstata13)
#Pieter_w1 <- read.dta13("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\SHARE\\2004.dta")
#Pieter_w2 <- read.dta13("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\SHARE\\2006.dta")
#Pieter_w3 <- read.dta13("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\SHARE\\2008.dta")
#Pieter_w4 <- read.dta13("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\SHARE\\2010.dta")
#Pieter_w5 <- read.dta13("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\SHARE\\2013.dta")

