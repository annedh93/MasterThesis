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
  cv_w4 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew4_rel6-0-0_ALL_datasets_spss\\sharew4_rel6-0-0_cv_r.sav",
                    to.data.frame=TRUE)
  cv_w5 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew5_rel6-0-0_ALL_datasets_spss\\sharew5_rel6-0-0_cv_r.sav",
                        to.data.frame=TRUE)
  cv_w6 = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharew6_rel6-0-0_ALL_datasets_spss\\sharew6_rel6-0-0_cv_r.sav",
                    to.data.frame=TRUE)
  coverscreen = read.spss("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\sharewX_rel6-0-0_gv_allwaves_cv_r_spss\\sharewX_rel6-0-0_gv_allwaves_cv_r.sav",
                          to.data.frame=TRUE)
  
## try to merge the four waves on basis of ph, with attrition from first wave onwards
ph_mergeid <- list(mergeid1 = ph_w1[,"mergeid"], mergeid2 = ph_w2[,"mergeid"],
                mergeid4 = ph_w4[,"mergeid"], mergeid5 = ph_w5[,"mergeid"], mergeid6 = ph_w6[,"mergeid"])
cv_mergeid <- list(mergeid1 = cv_w1[,"mergeid"], mergeid2 = cv_w2[,"mergeid"],
                   mergeid4 = cv_w4[,"mergeid"], mergeid5 = cv_w5[,"mergeid"], mergeid6 = cv_w6[,"mergeid"])
health_mergeid <- list(mergeid1 = health_w1[,"mergeid"], mergeid2 = health_w2[,"mergeid"],
                       mergeid4 = health_w4[,"mergeid"], mergeid5 = health_w5[,"mergeid"], mergeid6 = health_w6[,"mergeid"])
phvec <- unlist(ph_mergeid)
ph_unique <- sort(unique(phvec))
ph_retention <- matrix(0,length(ph_unique),5)
cv_retention <- matrix(0,length(ph_unique),5)
health_retention <- matrix(0,length(ph_unique),5)
for(i in 1:5){ 
  ph_retention[,i] <- match(ph_unique,ph_mergeid[[i]]) 
  cv_retention[,i] <- match(ph_unique,cv_mergeid[[i]])
  health_retention[,i] <- match(ph_unique,health_mergeid[[i]])
}

# construct panels and tables
cv_panel <- retention_cv(ph_unique,cv_w1,cv_w2,cv_w4,cv_w5,cv_w6,cv_retention)
diseases_panel <- retention_diseases(ph_unique,ph_w1,ph_w2,ph_w4,ph_w5,ph_w6,ph_retention)
health_panel <- retention_health(ph_unique,health_w1,health_w2,health_w4,health_w5,health_w6,health_retention)

# descriptives general
output1(ph_unique, ph_retention, cv_panel, health_panel)
# descriptives prevalence chronic illneses
output2(ph_unique,ph_retention,diseases_panel)

# retention rates
retention12 <- 100*(sum(ph_mergeid[[1]]%in%ph_mergeid[[2]])/length(ph_mergeid[[1]]))
retention24 <- 100*(sum(ph_mergeid[[2]]%in%ph_mergeid[[3]])/length(ph_mergeid[[2]]))
retention45 <- 100*(sum(ph_mergeid[[3]]%in%ph_mergeid[[4]])/length(ph_mergeid[[3]]))
retention56 <- 100*(sum(ph_mergeid[[4]]%in%ph_mergeid[[5]])/length(ph_mergeid[[4]]))

######################################################################################################################################
c_cv_panel <- cv_panel[complete.cases(cv_panel[,10:13]),] # retain complete cases on the chronically ill variable

waves <- cbind(1:4)
c_cv_panel <- cv_panel[complete.cases(cv_panel),]
# make plot of general health retained individuals over waves
samp <- sort(sample(1:nrow(c_cv_panel),100))
cret <- c_cv_panel[samp,]
colfunc <- colorRampPalette(c("blue","red"))
colfunc <- colfunc(length(samp))
sorted <- cret[sort.list(cret[,6]),]
# plot for general health over four waves
plot(waves,sorted[1,2:5],type="l",ylim = range(c_cv_panel[,2:5]),col=colfunc[1], 
     xlab = "waves",ylab = "general health (3 = Excellent)")
for(i in 2:length(samp)){lines(waves,jitter(as.numeric(sorted[i,2:5])),col=colfunc[i])}
# plot for chronic illness over four waves
plot(waves,sorted[1,6:9],type="l",ylim = range(c_cv_panel[,6:9]),col=colfunc[1], 
     xlab = "waves",ylab = "chronic illness (3 = Yes)")
for(i in 2:length(samp)){lines(waves,jitter(as.numeric(sorted[i,6:9])),col=colfunc[i])}

# chronic variable gv_health
# plot chronic health over time
mergeid <- list(mergeid1 = health_w1[,"mergeid"], mergeid2 = health_w2[,"mergeid"],
                mergeid4 = health_w4[,"mergeid"], mergeid5 = health_w5[,"mergeid"])
retention <- matrix(0,length(mergeid[[1]]),4)
for(i in 1:4){ retention[,i] <- match(mergeid[[1]],mergeid[[i]]) }
retention <- retention[complete.cases(retention),]
health_panel <- cbind.data.frame(as.numeric(health_w1[retention[,1],"chronicw1"]>0),
                                 as.numeric(health_w2[retention[,2],"chronicw2"]>0),
                                 as.numeric(health_w4[retention[,3],"chronicw4"]>0),
                                 as.numeric(health_w5[retention[,4],"chronicw5"]>0))
colnames(health_panel) <- cbind("chr 1", "chr 2", "chr 4", "chr 5")

nc_health <- health_panel[which(health_panel[,1]==0),]
ill_health <- nc_health[which(rowSums(nc_health[,2:3])>0),]
cret <- ill_health[complete.cases(ill_health),]
samp <- sort(sample(1:nrow(cret),100))
cret <- cret[samp,]
colfunc <- colorRampPalette(c("blue","red"))
colfunc <- colfunc(length(samp))
sorted <- cret[sort.list(cret[,1]),]
plot(waves, sorted[1,],type="l", xaxt = "n", ylim = range(c(0,1)),col=colfunc[1], 
     xlab = "waves",ylab = "chronic illness (1 = Yes)")
for(i in 2:length(samp)){lines(waves,jitter(as.numeric(sorted[i,])),col=colfunc[i])}
axis(1, at = 1:4, labels = c(1,2,4,5))
