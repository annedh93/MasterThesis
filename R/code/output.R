output1 <- function(ph_unique, ph_retention, cv_panel, health_panel){

self_health1 <- as.numeric(health_panel[complete.cases(health_panel[,"self-perceived health 1"]),"self-perceived health 1"])-2
self_health2 <- as.numeric(health_panel[complete.cases(health_panel[,"self-perceived health 2"]),"self-perceived health 2"])-2
self_health4 <- as.numeric(health_panel[complete.cases(health_panel[,"self-perceived health 4"]),"self-perceived health 4"])-2
self_health5 <- as.numeric(health_panel[complete.cases(health_panel[,"self-perceived health 5"]),"self-perceived health 5"])-2
self_health6 <- as.numeric(health_panel[complete.cases(health_panel[,"self-perceived health 6"]),"self-perceived health 6"])-2

total <- c((length(ph_unique)-sum(is.na(ph_retention[,1]))),(length(ph_unique)-sum(is.na(ph_retention[,2]))),
           (length(ph_unique)-sum(is.na(ph_retention[,3]))),(length(ph_unique)-sum(is.na(ph_retention[,4]))),
           (length(ph_unique)-sum(is.na(ph_retention[,5]))))
output <- t(cbind(c((table(cv_panel[,"gender2004"])[4]/total[1])*100,NA,NA,NA),
                   c((table(cv_panel[,"gender2007"])[4]/total[2])*100,NA,NA,NA),
                   c((table(cv_panel[,"gender2011"])[4]/total[3])*100,NA,NA,NA),
                   c((table(cv_panel[,"gender2013"])[4]/total[4])*100,NA,NA,NA),
                   c((table(cv_panel[,"gender2015"])[4]/total[5])*100,NA,NA,NA),
                   c(mean(cv_panel[complete.cases(cv_panel[,"age2004"]),"age2004"]),
                     sd(cv_panel[complete.cases(cv_panel[,"age2004"]),"age2004"]),
                     23,
                     max(cv_panel[complete.cases(cv_panel[,"age2004"]),"age2004"])),
                   c(mean(cv_panel[complete.cases(cv_panel[,"age2007"]),"age2007"]),
                     sd(cv_panel[complete.cases(cv_panel[,"age2007"]),"age2007"]),
                     15,
                     max(cv_panel[complete.cases(cv_panel[,"age2007"]),"age2007"])),
                   c(mean(cv_panel[complete.cases(cv_panel[,"age2011"]),"age2011"]),
                     sd(cv_panel[complete.cases(cv_panel[,"age2011"]),"age2011"]),
                     12,
                     max(cv_panel[complete.cases(cv_panel[,"age2011"]),"age2011"])),
                   c(mean(cv_panel[complete.cases(cv_panel[,"age2013"]),"age2013"]),
                     sd(cv_panel[complete.cases(cv_panel[,"age2013"]),"age2013"]),
                     22,
                     max(cv_panel[complete.cases(cv_panel[,"age2013"]),"age2013"])),
                   c(mean(cv_panel[complete.cases(cv_panel[,"age2015"]),"age2015"]),
                     sd(cv_panel[complete.cases(cv_panel[,"age2015"]),"age2015"]),
                     22,
                     max(cv_panel[complete.cases(cv_panel[,"age2015"]),"age2015"])),
                   c(mean(health_panel[complete.cases(health_panel[,"casp 1"]),"casp 1"]),
                     sd(health_panel[complete.cases(health_panel[,"casp 1"]),"casp 1"]),
                     min(health_panel[complete.cases(health_panel[,"casp 1"]),"casp 1"]),
                     max(health_panel[complete.cases(health_panel[,"casp 1"]),"casp 1"])),
                   c(mean(health_panel[complete.cases(health_panel[,"casp 2"]),"casp 2"]),
                     sd(health_panel[complete.cases(health_panel[,"casp 2"]),"casp 2"]),
                     min(health_panel[complete.cases(health_panel[,"casp 2"]),"casp 2"]),
                     max(health_panel[complete.cases(health_panel[,"casp 2"]),"casp 2"])),
                   c(mean(health_panel[complete.cases(health_panel[,"casp 4"]),"casp 4"]),
                     sd(health_panel[complete.cases(health_panel[,"casp 4"]),"casp 4"]),
                     min(health_panel[complete.cases(health_panel[,"casp 4"]),"casp 4"]),
                     max(health_panel[complete.cases(health_panel[,"casp 4"]),"casp 4"])),
                   c(mean(health_panel[complete.cases(health_panel[,"casp 5"]),"casp 5"]),
                     sd(health_panel[complete.cases(health_panel[,"casp 5"]),"casp 5"]),
                     min(health_panel[complete.cases(health_panel[,"casp 5"]),"casp 5"]),
                     max(health_panel[complete.cases(health_panel[,"casp 5"]),"casp 5"])),
                   c(mean(health_panel[complete.cases(health_panel[,"casp 6"]),"casp 6"]),
                     sd(health_panel[complete.cases(health_panel[,"casp 6"]),"casp 6"]),
                     min(health_panel[complete.cases(health_panel[,"casp 6"]),"casp 6"]),
                     max(health_panel[complete.cases(health_panel[,"casp 6"]),"casp 6"])),
                   c(mean(self_health1),sd(self_health1),NA,NA),
                   c(mean(self_health2),sd(self_health2),NA,NA),
                   c(mean(self_health4),sd(self_health4),NA,NA),
                   c(mean(self_health5),sd(self_health5),NA,NA),
                   c(mean(self_health6),sd(self_health6),NA,NA)))
colnames(output) <- c("Mean","Std. deviation","Min.","Max")
rownames(output) <- c("Gender 2004","Gender 2007","Gender 2011","Gender 2013","Gender 2015",
                       "Age 2004","Age 2007", "Age 2011", "Age 2013", "Age 2015",
                       "casp 2004","casp 2007", "casp 2011","casp 2013","casp 2015",
                       "self-perceived health 1","self-perceived health 2","self-perceived health 4","self-perceived health 5","self-perceived health 6")
return(output)
}

output2 <- function(ph_unique,ph_retention,diseases_panel){
    total <- c((length(ph_unique)-sum(is.na(ph_retention[,1]))),(length(ph_unique)-sum(is.na(ph_retention[,2]))),
               (length(ph_unique)-sum(is.na(ph_retention[,3]))),(length(ph_unique)-sum(is.na(ph_retention[,4]))),
               (length(ph_unique)-sum(is.na(ph_retention[,5]))))
    chronic <- c(table(diseases_panel[,2])[3], table(diseases_panel[,3])[3], 
                 table(diseases_panel[,4])[3], table(diseases_panel[,5])[3], table(diseases_panel[,6])[3])
    notchronic <- c(table(diseases_panel[,2])[4], table(diseases_panel[,3])[4], 
                    table(diseases_panel[,4])[4], table(diseases_panel[,5])[4], table(diseases_panel[,6])[4])
    output2 <- t(cbind(total,
                       100*(chronic/total),
                       100*(notchronic/total),
                       c((table(diseases_panel[,7])[4]/chronic[1])*100, (table(diseases_panel[,8])[4]/chronic[2])*100,
                         (table(diseases_panel[,9])[4]/chronic[3])*100, (table(diseases_panel[,10])[4]/chronic[4])*100,(table(diseases_panel[,11])[4]/chronic[5])*100),
                       c((table(diseases_panel[,12])[4]/chronic[1])*100, (table(diseases_panel[,13])[4]/chronic[2])*100,
                         (table(diseases_panel[,14])[4]/chronic[3])*100, (table(diseases_panel[,15])[4]/chronic[4])*100,(table(diseases_panel[,16])[4]/chronic[5])*100),
                       c((table(diseases_panel[,17])[4]/chronic[1])*100, (table(diseases_panel[,18])[4]/chronic[2])*100,
                         (table(diseases_panel[,19])[4]/chronic[3])*100, (table(diseases_panel[,20])[4]/chronic[4])*100,(table(diseases_panel[,21])[4]/chronic[5])*100),                   
                       c((table(diseases_panel[,22])[4]/chronic[1])*100, (table(diseases_panel[,23])[4]/chronic[2])*100,
                         (table(diseases_panel[,24])[4]/chronic[3])*100, (table(diseases_panel[,25])[4]/chronic[4])*100,(table(diseases_panel[,26])[4]/chronic[5])*100),
                       c((table(diseases_panel[,27])[4]/chronic[1])*100, (table(diseases_panel[,28])[4]/chronic[2])*100,
                         (table(diseases_panel[,29])[4]/chronic[3])*100, (table(diseases_panel[,30])[4]/chronic[4])*100,(table(diseases_panel[,31])[4]/chronic[5])*100),
                       c((table(diseases_panel[,32])[4]/chronic[1])*100, (table(diseases_panel[,33])[4]/chronic[2])*100,
                         (table(diseases_panel[,34])[4]/chronic[3])*100, (table(diseases_panel[,35])[4]/chronic[4])*100,(table(diseases_panel[,36])[4]/chronic[5])*100),
                       c((table(diseases_panel[,37])[4]/chronic[1])*100, (table(diseases_panel[,38])[4]/chronic[2])*100,
                         (table(diseases_panel[,39])[4]/chronic[3])*100, NA,NA),
                       c((table(diseases_panel[,40])[4]/chronic[1])*100, (table(diseases_panel[,41])[4]/chronic[2])*100,
                         (table(diseases_panel[,42])[4]/chronic[3])*100, NA,NA),
                       c((table(diseases_panel[,43])[4]/chronic[1])*100, (table(diseases_panel[,44])[4]/chronic[2])*100,
                         (table(diseases_panel[,45])[4]/chronic[3])*100, NA,NA),
                       c((table(diseases_panel[,46])[4]/chronic[1])*100, (table(diseases_panel[,47])[4]/chronic[2])*100,
                         (table(diseases_panel[,48])[4]/chronic[3])*100, (table(diseases_panel[,49])[4]/chronic[4])*100,(table(diseases_panel[,50])[4]/chronic[5])*100),
                       c((table(diseases_panel[,51])[4]/chronic[1])*100, (table(diseases_panel[,52])[4]/chronic[2])*100,
                         (table(diseases_panel[,53])[4]/chronic[3])*100, (table(diseases_panel[,54])[4]/chronic[4])*100,(table(diseases_panel[,55])[4]/chronic[5])*100),
                       c((table(diseases_panel[,56])[4]/chronic[1])*100, (table(diseases_panel[,57])[4]/chronic[2])*100,
                         (table(diseases_panel[,58])[4]/chronic[3])*100, (table(diseases_panel[,59])[4]/chronic[4])*100,(table(diseases_panel[,60])[4]/chronic[5])*100),
                       c((table(diseases_panel[,61])[4]/chronic[1])*100, (table(diseases_panel[,62])[4]/chronic[2])*100,
                         (table(diseases_panel[,63])[4]/chronic[3])*100, (table(diseases_panel[,64])[4]/chronic[4])*100,(table(diseases_panel[,65])[4]/chronic[5])*100),
                       c((table(diseases_panel[,66])[4]/chronic[1])*100, (table(diseases_panel[,67])[4]/chronic[2])*100,
                         (table(diseases_panel[,68])[4]/chronic[3])*100, (table(diseases_panel[,69])[4]/chronic[4])*100,(table(diseases_panel[,70])[4]/chronic[5])*100)))
    rownames(output2) <- c("Total","Chronically ill","Not chronically ill",
                           "A heart attack", "High blood pressure or hypertension", "High cholesterol","A stroke or cerebral vascular disease",
                           "Diabetes or high blood sugar","Chronic lung disease","Asthma","Arthritis",
                           "Osteoperosis", "Cancer or malignant tumor","Stomach or duodenal ulcer",
                           "Parkinson disease","Cataracts","Hip fracture")
    colnames(output2) <- c("2004","2007","2011","2013","2015")
    return(output2)
}