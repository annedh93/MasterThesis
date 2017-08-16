retention <- function(id,
                         cv_w1,cv_w2,cv_w3,cv_w4,cv_w5,
                         dn_w1,dn_w2,dn_w3,dn_w4,dn_w5,
                         ch_w1,ch_w2,ch_w3,ch_w4,ch_w5,
                         ph_w1,ph_w2,ph_w3,ph_w4,ph_w5,
                         ep_w1,ep_w2,ep_w3,ep_w4,ep_w5,
                         ac_w1,ac_w2,ac_w3,ac_w4,ac_w5,
                         ho_w1,ho_w2,ho_w3,ho_w4,ho_w5,
                         health_w1,health_w2,health_w3,health_w4,health_w5,
                         isced_w1,isced_w2,isced_w3,isced_w4,isced_w5,
                         tech_w1,tech_w2,tech_w3,tech_w4,tech_w5,
                         cv_retention,dn_retention,ch_retention,ph_retention,ep_retention,ac_retention,ho_retention,
                         health_retention,isced_retention,tech_retention,
                         index1,index2,index3,index4,index5){
  # define variables that will be put in the data frame
  #names <- cbind("id","wave","mergeid","coupleid","country","firstwave",
  #               "gender","age","partnerinhh","mstatc","mstat","mstat2","yedu","isced","child","child2","hhsize",
  #               "casp","sphus","chronic","chronicn","limac","ADL","IADL",
  #               "maxgrip","maxgrip2","employment","lifesat","tenure","nursing")
  names <- cbind("id","numberid","wave","mergeid","coupleid","numbercid","country","firstwave",
                 "gender","age","partnerinhh","mstatc","mstat","isced","child","hhsize",
                 "casp","sphus","chronic","chronicn","limac","ADL","IADL",
                 "employment","lifesat","tenure","nursing")
  output <- as.data.frame(matrix(NA,length(id),length(names)))
  colnames(output) <- names
  age <- c("age2004","age2007","age2011","age2013","age2015")
  coupleid <- c("coupleid1","coupleid2","coupleid4","coupleid5","coupleid6")
  isced <- c("isced1997_r","isced1997_r","isced1997_r","isced1997_r")
  chronicn <- c("chronicw1","chronicw2","chronicw4","chronicw5","chronicw6c")
  
  # fill the data frame
  output[,"id"] <- as.character(id)
  output[,"wave"] <- rep(c(1,2,4,5,6),length(id)/5)
  for(i in 1:5){ 
    index <- eval(parse(text=paste("index",i,sep="")))
    cv <- eval(parse(text=paste("cv_w",i,sep="")))
    dn <- eval(parse(text=paste("dn_w",i,sep="")))
    ch <- eval(parse(text=paste("ch_w",i,sep="")))
    ph <- eval(parse(text=paste("ph_w",i,sep="")))
    ep <- eval(parse(text=paste("ep_w",i,sep="")))
    ac <- eval(parse(text=paste("ac_w",i,sep="")))
    ho <- eval(parse(text=paste("ho_w",i,sep="")))
    health <- eval(parse(text=paste("health_w",i,sep="")))
    #imp <- eval(parse(text=paste("imp_w",i,sep="")))
    isced <- eval(parse(text=paste("isced_w",i,sep="")))
    tech <- eval(parse(text=paste("tech_w",i,sep="")))
    output[index,"mergeid"] <- as.character(cv[cv_retention[index],"mergeid"])
    output[index,"coupleid"] <- as.character(cv[cv_retention[index],coupleid[i]])
    output[index,"firstwave"] <- cv[cv_retention[index],"firstwave"]
    output[index,"country"] <- as.character(cv[cv_retention[index],"country"])
    output[index,"gender"] <- as.character(cv[cv_retention[index],"gender"])
    output[index,"age"] <- cv[cv_retention[index],age[i]]
    output[index,"partnerinhh"] <- cv[cv_retention[index],"partnerinhh"]
    if(i > 1){
      output[index,"mstatc"] <- dn[dn_retention[index],"dn044_"]  
    }
    output[index,"mstat"] <- dn[dn_retention[index],"dn014_"]
    if(i==1){
      output[index,"isced"] <- isced[isced_retention[index],"isced1997y_r"]
    }
    else{
      output[index,"isced"] <- dn[dn_retention[index],"dn041_"]
    }
    output[index,"child"] <- ch[ch_retention[index],"ch001_"]
    output[index,"hhsize"] <- cv[cv_retention[index],"hhsize"]
    output[index,"casp"] <- health[health_retention[index],"casp"]
    output[index,"sphus"] <- health[health_retention[index],"sphus"]
    output[index,"chronic"] <- ph[ph_retention[index],"ph004_"]
    output[index,"chronicn"] <- health[health_retention[index],chronicn[i]]
    output[index,"limac"] <- ph[ph_retention[index],"ph005_"]
    output[index,"ADL"] <- health[health_retention[index],"adl"]
    output[index,"IADL"] <- health[health_retention[index],"iadl"]
    output[index,"employment"] <- ep[ep_retention[index],"ep005_"]
    if(i>1){
      output[index,"lifesat"] <- ac[ac_retention[index],"ac012_"]
      output[index,"nursing"] <- tech[tech_retention[index],"mn024_"]
    }
    output[index,"tenure"] <- ho[ho_retention[index],"ho002_"]
  }
  
  output[which(output[,"age"]<0 & output[,"age"]!=-9),"age"] <- NA
  output[which(output[,"partnerinhh"]<4),"partnerinhh"] <- NA
  output[which(output[,"mstatc"]<3),"mstatc"] <- NA
  output[which(output[,"mstat"]<3),"mstat"] <- NA
  output[which(output[,"isced"]<0 | output[,"isced"]>94),"isced"] <- NA
  output[which(output[,"child"]<0),"child"] <- NA
  output[which(output[,"sphus"]<3),"sphus"] <- NA
  output[which(output[,"chronic"]<3),"chronic"] <- NA
  output[which(output[,"chronicn"]<0),"chronicn"] <- NA
  output[which(output[,"limac"]<3),"limac"] <- NA
  output[which(output[,"ADL"]<0),"ADL"] <- NA
  output[which(output[,"IADL"]<0),"IADL"] <- NA
  output[which(output[,"employment"]<3),"employment"] <- NA
  output[which(output[,"lifesat"]<0),"lifesat"] <- NA
  output[which(output[,"tenure"]<3),"tenure"] <- NA
  output[which(output[,"nursing"]<3),"nursing"] <- NA
  
  output[,"numberid"] <- as.numeric(as.factor(output[,"id"]))
  output[,"numbercid"] <- as.numeric(as.factor(output[,"coupleid"]))
  return(output)
}

retention_diseases <- function(ph_waves,ph_w1,ph_w2,ph_w4,ph_w5,ph_w6,ph_retention,index1,index2,index3,index4,index5){
    # define variables that will be put in the data frame
    names <- c("id","wave","chronically ill",
               "heart attack",
               "high blood pressure or hypertension",
               "high blood cholesterol",
               "stroke or cerebral vascular disease",
               "diabetes or high blood sugar",
               "chronic lung disease",
               "Asthma",
               "Arthritis",
               "Osteoporosis",
               "cancer or malignant tumor",
               "stomach or duodenal ulcer",
               "Parkinson disease",
               "cataracts",
               "hip fracture",
               "other fractures",
               "Alzheimer's disease",
               "benign tumor",
               "Other affective or emotional disorders",
               "Rheumatoid Arthritis",
               "Osteoarthritis or other rheumatism",
               "Chronic kidney disease")
    output <- as.data.frame(matrix(NA,length(ph_waves),length(names)))
    colnames(output) <- names

    # fill the data frame
    output[,"id"] <- as.character(ph_waves)
    output[,"wave"] <- rep(c(1,2,4,5,6),length(ph_waves)/5)
    for(i in 1:5){ 
      index <- eval(parse(text=paste("index",i,sep="")))
      ph <- eval(parse(text=paste("ph_w",i,sep="")))
      output[index,"chronically ill"] <- ph[ph_retention[index],"ph004_"]
      output[index,"heart attack"] <- ph[ph_retention[index],"ph006d1"]
      output[index,"high blood pressure or hypertension"] <- ph[ph_retention[index],"ph006d2"]
      output[index,"high blood cholesterol"] <- ph[ph_retention[index],"ph006d3"]
      output[index,"stroke or cerebral vascular disease"] <- ph[ph_retention[index],"ph006d4"]
      output[index,"diabetes or high blood sugar"] <- ph[ph_retention[index],"ph006d5"]
      output[index,"chronic lung disease"] <- ph[ph_retention[index],"ph006d6"]
      if(i == 3){
        output[index,"Arthritis"] <- as.numeric(ph[ph_retention[index],"ph006d8"])
      }
      if(i <= 2){
        output[index,"Asthma"] <- as.numeric(ph[ph_retention[index],"ph006d7"])
        output[index,"Arthritis"] <- as.numeric(ph[ph_retention[index],"ph006d8"])
        output[index,"Osteoporosis"] <- as.numeric(ph[ph_retention[index],"ph006d9"])
      }
      output[index,"cancer or malignant tumor"] <- as.numeric(ph[ph_retention[index],"ph006d10"])
      output[index,"stomach or duodenal ulcer"] <- as.numeric(ph[ph_retention[index],"ph006d11"])
      output[index,"Parkinson disease"] <- as.numeric(ph[ph_retention[index],"ph006d12"])
      output[index,"cataracts"] <- as.numeric(ph[ph_retention[index],"ph006d13"])
      output[index,"hip fracture"] <- as.numeric(ph[ph_retention[index],"ph006d14"])
      if(i >= 2){
        output[index,"other fractures"] <- as.numeric(ph[ph_retention[index],"ph006d15"])
        output[index,"Alzheimer's disease"] <- as.numeric(ph[ph_retention[index],"ph006d16"])
      }
      if(i == 2){
        output[index,"benign tumor"] <- as.numeric(ph[ph_retention[index],"ph006d17"])
      }
      if(i >= 4){
        output[index,"Other affective or emotional disorders"] <- as.numeric(ph[ph_retention[index],"ph006d18"])
        output[index,"Rheumatoid Arthritis"] <- as.numeric(ph[ph_retention[index],"ph006d19"])
        output[index,"Osteoarthritis or other rheumatism"] <- as.numeric(ph[ph_retention[index],"ph006d20"])
      }
      if(i == 5){
        output[index,"Chronic kidney disease"] <- as.numeric(ph[ph_retention[index],"ph006d21"])
      }
    }
    
  return(output)
}