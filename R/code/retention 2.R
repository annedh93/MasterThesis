retention <- function(ph_waves,
                         cv_w1,cv_w2,cv_w3,cv_w4,cv_w5,
                         dn_w1,dn_w2,dn_w3,dn_w4,dn_w5,
                         ch_w1,ch_w2,ch_w3,ch_w4,ch_w5,
                         ph_w1,ph_w2,ph_w3,ph_w4,ph_w5,
                         ep_w1,ep_w2,ep_w3,ep_w4,ep_w5,
                         health_w1,health_w2,health_w3,health_w4,health_w5,
                         imp_w1,imp_w2,imp_w3,imp_w4,imp_w5,
                         cv_retention,dn_retention,ch_retention,ph_retention,ep_retention,health_retention,imp_retention,
                         index1,index2,index3,index4,index5){
  # define variables that will be put in the data frame
  names <- cbind("id","wave","mergeid","country","firstwave",
                 "gender","age","partner in household","marital status","number of children","household size",
                 "casp","self perceived health","chronic","chronic2","limited in activities","ADL","IADL",
                 "employment")
  output <- as.data.frame(matrix(0,length(ph_waves),length(names)))
  colnames(output) <- names
  age <- c("age2004","age2007","age2011","age2013","age2015")
  
  # fill the data frame
  output[,"id"] <- as.character(ph_waves)
  output[,"wave"] <- rep(c(1,2,4,5,6),length(ph_waves)/5)
  for(i in 1:5){ 
    index <- eval(parse(text=paste("index",i,sep="")))
    cv <- eval(parse(text=paste("cv_w",i,sep="")))
    dn <- eval(parse(text=paste("dn_w",i,sep="")))
    ch <- eval(parse(text=paste("ch_w",i,sep="")))
    ph <- eval(parse(text=paste("ph_w",i,sep="")))
    ep <- eval(parse(text=paste("ep_w",i,sep="")))
    health <- eval(parse(text=paste("health_w",i,sep="")))
    imp <- eval(parse(text=paste("imp_w",i,sep="")))
    output[index,"mergeid"] <- as.character(cv[cv_retention[index],"mergeid"])
    output[index,"firstwave"] <- cv[cv_retention[index],"firstwave"]
    output[index,"country"] <- as.character(cv[cv_retention[index],"country"])
    output[index,"gender"] <- as.character(cv[cv_retention[index],"gender"])
    output[index,"age"] <- cv[cv_retention[index],age[i]]
    output[index,"partner in household"] <- cv[cv_retention[index],"partnerinhh"]
    output[index,"marital status"] <- as.character(dn[dn_retention[index],"dn014_"])
    output[index,"number of children"] <- ch[ch_retention[index],"ch001_"]
    output[index,"household size"] <- cv[cv_retention[index],"hhsize"]
    output[index,"casp"] <- health[health_retention[index],"casp"]
    output[index,"self perceived health"] <- health[health_retention[index],"sphus"]
    output[index,"chronic"] <- ph[ph_retention[index],"ph004_"]
    output[index,"limited in activities"] <- ph[ph_retention[index],"ph005_"]
    output[index,"ADL"] <- health[health_retention[index],"adl"]
    output[index,"IADL"] <- health[health_retention[index],"iadl"]
    output[index,"chronic2"] <- imp[imp_retention[index],"chronic"]
    output[index,"employment"] <- ep[ep_retention[index],"ep005_"]
  }

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