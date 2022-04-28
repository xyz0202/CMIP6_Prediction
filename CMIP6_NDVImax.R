library(lubridate)

#============ read in dataset =================
setwd("./Desktop/work/CMIP6data/Australia/CMIP6demo")
pr    <- read.csv("pr_demo.csv", header = F, check.names = F)
tas   <- read.csv("tas_demo.csv", header = F, check.names = F)
ndvi  <- read.csv("ndvi_demo.csv", header = F, check.names = F)
par   <- read.csv("AttributionResults_demo.csv", header = TRUE, check.names = F)
max_n <- read.csv("Max_NDVI.csv", header = TRUE, check.names = F)

#=============================
#ndvi_ts <- ts(ndvi,start = c(1982,1), end=c(2014,12))
#pr_ts   <- ts(pr,start = c(1960,1), end=c(2014,12))
#tas_ts  <- ts(tas,start = c(1960,1), end=c(2014,12))
#yst <- start(ndvi_ts)[1]
#mst <- start(ndvi_ts)[2]
#============ read in parameters =================
osp  <- par$osp
acp  <- par$acp
tosp <- par$tosp
tacp <- par$tacp
#============ calculate optimal acc p =================

  #========== get the start/max location ==========
    str_loc <- c()
    max_loc <- data.frame()
  #loc <- function(str_loc,max_loc,max_mon)
  
    for (i in (seq(2,dim(ndvi)[2],by=12))){
      for (j in (seq(2,dim(pr)[2],by=12))){
        if (year(pr[1,j])==year(ndvi[1,i])){
          str_loc <- c(str_loc,j)
        }
      } 
    }
    for (n in 1:length(str_loc)){
      for (k in 1:dim(max_n)[1]) {
        max_loc[k,n] <- str_loc[n]+max_n[k,2]-1
      }
    }
  #============ calculate the sum of p =====================
  #oap_calculator <- function(op_acp, op_avt, stan_p, stan_t)
  op_acp <- data.frame()
 
  for (pn in 2:dim(pr)[1]) {
    for (mr in 1:dim(max_loc)[2]){
      str_acp       <- max_loc[pn,mr]-osp[pn-1]-acp[pn-1]
      end_acp       <- max_loc[pn,mr]-osp[pn-1]
      if (!is.na(str_acp) & !is.na(end_acp)){
      ocp             <- as.numeric(pr[pn,str_acp:end_acp])
      #ocp_mean        <- mean(ocp)
      #ocp_standev     <- sd(ocp)
      op_acp[pn-1,mr] <- sum(ocp)
      #std_p[pn-1,mr]  <- scale(ocp, center = T, scale = T )
      } else{
      op_acp[pn-1,mr] <- NA
      }
    }
  }
  
  #============ calculate the mean of t ==================
  op_avt <- data.frame()
   
  for (tn in 2:dim(tas)[1]) {
    for (mr in 1:dim(max_loc)[2]) {
        str_avt      <- max_loc[tn,mr]-tosp[tn-1]-tacp[tn-1]
        end_avt      <- max_loc[tn,mr]-tosp[tn-1]
      }
        if(!is.na(str_avt) & !is.na(end_avt)){
        avt             <- as.numeric(tas[tn,str_avt:end_avt])
        op_avt[tn-1,mr] <- mean(avt)
        #std_t[pn,mr]  <- scale(as.numeric(tas[tn,str_avt:end_avt]), center = T, scale = T )
      }else{
        op_avt[tn-1,mr] <- NA
      }
    }
 
#=========== if significant break points, standard score will be calculated ===
    #============ standard score =================
  BH  <- as.numeric(par$Break.Height)
  SC  <- as.numeric(par$Slope.Change)
  SCT <- as.numeric(par$Slope.ChangeTmp)
  #============= p standard score================
  std_p   <- data.frame()
  std_t   <- data.frame()
  for (p in 1:dim(op_acp)[1]) {
    for (s in 1:dim(op_acp)[2]) {
      if (!is.na(BH[p])){
        if (!is.na(op_acp[p,s])){
          mean_p     <- mean(as.numeric(op_acp[p,]))
          dev_p      <- sd(as.numeric(op_acp[p,]))
          std_p[p,s] <- (as.numeric(op_acp[p,s])-mean_p)/dev_p
          #std_p[p,s] <- scale(op_acp[p,s],center = op_acp[p,],scale = op_acp[p,]) 
        }else{
          std_p[p,s]  <- NA
        }
      }else{
        std_p[p,s] <- op_acp[p,s]
      }
    }
  }
  #============= t standard score================
  for (t in 1:dim(op_avt)[1]) {
    for (st in 1:dim(op_avt)[2]) {
      if (!is.na(BH[t])){
        if (!is.na(op_avt[t,st])){
          mean_t      <- mean(as.numeric(op_avt[t,]))
          dev_t       <- sd(as.numeric(op_avt[t,]))
          std_t[t,st] <- (as.numeric(op_avt[t,st])-mean_t)/dev_t
          #std_p[p,s] <- scale(op_acp[p,s],center = op_acp[p,],scale = op_acp[p,]) 
        }else{
          std_t[t,st]  <- NA
        }
      }else{
        std_t[t,st] <- op_avt[t,st]
      }
    }
  }

  #=========== use parameter calculate NDVImax ===
    sl_p  <- as.numeric(par$pre.slope)
    sl_t  <- as.numeric(par$temp.slope)  
    b     <- as.numeric(par$intercept)  
    #parm  <- c(sl_p,sl_t,BH,b,SC,SCT)
    ndvi_cm6 <- data.frame()
    
    for (x in 1:dim(std_p)[1]) {
      for (y in 1:dim(std_p)[2]) {
        if (is.na(std_p[x,y])){
          ndvi_cm6[x,y] <- NA
        }else  if (is.na(BH[x])){
          if(is.na(sl_t[x])|is.na(std_t[x,y])){
            ndvi_cm6[x,y] <- b[x] + sl_p[x]*as.numeric(std_p[x,y]) 
          }else{
            ndvi_cm6[x,y] <- b[x] + sl_p[x]*std_p[x,y] + sl_t[x]*as.numeric(std_t[x,y]) 
          }
        }else{
          if(is.na(sl_t[x])|is.na(sl_t[x])){
            ndvi_cm6[x,y] <- b[x] + sl_p[x]*as.numeric(std_p[x,y]) +  BH[x] + SC[x]*as.numeric(std_p[x,y]) 
          }else{
            ndvi_cm6[x,y] <- b[x] + sl_p[x]*as.numeric(std_p[x,y]) + sl_t[x]*as.numeric(std_t[x,y])+ BH[x] + SC[x]*as.numeric(std_p[x,y]) + SCT[x]*sl_t[x]*as.numeric(std_t[x,y]) 
          }
        }
      }
      
    }
          
    #=========== write out NDVImax ======== 
    #date   <- seq.Date(from = as.Date("1982",format="%Y"), by="year", length.out = 33)
    row.names(ndvi_cm6) <- ndvi[-1,1]
    colnames(ndvi_cm6)  <- seq(1982,2014,1)
    fnout <- "./CMIP6_NDVImax.csv"
    write.csv(ndvi_cm6, fnout)



 




