## Packages ################################
require(TropFishR)
require(LBSPR)
require(LIME)
require(fishdynr)
require(fishmethods)





## True Parameter Input ################################
# setwd("D:/DPLBM/Pub1/Longlived")
# Longmodel <- readRDS("D:/DPLBM/Pub1/Longlived/LFQlongmodel.rds")
Longmodel <- readRDS("LFQlongmodel.rds")
list <- readRDS("Longmodel_list.rds")
iters <- 300

Longmodel1 <- Longmodel
# max(list$inds[[1]][[1]][,5])
# list$inds[[1]][[1]][,6][max(list$inds[[1]][[1]][,5])]


for(i in 1:iters){
  Longmodel1[[i]]$Linf <- 90
  Longmodel1[[i]]$K <- 0.13
  Longmodel1[[i]]$t0 <- -0.01
  Longmodel1[[i]]$M <- 0.18
  Longmodel1[[i]]$Lmat <- 50
  Longmodel1[[i]]$tmax <- 26
}


setwd("D:/DPLBM/Pub1/Longlived/True")
saveRDS(Longmodel1, file = "LFQlongmodel.rds")








## Thompson and Bell ################################
modpath = "D:/DPLBM/Pub1/Longlived/True"
setwd(modpath)
LFQlongmodel <- readRDS("LFQlongmodel.rds")
iters <- 300


#' LFQ conversions
LFQlongmodel1 <- LFQlongmodel
for (i in 1:iters){
  LFQlongmodel1[[i]] <- lfqModify(LFQlongmodel1[[i]], bin_size = 3)
  LFQlongmodel1[[i]] <- lfqModify(LFQlongmodel1[[i]], vectorise_catch = T)
}

#' CC analysis
LFQlongmodel2 <- LFQlongmodel1
res_cc <- list()
for (i in 1:iters){
  res_cc[[i]] <- catchCurve(LFQlongmodel2[[i]], auto = TRUE, calc_ogive = T)
  res_cc[[i]]$wqs <- (res_cc[[i]]$L75 - res_cc[[i]]$L50)*2
}


LFQlongmodel3 <- LFQlongmodel2

for(i in 1:iters){
  LFQlongmodel3[[i]] <- lfqModify(LFQlongmodel3[[i]], plus_group = "Linf")
  LFQlongmodel3[[i]]$a <- 0.0123
  LFQlongmodel3[[i]]$b <- 3.035
  LFQlongmodel3[[i]]$FM <- res_cc[[i]]$Z - LFQlongmodel3[[i]]$M  
  LFQlongmodel3[[i]]$E <- LFQlongmodel3[[i]]$FM / res_cc[[i]]$Z 
  LFQlongmodel3[[i]]$wmat <- 0.15 * LFQlongmodel3[[i]]$Lmat # default definition of wmat
}

## selectivity from CC
for(i in 1:iters){
  LFQlongmodel3[[i]]$FM <- LFQlongmodel3[[i]]$FM * logisticSelect(LFQlongmodel3[[i]]$midLengths,
                                                                    res_cc[[i]]$L50, res_cc[[i]]$wqs)
}



#' TB
res_TB2 <- list()
for (i in 1:iters){
  res_TB2[[i]] <- predict_mod(LFQlongmodel3[[i]],
                              type = "ThompBell",
                              FM_change = seq(0, 2.5, 0.01),
                              stock_size_1 = 1,
                              curr.E = LFQlongmodel3[[i]]$E,
                              plot = FALSE,
                              hide.progressbar = TRUE)
}


saveRDS(res_cc, file = "Longmodel_res_cc.rds")
saveRDS(res_TB2, file = "Longmodel_res_TB.rds")

rm(list = ls())



## Length-Based Risk Analysis ####
## functions

calc_abund_SB <- function(ages, Sl_a, M, FM, R0){
  
  N_a <- R0
  
  Fmat <- NA
  for(i in 1:length(Sl_a)){
    Fmat[i] <- Sl_a[i] * FM
  }
  
  dt <- diff(ages)
  
  for (i in 2:length(ages)) {
    if (i < length(ages))
      N_a[i] <- N_a[i - 1] * exp(-(M + Fmat[i-1])*dt[i - 1])
    if (i == length(ages))
      N_a[i] <- N_a[i - 1] * exp(-(M + Fmat[i-1])*dt[i - 1])/(1 - exp(-(M + Fmat[i-1])*dt[i - 1]))
  }
  return(N_a)
}


SPR_SB <- function(ages, Sl_a, Mat_a, W_a, M, FM){
  
  dt <- diff(ages)[1]
  
  Na0 <- calc_abund_SB(ages = ages, Sl_a, M = M, FM = 0, R0 = R0)
  Naf <- calc_abund_SB(ages = ages, Sl_a, M = M, FM = FM, R0 = R0)
  
  # fished and unfished
  SB0 <- sum(Na0*Mat_a*W_a)*dt
  SBf <- sum(Naf*Mat_a*W_a)*dt
  
  # Compute and return SPR value
  SPR <- SBf/SB0
}

YPR_SB <- function(FM, ages, M, R0, W_a, Sl_a) {
  
  dt <- diff(ages)[1]
  
  Nage <- calc_abund_SB(ages = ages, Sl_a = Sl_a,  M = M, FM = FM, R0 = R0)
  YPR <- FM * sum(Nage * W_a * Sl_a) * dt
  
  return(YPR)
}

SSB_SB <- function(ages, Sl_a, Mat_a, W_a, M, FM){
  
  dt <- diff(ages)[1]
  
  Nage <- calc_abund_SB(ages = ages, Sl_a = Sl_a, M = M, FM = FM, R0 = R0)
  SSB <- sum(Nage * Mat_a * W_a) * dt
  
  return(SSB)
}


modpath <- "D:/DPLBM/Pub1/Longlived/True"
setwd(modpath)
LFQlongmodel <- readRDS("LFQlongmodel.rds")
iters <- 300


#' BH eq - fishmethods
Longlived_res_bheq <- list()
LFQlongmodel1 <- LFQlongmodel 

len <- list()
for(i in 1:iters){ 
  len[[i]] <- rep(LFQlongmodel1[[i]]$midLengths, rowSums(LFQlongmodel1[[i]]$catch))
}

max <- list()
for (i in 1:iters){
  max[[i]] <- max(LFQlongmodel1[[i]]$midLengths)
  if (max[[i]] > LFQlongmodel1[[i]]$Linf){
    max[[i]] = LFQlongmodel1[[i]]$Linf
  }
}

res_cc <- readRDS("Longmodel_res_cc.rds")
for(i in 1:iters){
  LFQlongmodel1[[i]]$L50 <- res_cc[[i]]$L50
}

set.seed(1)

source("~/DPLBM/bheq_mod.R")

for(i in 1:iters){
  Longlived_res_bheq[[i]] <- bheq_LBRA(len[[i]], type = 2, K = LFQlongmodel1[[i]]$K, Linf = LFQlongmodel1[[i]]$Linf, Lc = LFQlongmodel1[[i]]$L50, La = max[[i]], nboot = 200)
}

saveRDS(Longlived_res_bheq, file = "Longmodel_res_bheq.rds")


## beta --> find M
## 1/beta = -log(0.001)/delta_a_lambda
## delta_a_lambda = theor_a - obs_a
del_a <- NA
dev <- NA
for(i in 1:iters){
  del_a[i] <-  LFQlongmodel1[[i]]$tmax - ceiling(-log(0.001)/LFQlongmodel1[[i]]$M)
  dev[i] <- ceiling(-log(0.001)/LFQlongmodel1[[i]]$M) + (del_a[i] / -log(0.001))
}



lwa = 0.0123
lwb = 3.035
t0 <- -0.01
R0 <- 1
tincr <- 1/12 # monthly
binwidth <- 3
CV <- 0.07


linf <- list()
vbk <- list()
ages <- list()
L_a <- list()
W_a <- list()
Mat_a <- list()
M <- list()
FM <- list()
Sl_a <- list()
mids <- list()
highs <- list()
lows <- list()
plba_a <- list()
Mat <- list()

lbprobs <- function(mnl, sdl) return(pnorm(highs[[i]], mnl, sdl) - pnorm(lows[[i]], mnl, sdl))
vlprobs <- Vectorize(lbprobs, vectorize.args = c("mnl", "sdl"))

for (i in 1:iters){
  linf[[i]] <- LFQlongmodel1[[i]]$Linf
  vbk[[i]] <- LFQlongmodel1[[i]]$K
  ages[[i]] <- seq(0, -log(0.001)/LFQlongmodel1[[i]]$M, tincr)
  L_a[[i]] <- linf[[i]]*(1-exp(-vbk[[i]]*(ages[[i]] - t0)))
  
  mids[[i]] <- seq((binwidth/2), linf[[i]]*1.5, binwidth)
  highs[[i]] <- mids[[i]] + (binwidth/2)
  lows[[i]] <- mids[[i]] - (binwidth/2)
  plba_a[[i]] <- t(vlprobs(L_a[[i]], L_a[[i]] * CV))
  plba_a[[i]] <- plba_a[[i]] / rowSums(plba_a[[i]])
  
  W_a[[i]] <- lwa*L_a[[i]]^lwb
  
  M[[i]] <- -log(0.001)/dev[i]
  # M[[i]] <- LFQlongmodel1[[i]]$M
  
  # instead of knife-edge --> logistic
  Mat_a[[i]] <- 1 / (1 + exp(-(mids[[i]] - LFQlongmodel1[[i]]$Lmat) /
                               ((LFQlongmodel1[[i]]$Lmat * 0.15) / ( log(0.75/(1-0.75)) - log(0.25/(1-0.25)) ))))
  Mat_a[[i]] <- apply(t(plba_a[[i]])*Mat_a[[i]], 2, sum)
  Sl_a[[i]] <- 1 / (1 + exp(-(mids[[i]] - LFQlongmodel1[[i]]$L50) /
                              (res_cc[[i]]$wqs / ( log(0.75/(1-0.75)) - log(0.25/(1-0.25)) ))))
  Sl_a[[i]] <- apply(t(plba_a[[i]])*Sl_a[[i]], 2, sum)
}

for(i in 1:iters){
  FM[[i]] <- Longlived_res_bheq[[i]]$z - M[[i]]
  if(FM[[i]] < 0){
    FM[[i]] = 0
  }
}


#' run model
Longlived_res_SB <- list(SPR = NA, YPR = NA, Fmsy = NA, Bmsy = NA, FFmsy = NA, BBmsy = NA, FM = NA)
for (i in 1:iters){
  Longlived_res_SB$SPR[i] <- SPR_SB(ages[[i]], Sl_a[[i]], Mat_a[[i]], W_a[[i]], M[[i]], FM[[i]])
  Longlived_res_SB$YPR[i] <- YPR_SB(FM[[i]], ages[[i]], M[[i]], R0, W_a[[i]], Sl_a[[i]])
  Longlived_res_SB$Fmsy[i] <- optimize(YPR_SB, ages = ages[[i]], M = M[[i]], R0 = R0, W_a = W_a[[i]], Sl_a = Sl_a[[i]], lower = 0, upper = 10, maximum = TRUE)$maximum
  Longlived_res_SB$Bmsy[i] <- SSB_SB(ages[[i]], Sl_a[[i]], Mat_a[[i]], W_a[[i]], M[[i]], Longlived_res_SB$Fmsy[i])
  Longlived_res_SB$FFmsy[i] <- FM[[i]]/Longlived_res_SB$Fmsy[i]
  Longlived_res_SB$BBmsy[i] <- SSB_SB(ages[[i]], Sl_a[[i]], Mat_a[[i]], W_a[[i]], M[[i]], FM[[i]]) / Longlived_res_SB$Bmsy[i]
  Longlived_res_SB$FM[i] <- FM[[i]]
  Longlived_res_SB$SPRmsy[i] <- SPR_SB(ages[[i]], Sl_a[[i]], Mat_a[[i]], W_a[[i]], M[[i]], Longlived_res_SB$Fmsy[i])
}


#' save SB
saveRDS(Longlived_res_SB, file = "Longmodel_res_LBRA.rds")

rm(list = ls())





## LBSPR ####
modpath <- "D:/DPLBM/Pub1/Longlived/True"
setwd(modpath)
LFQlongmodel <- readRDS("LFQlongmodel.rds")
iters <- 300

LFQlongmodel1 <- LFQlongmodel
a = 0.0123
b = 3.035


for(i in 1:iters){
  LFQlongmodel1[[i]] <- lfqModify(LFQlongmodel1[[i]], bin_size = 3)
}


pg <- list()
for(i in 1:iters){
  if(max(LFQlongmodel1[[i]]$midLengths) < LFQlongmodel1[[i]]$Linf){
    closest<-function(xv,sv){
      xv[which(abs(xv-sv)==min(abs(xv-sv)))] 
    }
    
    k <- seq(min(LFQlongmodel1[[i]]$midLengths),120,3)
    l <- closest(k, LFQlongmodel1[[i]]$Linf+2)
    
    pg[[i]] <- c(LFQlongmodel1[[i]]$midLengths[-length(LFQlongmodel1[[i]]$midLengths)], seq(max(LFQlongmodel1[[i]]$midLengths),l,3))
  }else{
    pg[[i]] <- LFQlongmodel1[[i]]$midLengths
  }
}

catch <- list()
for(i in 1:iters){
  LFQlongmodel1[[i]]$catch <- rowSums(LFQlongmodel1[[i]]$catch)
  catch[[i]] <- c(LFQlongmodel1[[i]]$catch,rep(0,length(LFQlongmodel1[[i]]$midLengths[pg[[i]] > max(LFQlongmodel1[[i]]$midLengths)])))
}


res_cc <- readRDS("Longmodel_res_cc.rds")
for(i in 1:iters){
  LFQlongmodel1[[i]]$SL50 <- res_cc[[i]]$L50
  LFQlongmodel1[[i]]$SL95 <- res_cc[[i]]$L95
}


# fmsy calc
OptYield <- function(logF, LHpars) {
  LHpars@FM <- exp(logF)/LHpars@M
  Sim <- LBSPRsim(LHpars)
  -Sim@Yield
}




#' run model
lbspr_res <- list()
FMSY <- NA
SPRmsy <- list()


#' run model
lbspr_res <- list()
for(i in 1:iters){
  PARSlong <- new("LB_pars")
  PARSlong@Linf <- LFQlongmodel1[[i]]$Linf
  PARSlong@MK <- LFQlongmodel1[[i]]$M/LFQlongmodel1[[i]]$K
  PARSlong@L_units <- "cm"
  PARSlong@L50 <- LFQlongmodel[[i]]$Lmat
  PARSlong@L95 <- LFQlongmodel[[i]]$Lmat * 1.3 # definition
  PARSlong@Walpha <- a
  PARSlong@Wbeta <- b
  PARSlong@BinWidth <- 3
  
  SPRlong <- new("LB_lengths", LB_pars = PARSlong)
  SPRlong@LMids <- pg[[i]]
  SPRlong@LData <- as.matrix(catch[[i]])
  SPRlong@L_units <- "cm"
  SPRlong@Years <- 1
  SPRlong@NYears <- 1
  
  lbspr_res[[i]] <- tryCatch(LBSPRfit(LB_pars = PARSlong, LB_lengths = SPRlong, yrs = 1, Control = list(modtype = "GTG")))
  
  # Fmsy
  PARSlong@Steepness <- 0.99
  PARSlong@SL50 <- lbspr_res[[i]]@Ests[,"SL50"]
  PARSlong@SL95 <- lbspr_res[[i]]@Ests[,"SL95"]
  
  PARSlong@M <- LFQlongmodel1[[i]]$M
  PARSlong@BinMax <- 1.3 * PARSlong@Linf
  PARSlong@BinMin <- 0
  PARSlong@BinWidth <- 3
  
  opt <- optimise(OptYield, interval=log(c(0.001, 0.8)), LHpars= PARSlong)
  
  FMSY[i] <- exp(opt$minimum)
  
  
  # SPR msy
  PARSlong@FM <- FMSY[i]/LFQlongmodel1[[i]]$M
  
  SPRmsy[[i]] <- LBSPRsim(PARSlong)
  
}

#' save data
LBSPR_outs <- list(pLCatch = list(rep(NA, 100)), SL50 = NA, SL95 = NA, FM = NA, SPR = NA, SPR_Var = NA, SL50_Var = NA, SL95_Var = NA, FM_Var = NA)
for (i in 1:iters){
  LBSPR_outs$pLCatch[[i]] <- lbspr_res[[i]]@pLCatch
  LBSPR_outs$SL50[[i]] <- lbspr_res[[i]]@Ests[,"SL50"]
  LBSPR_outs$SL95[[i]] <- lbspr_res[[i]]@Ests[,"SL95"]
  LBSPR_outs$FM[[i]] <- lbspr_res[[i]]@Ests[,"FM"]
  LBSPR_outs$SPR[[i]] <- lbspr_res[[i]]@SPR
  LBSPR_outs$SPR_Var[[i]] <- lbspr_res[[i]]@Vars[,"SPR"]
  LBSPR_outs$SL50_Var[[i]] <- lbspr_res[[i]]@Vars[,"SL50"]
  LBSPR_outs$SL95_Var[[i]] <- lbspr_res[[i]]@Vars[,"SL95"]
  LBSPR_outs$FM_Var[[i]] <- lbspr_res[[i]]@Vars[,"FM"]
  LBSPR_outs$Fmort[[i]] <- lbspr_res[[i]]@Ests[,"FM"] * LFQlongmodel[[i]]$M
  LBSPR_outs$Fmsy[[i]] <- FMSY[i]
  LBSPR_outs$FFmsy[[i]] <- LBSPR_outs$Fmort[[i]]/FMSY[i]
  LBSPR_outs$SPRmsy[[i]] <- SPRmsy[[i]]@SPR
}
saveRDS(LBSPR_outs, file = "Longmodel_res_LBSPR.rds")

rm(list = ls())











## LIME ####
modpath <- "D:/DPLBM/Pub1/Longlived/True"
setwd(modpath)
LFQlongmodel <- readRDS("LFQlongmodel.rds")
iters <- 300

#' add SL50 and 95
LFQlongmodel1 <- LFQlongmodel
a = 0.0123
b = 3.035

res_cc <- readRDS("Longmodel_res_cc.rds")
for(i in 1:iters){
  LFQlongmodel1[[i]]$SL50 <- res_cc[[i]]$L50
  LFQlongmodel1[[i]]$SL95 <- res_cc[[i]]$L95
}

#' Input parameters for lh list
lh <- list()
for(i in 1:iters){
  lh[[i]] <- create_lh_list(vbk = LFQlongmodel1[[i]]$K, 
                            linf = LFQlongmodel1[[i]]$Linf,
                            t0 = -0.01,
                            lwa = a,
                            lwb = b,
                            S50 = LFQlongmodel1[[i]]$SL50,
                            S95 = LFQlongmodel1[[i]]$SL95,
                            selex_input = "length",
                            selex_type = "logistic",
                            M50 = LFQlongmodel1[[i]]$Lmat,
                            M95 = NULL,
                            maturity_input = "length",
                            M = LFQlongmodel1[[i]]$M,
                            binwidth = 3,
                            R0 = 1,
                            nseasons = 12)
}



#' set up data
for(i in 1:iters){
  LFQlongmodel1[[i]] <- lfqModify(LFQlongmodel1[[i]], bin_size = 3)
  LFQlongmodel1[[i]]$catch <- rowSums(LFQlongmodel1[[i]]$catch)
  LFQlongmodel1[[i]]$catch <- as.matrix(LFQlongmodel1[[i]]$catch)
}

#' LF list
LF <- list()
for(i in 1:iters){
  LF[[i]] <- t(LFQlongmodel1[[i]]$catch)
  rownames(LF[[i]]) <- as.numeric(1)
  colnames(LF[[i]]) <- LFQlongmodel1[[i]]$midLengths
}

#' years with length data
data_LF <- list()
for (i in 1:iters){
  data_LF[[i]] <- list("years" = 1, "LF" = LF[[i]])
}

#' Run LIME
Longlived_res_LIME <- list()
data_all <- list()

for(i in 1:iters){
  data_all[[i]] <- create_inputs(lh = lh[[i]], input_data = data_LF[[i]])
}

# a <- c(which(p %in% "The model is likely not converged"))
for (i in 1:300){
  Longlived_res_LIME[[i]] <- run_LIME(modpath = NULL,
                                        input = data_all[[i]],
                                        # lh = lh[[i]],
                                        # input_data = data_LF[[i]],
                                        # est_sigma = "log_sigma_R",
                                        data_avail = "LC",
                                        derive_quants = TRUE
  )
}

# calc_derived_quants(Longlived_res_LIME[[i]]$obj, lh[[i]])

for(i in 1:300){
  Longlived_res_LIME[[i]]$Derived <- Longlived_res_LIME2[[i]]$Derived
}

saveRDS(Longlived_res_LIME, file = "Longmodel_res_LIME.rds")

# devtools::install_github("merrillrudd/LIME")
# devtools::install_github("merrillrudd/LIME@b135ea8d2b9317c3e1a44e90f92a7cd9776938d1")

rm(list = ls())