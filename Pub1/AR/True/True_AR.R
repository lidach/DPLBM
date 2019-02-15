## Packages ################################
require(TropFishR)
require(LBSPR)
require(LIME)
require(fishdynr)
require(fishmethods)



## True Parameter Input ################################
# setwd("D:/DPLBM/Pub1/ARlived")
# ARlived <- readRDS("D:/DPLBM/Pub1/ARlived/LFQARmodel.rds")
ARlived <- readRDS("LFQARmodel.rds")
iters <- 300

ARlived1 <- ARlived

for(i in 1:iters){
  ARlived1[[i]]$Linf <- 64.6
  ARlived1[[i]]$K <- 0.21
  ARlived1[[i]]$t0 <- -0.01
  ARlived1[[i]]$M <- 0.32
  ARlived1[[i]]$Lmat <- 34
  ARlived1[[i]]$tmax <- 18
}


setwd("D:/DPLBM/Pub1/ARlived/True")
saveRDS(ARlived1, file = "LFQARmodel.rds")

rm(list = ls())







## Thompson and Bell ################################
LFQARmodel <- readRDS("LFQARmodel.rds")
iters <- 300

#' LFQ conversions
LFQARmodel1 <- LFQARmodel
for (i in 1:iters){
  LFQARmodel1[[i]] <- lfqModify(LFQARmodel1[[i]], bin_size = 2)
  LFQARmodel1[[i]] <- lfqModify(LFQARmodel1[[i]], vectorise_catch = T)
}

#' CC analysis
LFQARmodel2 <- LFQARmodel1
res_cc <- list()
for (i in 1:iters){
  res_cc[[i]] <- catchCurve(LFQARmodel2[[i]], auto=TRUE, calc_ogive = TRUE)
  res_cc[[i]]$wqs <- (res_cc[[i]]$L75 - res_cc[[i]]$L50)*2
}

LFQARmodel3 <- LFQARmodel2

for(i in 1:iters){
  LFQARmodel3[[i]] <- lfqModify(LFQARmodel3[[i]], plus_group = "Linf")
  LFQARmodel3[[i]]$a <- 0.018
  LFQARmodel3[[i]]$b <- 2.895
  LFQARmodel3[[i]]$FM <- res_cc[[i]]$Z - LFQARmodel3[[i]]$M  
  LFQARmodel3[[i]]$E <- LFQARmodel3[[i]]$FM / res_cc[[i]]$Z 
  LFQARmodel3[[i]]$wmat <- 0.15 * LFQARmodel3[[i]]$Lmat # default definition of wmat
}

## selectivity from CC
for(i in 1:iters){
  LFQARmodel3[[i]]$FM <- LFQARmodel3[[i]]$FM * logisticSelect(LFQARmodel3[[i]]$midLengths,
                                                                    res_cc[[i]]$L50, res_cc[[i]]$wqs)
}



#' TB
res_TB2 <- list()
for (i in 1:iters){
  res_TB2[[i]] <- predict_mod(LFQARmodel3[[i]],
                              type = "ThompBell",
                              FM_change = seq(0, 2.5, 0.01),
                              stock_size_1 = 1,
                              curr.E = LFQARmodel3[[i]]$E,
                              plot = FALSE,
                              hide.progressbar = TRUE)
}

saveRDS(res_cc, file = "ARmodel_res_cc.rds")
saveRDS(res_TB2, file = "ARmodel_res_TB.rds")
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


modpath <- "D:/DPLBM/Pub1/AR/True"
setwd(modpath)
LFQARmodel <- readRDS("LFQARmodel.rds")
iters <- 300


#' BH eq - fishmethods
ARlived_res_bheq <- list()
LFQARmodel1 <- LFQARmodel 

len <- list()
for(i in 1:iters){ 
  len[[i]] <- rep(LFQARmodel1[[i]]$midLengths, rowSums(LFQARmodel1[[i]]$catch))
}

max <- list()
for (i in 1:iters){
  max[[i]] <- max(LFQARmodel1[[i]]$midLengths)
  if (max[[i]] > LFQARmodel1[[i]]$Linf){
    max[[i]] = LFQARmodel1[[i]]$Linf
  }
}

res_cc <- readRDS("ARmodel_res_cc.rds")
for(i in 1:iters){
  LFQARmodel1[[i]]$L50 <- res_cc[[i]]$L50
}

set.seed(1)

source("~/DPLBM/bheq_mod.R")

for(i in 1:iters){
  ARlived_res_bheq[[i]] <- bheq_LBRA(len[[i]], type = 2, K = LFQARmodel1[[i]]$K, Linf = LFQARmodel1[[i]]$Linf, Lc = LFQARmodel1[[i]]$L50, La = max[[i]], nboot = 200)
}

saveRDS(ARlived_res_bheq, file = "ARmodel_res_bheq.rds")


## beta --> find M
## 1/beta = -log(0.001)/delta_a_lambda
## delta_a_lambda = theor_a - obs_a
del_a <- NA
dev <- NA
for(i in 1:iters){
  del_a[i] <-  LFQARmodel1[[i]]$tmax - ceiling(-log(0.001)/LFQARmodel1[[i]]$M)
  dev[i] <- ceiling(-log(0.001)/LFQARmodel1[[i]]$M) + (del_a[i] / -log(0.001))
}



lwa <- 0.018
lwb <- 2.895
t0 <- -0.01
R0 <- 1
tincr <- 1/12 # monthly
binwidth <- 2
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
  linf[[i]] <- LFQARmodel1[[i]]$Linf
  vbk[[i]] <- LFQARmodel1[[i]]$K
  ages[[i]] <- seq(0, LFQARmodel1[[i]]$tmax, tincr)
  L_a[[i]] <- linf[[i]]*(1-exp(-vbk[[i]]*(ages[[i]] - t0)))
  
  mids[[i]] <- seq((binwidth/2), linf[[i]]*1.5, binwidth)
  highs[[i]] <- mids[[i]] + (binwidth/2)
  lows[[i]] <- mids[[i]] - (binwidth/2)
  plba_a[[i]] <- t(vlprobs(L_a[[i]], L_a[[i]] * CV))
  plba_a[[i]] <- plba_a[[i]] / rowSums(plba_a[[i]])
  
  W_a[[i]] <- lwa*L_a[[i]]^lwb
  
  M[[i]] <- -log(0.001)/dev[i]
  # M[[i]] <- LFQARmodel1[[i]]$M
  
  # instead of knife-edge --> logistic
  Mat_a[[i]] <- 1 / (1 + exp(-(mids[[i]] - LFQARmodel1[[i]]$Lmat) /
                               ((LFQARmodel1[[i]]$Lmat * 0.15) / ( log(0.75/(1-0.75)) - log(0.25/(1-0.25)) ))))
  Mat_a[[i]] <- apply(t(plba_a[[i]])*Mat_a[[i]], 2, sum)
  Sl_a[[i]] <- 1 / (1 + exp(-(mids[[i]] - LFQARmodel1[[i]]$L50) /
                              (res_cc[[i]]$wqs / ( log(0.75/(1-0.75)) - log(0.25/(1-0.25)) ))))
  Sl_a[[i]] <- apply(t(plba_a[[i]])*Sl_a[[i]], 2, sum)
}

for(i in 1:iters){
  FM[[i]] <- ARlived_res_bheq[[i]]$z - M[[i]]
  if(FM[[i]] < 0){
    FM[[i]] = 0
  }
}


#' run model
ARlived_res_SB <- list(SPR = NA, YPR = NA, Fmsy = NA, Bmsy = NA, FFmsy = NA, BBmsy = NA, FM = NA)
for (i in 1:iters){
  ARlived_res_SB$SPR[i] <- SPR_SB(ages[[i]], Sl_a[[i]], Mat_a[[i]], W_a[[i]], M[[i]], FM[[i]])
  ARlived_res_SB$YPR[i] <- YPR_SB(FM[[i]], ages[[i]], M[[i]], R0, W_a[[i]], Sl_a[[i]])
  ARlived_res_SB$Fmsy[i] <- optimize(YPR_SB, ages = ages[[i]], M = M[[i]], R0 = R0, W_a = W_a[[i]], Sl_a = Sl_a[[i]], lower = 0, upper = 10, maximum = TRUE)$maximum
  ARlived_res_SB$Bmsy[i] <- SSB_SB(ages[[i]], Sl_a[[i]], Mat_a[[i]], W_a[[i]], M[[i]], ARlived_res_SB$Fmsy[i])
  ARlived_res_SB$FFmsy[i] <- FM[[i]]/ARlived_res_SB$Fmsy[i]
  ARlived_res_SB$BBmsy[i] <- SSB_SB(ages[[i]], Sl_a[[i]], Mat_a[[i]], W_a[[i]], M[[i]], FM[[i]]) / ARlived_res_SB$Bmsy[i]
  ARlived_res_SB$FM[i] <- FM[[i]]
  ARlived_res_SB$SPRmsy[i] <- SPR_SB(ages[[i]], Sl_a[[i]], Mat_a[[i]], W_a[[i]], M[[i]], ARlived_res_SB$Fmsy[i])
}


#' save SB
saveRDS(ARlived_res_SB, file = "ARmodel_res_LBRA.rds")

rm(list = ls())





## LBSPR ####
modpath <- "D:/DPLBM/Pub1/AR/True"
setwd(modpath)
LFQARmodel <- readRDS("LFQARmodel.rds")
iters <- 300

LFQARmodel1 <- LFQARmodel
a <- 0.018
b <- 2.895


for(i in 1:iters){
  LFQARmodel1[[i]] <- lfqModify(LFQARmodel1[[i]], bin_size = 2)
}


pg <- list()
for(i in 1:iters){
  if(max(LFQARmodel1[[i]]$midLengths) < LFQARmodel1[[i]]$Linf){
    closest<-function(xv,sv){
      xv[which(abs(xv-sv)==min(abs(xv-sv)))] 
    }
    
    k <- seq(min(LFQARmodel1[[i]]$midLengths),120,2)
    l <- closest(k, LFQARmodel1[[i]]$Linf+2)
    
    pg[[i]] <- c(LFQARmodel1[[i]]$midLengths[-length(LFQARmodel1[[i]]$midLengths)], seq(max(LFQARmodel1[[i]]$midLengths),l,2))
  }else{
    pg[[i]] <- LFQARmodel1[[i]]$midLengths
  }
}

catch <- list()
for(i in 1:iters){
  LFQARmodel1[[i]]$catch <- rowSums(LFQARmodel1[[i]]$catch)
  catch[[i]] <- c(LFQARmodel1[[i]]$catch,rep(0,length(LFQARmodel1[[i]]$midLengths[pg[[i]] > max(LFQARmodel1[[i]]$midLengths)])))
}


res_cc <- readRDS("ARmodel_res_cc.rds")
for(i in 1:iters){
  LFQARmodel1[[i]]$SL50 <- res_cc[[i]]$L50
  LFQARmodel1[[i]]$SL95 <- res_cc[[i]]$L95
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
  PARSAR <- new("LB_pars")
  PARSAR@Linf <- LFQARmodel1[[i]]$Linf
  PARSAR@MK <- LFQARmodel1[[i]]$M/LFQARmodel1[[i]]$K
  PARSAR@L_units <- "cm"
  PARSAR@L50 <- LFQARmodel[[i]]$Lmat
  PARSAR@L95 <- LFQARmodel[[i]]$Lmat * 1.3 # definition
  PARSAR@Walpha <- a
  PARSAR@Wbeta <- b
  PARSAR@BinWidth <- 2
  
  SPRAR <- new("LB_lengths", LB_pars = PARSAR)
  SPRAR@LMids <- pg[[i]]
  SPRAR@LData <- as.matrix(catch[[i]])
  SPRAR@L_units <- "cm"
  SPRAR@Years <- 1
  SPRAR@NYears <- 1
  
  lbspr_res[[i]] <- tryCatch(LBSPRfit(LB_pars = PARSAR, LB_lengths = SPRAR, yrs = 1, Control = list(modtype = "GTG")))
  
  # Fmsy
  PARSAR@Steepness <- 0.99
  PARSAR@SL50 <- lbspr_res[[i]]@Ests[,"SL50"]
  PARSAR@SL95 <- lbspr_res[[i]]@Ests[,"SL95"]
  
  PARSAR@M <- LFQARmodel1[[i]]$M
  PARSAR@BinMax <- 1.3 * PARSAR@Linf
  PARSAR@BinMin <- 0
  PARSAR@BinWidth <- 2
  
  opt <- optimise(OptYield, interval=log(c(0.001, 0.8)), LHpars= PARSAR)
  
  FMSY[i] <- exp(opt$minimum)
  
  
  # SPR msy
  PARSAR@FM <- FMSY[i]/LFQARmodel1[[i]]$M
  
  SPRmsy[[i]] <- LBSPRsim(PARSAR)
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
  LBSPR_outs$Fmort[[i]] <- lbspr_res[[i]]@Ests[,"FM"] * LFQARmodel[[i]]$M
  LBSPR_outs$Fmsy[[i]] <- FMSY[i]
  LBSPR_outs$FFmsy[[i]] <- LBSPR_outs$Fmort[[i]]/FMSY[i]
  LBSPR_outs$SPRmsy[[i]] <- SPRmsy[[i]]@SPR
}
saveRDS(LBSPR_outs, file = "ARmodel_res_LBSPR.rds")

rm(list = ls())








## LIME ####
modpath <- "D:/DPLBM/Pub1/ARlived/True"
setwd(modpath)
LFQARmodel <- readRDS("LFQARmodel.rds")
iters <- 300

#' add SL50 and 95
LFQARmodel1 <- LFQARmodel
a = 0.018
b = 2.895

res_cc <- readRDS("ARmodel_res_cc.rds")
for(i in 1:iters){
  LFQARmodel1[[i]]$SL50 <- res_cc[[i]]$L50
  LFQARmodel1[[i]]$SL95 <- res_cc[[i]]$L95
}

#' Input parameters for lh list
lh <- list()
for(i in 1:iters){
  lh[[i]] <- create_lh_list(vbk = LFQARmodel1[[i]]$K, 
                            linf = LFQARmodel1[[i]]$Linf,
                            t0 = -0.01,
                            lwa = a,
                            lwb = b,
                            S50 = LFQARmodel1[[i]]$SL50,
                            S95 = LFQARmodel1[[i]]$SL95,
                            selex_input = "length",
                            selex_type = "logistic",
                            M50 = LFQARmodel1[[i]]$Lmat,
                            M95 = NULL,
                            maturity_input = "length",
                            M = LFQARmodel1[[i]]$M,
                            binwidth = 2,
                            R0 = 1,
                            nseasons = 1)
}



#' set up data
for(i in 1:iters){
  LFQARmodel1[[i]] <- lfqModify(LFQARmodel1[[i]], bin_size = 2)
  LFQARmodel1[[i]]$catch <- rowSums(LFQARmodel1[[i]]$catch)
  LFQARmodel1[[i]]$catch <- as.matrix(LFQARmodel1[[i]]$catch)
}

#' LF list
LF <- list()
for(i in 1:iters){
  LF[[i]] <- t(LFQARmodel1[[i]]$catch)
  rownames(LF[[i]]) <- as.numeric(1)
  colnames(LF[[i]]) <- LFQARmodel1[[i]]$midLengths
}

#' years with length data
data_LF <- list()
for (i in 1:iters){
  data_LF[[i]] <- list("years" = 1, "LF" = LF[[i]])
}

#' Run LIME
ARlived_res_LIME <- list()
# data_all <- list()
# 
# for(i in 1:iters){
#   data_all[[i]] <- create_inputs(lh = lh[[i]], input_data = data_LF[[i]])
# }

# a <- c(which(p %in% "The model is likely not converged"))
for (i in 1:300){
  ARlived_res_LIME[[i]] <- run_LIME(modpath = NULL,
                                       # input = data_all[[i]],
                                       lh = lh[[i]],
                                       input_data = data_LF[[i]],
                                       est_sigma = "log_sigma_R",
                                       data_avail = "LC"
                                       # derive_quants = TRUE
  )
}

# calc_derived_quants(ARlived_res_LIME[[i]]$obj, lh[[i]])

for(i in 1:300){
  ARlived_res_LIME[[i]]$Derived <- ARlived_res_LIME2[[i]]$Derived
}

saveRDS(ARlived_res_LIME, file = "ARmodel_res_LIME.rds")

# devtools::install_github("merrillrudd/LIME")
# devtools::install_github("merrillrudd/LIME@b135ea8d2b9317c3e1a44e90f92a7cd9776938d1")

rm(list = ls())
