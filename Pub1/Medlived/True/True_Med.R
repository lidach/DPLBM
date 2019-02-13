## Packages ################################
require(TropFishR)
require(LBSPR)
require(LIME)
require(fishdynr)
require(fishmethods)






## True Parameter Input ################################
# setwd("D:/DPLBM/Pub1/Medlived")
# Medlived <- readRDS("D:/DPLBM/Pub1/Medlived/LFQmedmodel.rds")
Medlived <- readRDS("LFQmedmodel.rds")
list <- readRDS("Medlist.rds")
iters <- 300

Medlived1 <- Medlived

for(i in 1:iters){
  Medlived1[[i]]$Linf <- 64.6
  Medlived1[[i]]$K <- 0.21
  Medlived1[[i]]$t0 <- -0.01
  Medlived1[[i]]$M <- 0.32
  Medlived1[[i]]$Lmat <- 34
  Medlived1[[i]]$tmax <- 18
}


setwd("D:/DPLBM/Pub1/Medlived/True")
saveRDS(Medlived1, file = "LFQmedmodel.rds")

rm(list = ls())







## Thompson and Bell ################################
LFQmedmodel <- readRDS("LFQmedmodel.rds")

#' LFQ conversions
LFQmedmodel1 <- LFQmedmodel
for (i in 1:iters){
  LFQmedmodel1[[i]] <- lfqModify(LFQmedmodel1[[i]], bin_size = 2)
  LFQmedmodel1[[i]] <- lfqModify(LFQmedmodel1[[i]], vectorise_catch = T)
}

#' CC analysis
LFQmedmodel2 <- LFQmedmodel1
res_cc <- list()
for (i in 1:iters){
  res_cc[[i]] <- catchCurve(LFQmedmodel2[[i]], auto=TRUE, calc_ogive = TRUE)
  res_cc[[i]]$wqs <- (res_cc[[i]]$L75 - res_cc[[i]]$L50)*2
}

LFQmedmodel3 <- LFQmedmodel2

for(i in 1:iters){
  LFQmedmodel3[[i]] <- lfqModify(LFQmedmodel3[[i]], plus_group = "Linf")
  LFQmedmodel3[[i]]$a <- 0.018
  LFQmedmodel3[[i]]$b <- 2.895
  LFQmedmodel3[[i]]$FM <- res_cc[[i]]$Z - LFQmedmodel3[[i]]$M  
  LFQmedmodel3[[i]]$E <- LFQmedmodel3[[i]]$FM / res_cc[[i]]$Z 
  LFQmedmodel3[[i]]$wmat <- 0.15 * LFQmedmodel3[[i]]$Lmat # default definition of wmat
}

## selectivity from CC
for(i in 1:iters){
  LFQmedmodel3[[i]]$FM <- LFQmedmodel3[[i]]$FM * logisticSelect(LFQmedmodel3[[i]]$midLengths,
                                                                      res_cc[[i]]$L50, res_cc[[i]]$wqs)
}



#' TB
res_TB2 <- list()
for (i in 1:iters){
  res_TB2[[i]] <- predict_mod(LFQmedmodel3[[i]],
                              type = "ThompBell",
                              FM_change = seq(0, 2.5, 0.01),
                              stock_size_1 = 1,
                              curr.E = LFQmedmodel3[[i]]$E,
                              plot = FALSE,
                              hide.progressbar = TRUE)
}

saveRDS(res_cc, file = "Medmodel_res_cc.rds")
saveRDS(res_TB2, file = "Medmodel_res_TB.rds")
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


modpath <- "D:/DPLBM/Pub1/Medlived/True"
setwd(modpath)
LFQmedmodel <- readRDS("LFQmedmodel.rds")
iters <- 300


#' BH eq - fishmethods
Medlived_res_bheq <- list()
LFQmedmodel1 <- LFQmedmodel 

len <- list()
for(i in 1:iters){ 
  len[[i]] <- rep(LFQmedmodel1[[i]]$midLengths, rowSums(LFQmedmodel1[[i]]$catch))
  }

max <- list()
for (i in 1:iters){
  max[[i]] <- max(LFQmedmodel1[[i]]$midLengths)
  if (max[[i]] > LFQmedmodel1[[i]]$Linf){
    max[[i]] = LFQmedmodel1[[i]]$Linf
  }
}

res_cc <- readRDS("Medmodel_res_cc.rds")
for(i in 1:iters){
  LFQmedmodel1[[i]]$L50 <- res_cc[[i]]$L50
}

set.seed(1)

source("~/DPLBM/bheq_mod.R")

for(i in 1:iters){
  Medlived_res_bheq[[i]] <- bheq_LBRA(len[[i]], type = 2, K = LFQmedmodel1[[i]]$K, Linf = LFQmedmodel1[[i]]$Linf, Lc = LFQmedmodel1[[i]]$L50, La = max[[i]], nboot = 200)
}

saveRDS(Medlived_res_bheq, file = "Medmodel_res_bheq.rds")


## beta --> find M
## 1/beta = -log(0.001)/delta_a_lambda
## delta_a_lambda = theor_a - obs_a
del_a <- NA
dev <- NA
for(i in 1:iters){
  del_a[i] <-  LFQmedmodel1[[i]]$tmax - ceiling(-log(0.001)/LFQmedmodel1[[i]]$M)
  dev[i] <- ceiling(-log(0.001)/LFQmedmodel1[[i]]$M) + (del_a[i] / -log(0.001))
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
  linf[[i]] <- LFQmedmodel1[[i]]$Linf
  vbk[[i]] <- LFQmedmodel1[[i]]$K
  ages[[i]] <- seq(0, LFQmedmodel1[[i]]$tmax, tincr)
  L_a[[i]] <- linf[[i]]*(1-exp(-vbk[[i]]*(ages[[i]] - t0)))
  
  mids[[i]] <- seq((binwidth/2), linf[[i]]*1.5, binwidth)
  highs[[i]] <- mids[[i]] + (binwidth/2)
  lows[[i]] <- mids[[i]] - (binwidth/2)
  plba_a[[i]] <- t(vlprobs(L_a[[i]], L_a[[i]] * CV))
  plba_a[[i]] <- plba_a[[i]] / rowSums(plba_a[[i]])
  
  W_a[[i]] <- lwa*L_a[[i]]^lwb
  
  M[[i]] <- -log(0.001)/dev[i]
  # M[[i]] <- LFQmedmodel1[[i]]$M
  
  # instead of knife-edge --> logistic
  Mat_a[[i]] <- 1 / (1 + exp(-(mids[[i]] - LFQmedmodel1[[i]]$Lmat) /
                               ((LFQmedmodel1[[i]]$Lmat * 0.15) / ( log(0.75/(1-0.75)) - log(0.25/(1-0.25)) ))))
  Mat_a[[i]] <- apply(t(plba_a[[i]])*Mat_a[[i]], 2, sum)
  Sl_a[[i]] <- 1 / (1 + exp(-(mids[[i]] - LFQmedmodel1[[i]]$L50) /
                              (res_cc[[i]]$wqs / ( log(0.75/(1-0.75)) - log(0.25/(1-0.25)) ))))
  Sl_a[[i]] <- apply(t(plba_a[[i]])*Sl_a[[i]], 2, sum)
}

for(i in 1:iters){
  FM[[i]] <- Medlived_res_bheq[[i]]$z - M[[i]]
  if(FM[[i]] < 0){
    FM[[i]] = 0
  }
}


#' run model
Medlived_res_SB <- list(SPR = NA, YPR = NA, Fmsy = NA, Bmsy = NA, FFmsy = NA, BBmsy = NA, FM = NA)
for (i in 1:iters){
  Medlived_res_SB$SPR[i] <- SPR_SB(ages[[i]], Sl_a[[i]], Mat_a[[i]], W_a[[i]], M[[i]], FM[[i]])
  Medlived_res_SB$YPR[i] <- YPR_SB(FM[[i]], ages[[i]], M[[i]], R0, W_a[[i]], Sl_a[[i]])
  Medlived_res_SB$Fmsy[i] <- optimize(YPR_SB, ages = ages[[i]], M = M[[i]], R0 = R0, W_a = W_a[[i]], Sl_a = Sl_a[[i]], lower = 0, upper = 10, maximum = TRUE)$maximum
  Medlived_res_SB$Bmsy[i] <- SSB_SB(ages[[i]], Sl_a[[i]], Mat_a[[i]], W_a[[i]], M[[i]], Medlived_res_SB$Fmsy[i])
  Medlived_res_SB$FFmsy[i] <- FM[[i]]/Medlived_res_SB$Fmsy[i]
  Medlived_res_SB$BBmsy[i] <- SSB_SB(ages[[i]], Sl_a[[i]], Mat_a[[i]], W_a[[i]], M[[i]], FM[[i]]) / Medlived_res_SB$Bmsy[i]
  Medlived_res_SB$FM[i] <- FM[[i]]
  Medlived_res_SB$SPRmsy[i] <- SPR_SB(ages[[i]], Sl_a[[i]], Mat_a[[i]], W_a[[i]], M[[i]], Medlived_res_SB$Fmsy[i])
}


#' save SB
saveRDS(Medlived_res_SB, file = "Medmodel_res_LBRA.rds")

rm(list = ls())





## LBSPR ####
modpath <- "D:/DPLBM/Pub1/Medlived/True"
setwd(modpath)
LFQmedmodel <- readRDS("LFQmedmodel.rds")
iters <- 300

LFQmedmodel1 <- LFQmedmodel
a <- 0.018
b <- 2.895


for(i in 1:iters){
  LFQmedmodel1[[i]] <- lfqModify(LFQmedmodel1[[i]], bin_size = 2)
}


pg <- list()
for(i in 1:iters){
  if(max(LFQmedmodel1[[i]]$midLengths) < LFQmedmodel1[[i]]$Linf){
    closest<-function(xv,sv){
      xv[which(abs(xv-sv)==min(abs(xv-sv)))] 
    }
    
    k <- seq(min(LFQmedmodel1[[i]]$midLengths),120,2)
    l <- closest(k, LFQmedmodel1[[i]]$Linf+2)
    
    pg[[i]] <- c(LFQmedmodel1[[i]]$midLengths[-length(LFQmedmodel1[[i]]$midLengths)], seq(max(LFQmedmodel1[[i]]$midLengths),l,2))
  }else{
    pg[[i]] <- LFQmedmodel1[[i]]$midLengths
  }
}

catch <- list()
for(i in 1:iters){
  LFQmedmodel1[[i]]$catch <- rowSums(LFQmedmodel1[[i]]$catch)
  catch[[i]] <- c(LFQmedmodel1[[i]]$catch,rep(0,length(LFQmedmodel1[[i]]$midLengths[pg[[i]] > max(LFQmedmodel1[[i]]$midLengths)])))
}


res_cc <- readRDS("Medmodel_res_cc.rds")
for(i in 1:iters){
  LFQmedmodel1[[i]]$SL50 <- res_cc[[i]]$L50
  LFQmedmodel1[[i]]$SL95 <- res_cc[[i]]$L95
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
for(i in 1:iters){
  PARSmed <- new("LB_pars")
  PARSmed@Linf <- LFQmedmodel1[[i]]$Linf
  PARSmed@MK <- LFQmedmodel1[[i]]$M/LFQmedmodel1[[i]]$K
  PARSmed@L_units <- "cm"
  PARSmed@L50 <- LFQmedmodel[[i]]$Lmat
  PARSmed@L95 <- LFQmedmodel[[i]]$Lmat * 1.3 # definition
  PARSmed@SL50 <- res_cc[[i]]$L50
  PARSmed@SL95 <- res_cc[[i]]$L95
  PARSmed@Walpha <- a
  PARSmed@Wbeta <- b
  PARSmed@BinWidth <- 2
  
  SPRmed <- new("LB_lengths", LB_pars = PARSmed)
  SPRmed@LMids <- pg[[i]]
  SPRmed@LData <- as.matrix(catch[[i]])
  SPRmed@L_units <- "cm"
  SPRmed@Years <- 1
  SPRmed@NYears <- 1
  
  lbspr_res[[i]] <- tryCatch(LBSPRfit(LB_pars = PARSmed, LB_lengths = SPRmed, yrs = 1, Control = list(modtype = "GTG")))
  
  # Fmsy
  PARSmed@Steepness <- 0.99
  PARSmed@SL50 <- lbspr_res[[i]]@Ests[,"SL50"]
  PARSmed@SL95 <- lbspr_res[[i]]@Ests[,"SL95"]
  
  PARSmed@M <- LFQmedmodel1[[i]]$M
  PARSmed@BinMax <- 1.3 * PARSmed@Linf
  PARSmed@BinMin <- 0
  PARSmed@BinWidth <- 2
  
  opt <- optimise(OptYield, interval=log(c(0.001, 0.8)), LHpars= PARSmed)
  
  FMSY[i] <- exp(opt$minimum)
  
  
  # SPR msy
  PARSmed@FM <- FMSY[i]/LFQmedmodel1[[i]]$M
  
  SPRmsy[[i]] <- LBSPRsim(PARSmed)
  
  # FVec <- seq(0, FMSY*2, by=0.01)
  # Yields <- -sapply(log(FVec), OptYield, LHpars=LHpars)/MSY
  # plot(FVec, Yields, xlab="Fishing mortality (apical)", ylab="Yield/MSY", type="l")
  
}

#' save data
LBSPR_outs <- list(pLCatch = list(rep(NA, 100)), SL50 = NA, SL95 = NA, FM = NA, SPR = NA, SPR_Var = NA, SL50_Var = NA, SL95_Var = NA, FM_Var = NA, Fmort = NA, Fmsy = NA, FFmsy = NA, SPRmsy = NA)
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
  LBSPR_outs$Fmort[[i]] <- lbspr_res[[i]]@Ests[,"FM"] * 0.32
  LBSPR_outs$Fmsy[[i]] <- FMSY[i]
  LBSPR_outs$FFmsy[[i]] <- LBSPR_outs$Fmort[[i]]/FMSY[i]
  LBSPR_outs$SPRmsy[[i]] <- SPRmsy[[i]]@SPR
}
saveRDS(LBSPR_outs, file = "Medmodel_res_LBSPR.rds")

rm(list = ls())








## LIME ####
modpath <- "D:/DPLBM/Pub1/Medlived/True"
setwd(modpath)
LFQmedmodel <- readRDS("LFQmedmodel.rds")
iters <- 300

#' add SL50 and 95
LFQmedmodel1 <- LFQmedmodel
a = 0.018
b = 2.895

res_cc <- readRDS("Medmodel_res_cc.rds")
for(i in 1:iters){
  LFQmedmodel1[[i]]$SL50 <- res_cc[[i]]$L50
  LFQmedmodel1[[i]]$SL95 <- res_cc[[i]]$L95
}

#' Input parameters for lh list
lh <- list()
for(i in 1:iters){
  lh[[i]] <- create_lh_list(vbk = LFQmedmodel1[[i]]$K, 
                            linf = LFQmedmodel1[[i]]$Linf,
                            t0 = -0.01,
                            lwa = a,
                            lwb = b,
                            S50 = LFQmedmodel1[[i]]$SL50,
                            S95 = LFQmedmodel1[[i]]$SL95,
                            selex_input = "length",
                            selex_type = "logistic",
                            M50 = LFQmedmodel1[[i]]$Lmat,
                            M95 = NULL,
                            maturity_input = "length",
                            M = LFQmedmodel1[[i]]$M,
                            binwidth = 2,
                            R0 = 1,
                            nseasons = 1)
}



#' set up data
for(i in 1:iters){
  LFQmedmodel1[[i]] <- lfqModify(LFQmedmodel1[[i]], bin_size = 2)
  LFQmedmodel1[[i]]$catch <- rowSums(LFQmedmodel1[[i]]$catch)
  LFQmedmodel1[[i]]$catch <- as.matrix(LFQmedmodel1[[i]]$catch)
}

#' LF list
LF <- list()
for(i in 1:iters){
  LF[[i]] <- t(LFQmedmodel1[[i]]$catch)
  rownames(LF[[i]]) <- as.numeric(1)
  colnames(LF[[i]]) <- LFQmedmodel1[[i]]$midLengths
}

#' years with length data
data_LF <- list()
for (i in 1:iters){
  data_LF[[i]] <- list("years" = 1, "LF" = LF[[i]])
}

#' Run LIME
Medlived_res_LIME <- list()
# data_all <- list()
# 
# for(i in 1:iters){
#   data_all[[i]] <- create_inputs(lh = lh[[i]], input_data = data_LF[[i]])
# }

# a <- c(which(p %in% "The model is likely not converged"))
for (i in 1:300){
  Medlived_res_LIME[[i]] <- run_LIME(modpath = NULL,
                                     # input = data_all[[i]],
                                      lh = lh[[i]],
                                      input_data = data_LF[[i]],
                                      est_sigma = "log_sigma_R",
                                      data_avail = "LC"
                                      # derive_quants = TRUE
                                     )
}

# calc_derived_quants(Medlived_res_LIME[[i]]$obj, lh[[i]])

for(i in 1:300){
  Medlived_res_LIME[[i]]$Derived <- Medlived_res_LIME2[[i]]$Derived
}

saveRDS(Medlived_res_LIME, file = "Medmodel_res_LIME.rds")

# devtools::install_github("merrillrudd/LIME")
# devtools::install_github("merrillrudd/LIME@b135ea8d2b9317c3e1a44e90f92a7cd9776938d1")

rm(list = ls())