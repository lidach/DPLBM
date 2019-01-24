## Packages ################################
require(TropFishR)
require(LBSPR)
require(LIME)
require(fishdynr)
require(fishmethods)



## True Parameter Input ################################
# setwd("D:/DPLBM/Pub1/Errorlived")
# Errorlived <- readRDS("D:/DPLBM/Pub1/Errorlived/LFQerrormodel.rds")
Errorlived <- readRDS("LFQerrormodel.rds")
iters <- 300

Errorlived1 <- Errorlived

for(i in 1:iters){
  Errorlived1[[i]]$Linf <- 64.6
  Errorlived1[[i]]$K <- 0.21
  Errorlived1[[i]]$t0 <- -0.01
  Errorlived1[[i]]$M <- 0.32
  Errorlived1[[i]]$Lmat <- 34
  Errorlived1[[i]]$tmax <- 18
}


setwd("D:/DPLBM/Pub1/Errorlived/True")
saveRDS(Errorlived1, file = "LFQerrormodel.rds")

rm(list = ls())







## Thompson and Bell ################################
LFQerrormodel <- readRDS("LFQerrormodel.rds")
iters <- 300

#' LFQ conversions
LFQerrormodel1 <- LFQerrormodel
for (i in 1:iters){
  LFQerrormodel1[[i]] <- lfqModify(LFQerrormodel1[[i]], bin_size = 2)
  LFQerrormodel1[[i]] <- lfqModify(LFQerrormodel1[[i]], vectorise_catch = T)
}

#' CC analysis
LFQerrormodel2 <- LFQerrormodel1
res_cc <- list()
for (i in 1:iters){
  res_cc[[i]] <- catchCurve(LFQerrormodel2[[i]], auto=TRUE, calc_ogive = TRUE)
  res_cc[[i]]$wqs <- (res_cc[[i]]$L75 - res_cc[[i]]$L50)*2
}

LFQerrormodel3 <- LFQerrormodel2

for(i in 1:iters){
  LFQerrormodel3[[i]] <- lfqModify(LFQerrormodel3[[i]], plus_group = "Linf")
  LFQerrormodel3[[i]]$a <- 0.018
  LFQerrormodel3[[i]]$b <- 2.895
  LFQerrormodel3[[i]]$FM <- res_cc[[i]]$Z - LFQerrormodel3[[i]]$M  
  LFQerrormodel3[[i]]$E <- LFQerrormodel3[[i]]$FM / res_cc[[i]]$Z 
  LFQerrormodel3[[i]]$wmat <- 0.15 * LFQerrormodel3[[i]]$Lmat # default definition of wmat
}

## selectivity from CC
for(i in 1:iters){
  LFQerrormodel3[[i]]$FM <- LFQerrormodel3[[i]]$FM * logisticSelect(LFQerrormodel3[[i]]$midLengths,
                                                                  res_cc[[i]]$L50, res_cc[[i]]$wqs)
}



#' TB
res_TB2 <- list()
for (i in 1:iters){
  res_TB2[[i]] <- predict_mod(LFQerrormodel3[[i]],
                              type = "ThompBell",
                              FM_change = seq(0, 2.5, 0.01),
                              stock_size_1 = 1,
                              curr.E = LFQerrormodel3[[i]]$E,
                              plot = FALSE,
                              hide.progressbar = TRUE)
}

saveRDS(res_cc, file = "Errormodel_res_cc.rds")
saveRDS(res_TB2, file = "Errormodel_res_TB.rds")
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


modpath <- "D:/DPLBM/Pub1/Errorlived/True"
setwd(modpath)
LFQerrormodel <- readRDS("LFQerrormodel.rds")
iters <- 300


#' BH eq - fishmethods
Errorlived_res_bheq <- list()
LFQerrormodel1 <- LFQerrormodel 

len <- list()
for(i in 1:iters){ 
  len[[i]] <- rep(LFQerrormodel1[[i]]$midLengths, rowSums(LFQerrormodel1[[i]]$catch))
}

max <- list()
for (i in 1:iters){
  max[[i]] <- max(LFQerrormodel1[[i]]$midLengths)
  if (max[[i]] > LFQerrormodel1[[i]]$Linf){
    max[[i]] = LFQerrormodel1[[i]]$Linf
  }
}

res_cc <- readRDS("Errormodel_res_cc.rds")
for(i in 1:iters){
  LFQerrormodel1[[i]]$L50 <- res_cc[[i]]$L50
}

set.seed(1)

source("~/DPLBM/bheq_mod.R")

for(i in 1:iters){
  Errorlived_res_bheq[[i]] <- bheq_LBRA(len[[i]], type = 2, K = LFQerrormodel1[[i]]$K, Linf = LFQerrormodel1[[i]]$Linf, Lc = LFQerrormodel1[[i]]$L50, La = max[[i]], nboot = 200)
}

saveRDS(Errorlived_res_bheq, file = "Errormodel_res_bheq.rds")


## beta --> find M
## 1/beta = -log(0.001)/delta_a_lambda
## delta_a_lambda = theor_a - obs_a
del_a <- NA
dev <- NA
for(i in 1:iters){
  del_a[i] <-  LFQerrormodel1[[i]]$tmax - ceiling(-log(0.001)/LFQerrormodel1[[i]]$M)
  dev[i] <- ceiling(-log(0.001)/LFQerrormodel1[[i]]$M) + (del_a[i] / -log(0.001))
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
  linf[[i]] <- LFQerrormodel1[[i]]$Linf
  vbk[[i]] <- LFQerrormodel1[[i]]$K
  ages[[i]] <- seq(0, LFQerrormodel1[[i]]$tmax, tincr)
  L_a[[i]] <- linf[[i]]*(1-exp(-vbk[[i]]*(ages[[i]] - t0)))
  
  mids[[i]] <- seq((binwidth/2), linf[[i]]*1.5, binwidth)
  highs[[i]] <- mids[[i]] + (binwidth/2)
  lows[[i]] <- mids[[i]] - (binwidth/2)
  plba_a[[i]] <- t(vlprobs(L_a[[i]], L_a[[i]] * CV))
  plba_a[[i]] <- plba_a[[i]] / rowSums(plba_a[[i]])
  
  W_a[[i]] <- lwa*L_a[[i]]^lwb
  
  M[[i]] <- -log(0.001)/dev[i]
  # M[[i]] <- LFQerrormodel1[[i]]$M
  
  # instead of knife-edge --> logistic
  Mat_a[[i]] <- 1 / (1 + exp(-(mids[[i]] - LFQerrormodel1[[i]]$Lmat) /
                               ((LFQerrormodel1[[i]]$Lmat * 0.15) / ( log(0.75/(1-0.75)) - log(0.25/(1-0.25)) ))))
  Mat_a[[i]] <- apply(t(plba_a[[i]])*Mat_a[[i]], 2, sum)
  Sl_a[[i]] <- 1 / (1 + exp(-(mids[[i]] - LFQerrormodel1[[i]]$L50) /
                              (res_cc[[i]]$wqs / ( log(0.75/(1-0.75)) - log(0.25/(1-0.25)) ))))
  Sl_a[[i]] <- apply(t(plba_a[[i]])*Sl_a[[i]], 2, sum)
}

for(i in 1:iters){
  FM[[i]] <- Errorlived_res_bheq[[i]]$z - M[[i]]
  if(FM[[i]] < 0){
    FM[[i]] = 0
  }
}


#' run model
Errorlived_res_SB <- list(SPR = NA, YPR = NA, Fmsy = NA, Bmsy = NA, FFmsy = NA, BBmsy = NA, FM = NA)
for (i in 1:iters){
  Errorlived_res_SB$SPR[i] <- SPR_SB(ages[[i]], Sl_a[[i]], Mat_a[[i]], W_a[[i]], M[[i]], FM[[i]])
  Errorlived_res_SB$YPR[i] <- YPR_SB(FM[[i]], ages[[i]], M[[i]], R0, W_a[[i]], Sl_a[[i]])
  Errorlived_res_SB$Fmsy[i] <- optimize(YPR_SB, ages = ages[[i]], M = M[[i]], R0 = R0, W_a = W_a[[i]], Sl_a = Sl_a[[i]], lower = 0, upper = 10, maximum = TRUE)$maximum
  Errorlived_res_SB$Bmsy[i] <- SSB_SB(ages[[i]], Sl_a[[i]], Mat_a[[i]], W_a[[i]], M[[i]], Errorlived_res_SB$Fmsy[i])
  Errorlived_res_SB$FFmsy[i] <- FM[[i]]/Errorlived_res_SB$Fmsy[i]
  Errorlived_res_SB$BBmsy[i] <- SSB_SB(ages[[i]], Sl_a[[i]], Mat_a[[i]], W_a[[i]], M[[i]], FM[[i]]) / Errorlived_res_SB$Bmsy[i]
  Errorlived_res_SB$FM[i] <- FM[[i]]
}


#' save SB
saveRDS(Errorlived_res_SB, file = "Errormodel_res_LBRA.rds")

rm(list = ls())





## LBSPR ####
modpath <- "D:/DPLBM/Pub1/Errorlived/True"
setwd(modpath)
LFQerrormodel <- readRDS("LFQerrormodel.rds")
iters <- 300

LFQerrormodel1 <- LFQerrormodel
a <- 0.018
b <- 2.895


for(i in 1:iters){
  LFQerrormodel1[[i]] <- lfqModify(LFQerrormodel1[[i]], bin_size = 2)
}


pg <- list()
for(i in 1:iters){
  if(max(LFQerrormodel1[[i]]$midLengths) < LFQerrormodel1[[i]]$Linf){
    closest<-function(xv,sv){
      xv[which(abs(xv-sv)==min(abs(xv-sv)))] 
    }
    
    k <- seq(min(LFQerrormodel1[[i]]$midLengths),120,2)
    l <- closest(k, LFQerrormodel1[[i]]$Linf+2)
    
    pg[[i]] <- c(LFQerrormodel1[[i]]$midLengths[-length(LFQerrormodel1[[i]]$midLengths)], seq(max(LFQerrormodel1[[i]]$midLengths),l,2))
  }else{
    pg[[i]] <- LFQerrormodel1[[i]]$midLengths
  }
}

catch <- list()
for(i in 1:iters){
  LFQerrormodel1[[i]]$catch <- rowSums(LFQerrormodel1[[i]]$catch)
  catch[[i]] <- c(LFQerrormodel1[[i]]$catch,rep(0,length(LFQerrormodel1[[i]]$midLengths[pg[[i]] > max(LFQerrormodel1[[i]]$midLengths)])))
}


res_cc <- readRDS("Errormodel_res_cc.rds")
for(i in 1:iters){
  LFQerrormodel1[[i]]$SL50 <- res_cc[[i]]$L50
  LFQerrormodel1[[i]]$SL95 <- res_cc[[i]]$L95
}




#' run model
lbspr_res <- list()
for(i in 1:iters){
  PARSerror <- new("LB_pars")
  PARSerror@Linf <- LFQerrormodel1[[i]]$Linf
  PARSerror@MK <- LFQerrormodel1[[i]]$M/LFQerrormodel1[[i]]$K
  PARSerror@L_units <- "cm"
  PARSerror@L50 <- LFQerrormodel[[i]]$Lmat
  PARSerror@L95 <- LFQerrormodel[[i]]$Lmat * 1.3 # definition
  PARSerror@Walpha <- a
  PARSerror@Wbeta <- b
  PARSerror@BinWidth <- 2
  
  SPRerror <- new("LB_lengths", LB_pars = PARSerror)
  SPRerror@LMids <- pg[[i]]
  SPRerror@LData <- as.matrix(catch[[i]])
  SPRerror@L_units <- "cm"
  SPRerror@Years <- 1
  SPRerror@NYears <- 1
  
  lbspr_res[[i]] <- tryCatch(LBSPRfit(LB_pars = PARSerror, LB_lengths = SPRerror, yrs = 1, Control = list(modtype = "GTG")))
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
}
saveRDS(LBSPR_outs, file = "Errormodel_res_LBSPR.rds")

rm(list = ls())








## LIME ####
modpath <- "D:/DPLBM/Pub1/Errorlived/True"
setwd(modpath)
LFQerrormodel <- readRDS("LFQerrormodel.rds")
iters <- 300

#' add SL50 and 95
LFQerrormodel1 <- LFQerrormodel
a = 0.018
b = 2.895

res_cc <- readRDS("Errormodel_res_cc.rds")
for(i in 1:iters){
  LFQerrormodel1[[i]]$SL50 <- res_cc[[i]]$L50
  LFQerrormodel1[[i]]$SL95 <- res_cc[[i]]$L95
}

#' Input parameters for lh list
lh <- list()
for(i in 1:iters){
  lh[[i]] <- create_lh_list(vbk = LFQerrormodel1[[i]]$K, 
                            linf = LFQerrormodel1[[i]]$Linf,
                            t0 = -0.01,
                            lwa = a,
                            lwb = b,
                            S50 = LFQerrormodel1[[i]]$SL50,
                            S95 = LFQerrormodel1[[i]]$SL95,
                            selex_input = "length",
                            selex_type = "logistic",
                            M50 = LFQerrormodel1[[i]]$Lmat,
                            M95 = NULL,
                            maturity_input = "length",
                            M = LFQerrormodel1[[i]]$M,
                            binwidth = 2,
                            R0 = 1,
                            nseasons = 1)
}



#' set up data
for(i in 1:iters){
  LFQerrormodel1[[i]] <- lfqModify(LFQerrormodel1[[i]], bin_size = 2)
  LFQerrormodel1[[i]]$catch <- rowSums(LFQerrormodel1[[i]]$catch)
  LFQerrormodel1[[i]]$catch <- as.matrix(LFQerrormodel1[[i]]$catch)
}

#' LF list
LF <- list()
for(i in 1:iters){
  LF[[i]] <- t(LFQerrormodel1[[i]]$catch)
  rownames(LF[[i]]) <- as.numeric(1)
  colnames(LF[[i]]) <- LFQerrormodel1[[i]]$midLengths
}

#' years with length data
data_LF <- list()
for (i in 1:iters){
  data_LF[[i]] <- list("years" = 1, "LF" = LF[[i]])
}

#' Run LIME
Errorlived_res_LIME <- list()
# data_all <- list()
# 
# for(i in 1:iters){
#   data_all[[i]] <- create_inputs(lh = lh[[i]], input_data = data_LF[[i]])
# }

# a <- c(which(p %in% "The model is likely not converged"))
for (i in 1:300){
  Errorlived_res_LIME[[i]] <- run_LIME(modpath = NULL,
                                       # input = data_all[[i]],
                                       lh = lh[[i]],
                                       input_data = data_LF[[i]],
                                       est_sigma = "log_sigma_R",
                                       data_avail = "LC"
                                       # derive_quants = TRUE
  )
}

# calc_derived_quants(Errorlived_res_LIME[[i]]$obj, lh[[i]])

for(i in 1:300){
  Errorlived_res_LIME[[i]]$Derived <- Errorlived_res_LIME2[[i]]$Derived
}

saveRDS(Errorlived_res_LIME, file = "Errormodel_res_LIME.rds")

# devtools::install_github("merrillrudd/LIME")
# devtools::install_github("merrillrudd/LIME@b135ea8d2b9317c3e1a44e90f92a7cd9776938d1")

rm(list = ls())
