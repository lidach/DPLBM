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
# max(list$inds[[1]][[1]][,5])
# list$inds[[1]][[1]][,6][max(list$inds[[1]][[1]][,5])]


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








## Thompson and Bell ################################

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


saveRDS(res_TB2, file = "Medmodel_res_TB.rds")
rm(list = ls())






## Length-Based Risk Analysis ####
## functions
calc_abund_SB <- function (ages, M, F, R0)
{
  N_a <- rep(NA, length(ages))
  N_a[1] <- R0
  for (i in 2:length(ages)) {
    if (i < length(ages))
      N_a[i] <- N_a[i - 1] * exp(-(M + F))
    if (i == length(ages))
      N_a[i] <- N_a[i - 1] * exp(-(M + F))/(1 - exp(-(M + F)))
  }
  return(N_a)
}

SPR_SB <- function(ages, Mat_a, W_a, M, F){
  Na0 <- calc_abund_SB(ages = ages, M = M, F = 0, R0 = R0)
  Naf <- calc_abund_SB(ages = ages, M = M, F = F, R0 = R0)
  # fished and unfished
  SB0 <- sum(Na0*Mat_a*W_a)
  SBf <- sum(Naf*Mat_a*W_a)
  # Compute and return SPR value
  SPR <- SBf/SB0
}

YPR_SB <- function(F, ages, M, R0, W_a) {
  Nage <- calc_abund_SB(ages = ages, M = M, F = F, R0 = R0)
  YPR <- F * sum(Nage * W_a)
  return(YPR)
}




modpath <- "D:/DPLBM/Pub1/Medlived/True"
setwd(modpath)
LFQmedmodel <- readRDS("LFQmedmodel.rds")
LFQmedmodel_etc <- readRDS("Medlist.rds")
iters <- 300


#' BH eq - fishmethods
Medlived_res_bheq <- list()
LFQmedmodel1 <- LFQmedmodel 

max <- list()
for (i in 1:iters){
  max[[i]] <- max(LFQmedmodel_etc$inds[[i]][[1]][,2])
  if (max[[i]] > LFQmedmodel1[[i]]$Linf){
    max[[i]] = LFQmedmodel1[[i]]$Linf
  }
}

res_cc <- readRDS("Medmodel_res_cc.rds")
for(i in 1:iters){
  LFQmedmodel1[[i]]$L50 <- res_cc[[i]]$L50
}


for(i in 1:iters){
  Medlived_res_bheq[[i]] <- bheq(LFQmedmodel_etc$inds[[i]][[1]][,2], type = 2, K = LFQmedmodel1[[i]]$K, Linf = LFQmedmodel1[[i]]$Linf, Lc = LFQmedmodel1[[i]]$L50, La = max[[i]], nboot = 500)
}

saveRDS(Medlived_res_bheq, file = "Medlived_res_bheq.rds")


#' Extract inputs
list <- list()
for(i in 1:iters){
  list[[i]] <- LFQmedmodel_etc$inds[[i]]
}


lwa <- 0.018
lwb <- 2.895
t0 <- -0.01
R0 <- 1

linf <- list()
vbk <- list()
ages <- list()
L_a <- list()
W_a <- list()
Mat_a <- list()
M <- list()
F <- list()
for (i in 1:iters){
  linf[[i]] <- LFQmedmodel1[[i]]$Linf
  vbk[[i]] <- LFQmedmodel1[[i]]$K
  ages[[i]] <- 0:LFQmedmodel1[[i]]$tmax
  L_a[[i]] <- linf[[i]]*(1-exp(-vbk[[i]]*(ages[[i]] - t0)))
  W_a[[i]] <- lwa*L_a[[i]]^lwb
  Mat_a[[i]] <- as.numeric(sum(list[[i]][[1]][,4] > LFQmedmodel1[[i]]$Lmat))
  M[[i]] <- LFQmedmodel1[[i]]$M
}

for(i in 1:iters){
  F[[i]] <- Medlived_res_bheq[[i]]$z - M[[i]]
  if(F[[i]] < 0){
    F[[i]] = 0
  }
}


#' run model
Medlived_res_SB <- list(SPR = rep(NA, 100), YPR = rep(NA, 100))
for (i in 1:iters){
  Medlived_res_SB$SPR[[i]] <- SPR_SB(ages[[i]], Mat_a[[i]], W_a[[i]], M[[i]], F[[i]])
  Medlived_res_SB$YPR[[i]] <- YPR_SB(F[[i]], ages[[i]], M[[i]], R0, W_a[[i]])
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




#' run model
lbspr_res <- list()
for(i in 1:iters){
  PARSmed <- new("LB_pars")
  PARSmed@Linf <- LFQmedmodel1[[i]]$Linf
  PARSmed@MK <- LFQmedmodel1[[i]]$M/LFQmedmodel1[[i]]$K
  PARSmed@L_units <- "cm"
  PARSmed@L50 <- LFQmedmodel[[i]]$Lmat
  PARSmed@L95 <- LFQmedmodel[[i]]$Lmat * 1.3 # definition
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


# par(mfrow=c(2,2))
# plot(lh$L_a, type="l", lwd=3, col="forestgreen", xlab="Age", ylab="Length")
# plot(lh$W_a, type="l", lwd=3, col="forestgreen", xlab="Age", ylab="Weight")
# plot(lh$Mat_l, type="l", lwd=3, col="forestgreen", xlab="Length", ylab="Proportion mature")
# plot(lh$S_l, type="l", lwd=3, col="forestgreen", xlab="Length", ylab="Proportion vulnerable to gear")


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

# a <- c(which(p %in% "The model is likely not converged"))
for (i in 1:300){
  Medlived_res_LIME[[i]] <- run_LIME(modpath = NULL,
                                      lh = lh[[i]],
                                      input_data = data_LF[[i]],
                                      est_sigma = "log_sigma_R",
                                      data_avail = "LC")
}

# calc_derived_quants(Medlived_res_LIME[[i]]$obj, lh[[i]])

saveRDS(Medlived_res_LIME, file = "Medmodel_res_LIME.rds")
