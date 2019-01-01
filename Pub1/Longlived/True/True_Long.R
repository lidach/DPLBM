## Packages ################################
require(TropFishR)
require(LBSPR)
require(LIME)
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
  res_cc[[i]] <- catchCurve(LFQlongmodel2[[i]], reg_int = c(7,24), calc_ogive = T)
}


# res_cc[[i]] <- catchCurve(LFQlongmodel2[[i]], reg_int = NULL, calc_ogive = T)

p <- NA
for (i in 1:iters){
  if(res_cc[[i]]$FM[[1]] < 0) {
    p[i] <- F }
  else{
    p[i] <- 1}
}

which(p %in% 0)

for(i in 1:iters){
  p[[i]] <- res_cc[[i]]$FM
}


for(i in 1:iters){
  LFQlongmodel2[[i]]$Z <- res_cc[[i]]$Z
  LFQlongmodel2[[i]]$FM <- as.numeric(res_cc[[i]]$Z - res_cc[[i]]$M)
  LFQlongmodel2[[i]]$E <- res_cc[[i]]$FM[[1]]/res_cc[[i]]$Z
  LFQlongmodel2[[i]]$L50 <- res_cc[[i]]$L50
  LFQlongmodel2[[i]]$L75 <- res_cc[[i]]$L75
}

saveRDS(res_cc, file = "Longmodel_res_cc.rds")


#' CA
LFQlongmodel3 <- LFQlongmodel2
# a <- c(which(p %in% NA))
for(i in 1:iters){
  LFQlongmodel3[[i]] <- lfqModifydev(LFQlongmodel3[[i]], plus_group = TRUE)
  LFQlongmodel3[[i]]$a <- 0.0123
  LFQlongmodel3[[i]]$b <- 3.035
}


dev.off()
vpa_res <- list()
for(i in 1:iters){
  vpa_res[[i]] <- VPA(LFQlongmodel3[[i]], terminalF = LFQlongmodel3[[i]]$FM, analysis_type = "CA", plot = F)
}

for(i in 1:iters){
  LFQlongmodel3[[i]]$FM <- vpa_res[[i]]$FM_calc
}


p <- NA
for (i in 1:iters){
  if(vpa_res[[i]]$annualMeanNr[1] < 0) {
    p[i] <- F }
  else{
    p[i] <- 1}
}



saveRDS(vpa_res, file = "Longmodel_res_CA.rds")


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



#' save data
saveRDS(res_TB2, file = "Longmodel_res_TB2.rds")

rm(list = ls())








## LBRP ################################
## functions 
HCR.VBGF <- function(Linf, K, t0, ages){
  Lengths_exp <- Linf * (1 - exp(-K * (ages - t0)))
  return(Lengths_exp)
}

HCR.LtWt.fit <- function(a, b, Lts){
  exp.wts <- a*Lts^b
  return(exp.wts)
}

Calc_Pobjs <- function(catch_prop, midLengths, Lmat, Lopt0.9, Lopt1.1){
  Lmat_ind <- midLengths >= Lmat
  Lopt_ind <- midLengths > Lopt0.9 & midLengths < Lopt1.1
  Lmega_ind <- midLengths >= Lopt1.1
  Pmat <- sum(catch_prop[Lmat_ind])
  Popt <- sum(catch_prop[Lopt_ind])
  Pmega <- sum(catch_prop[Lmega_ind])
  Pobj <- Pmat + Popt + Pmega
  Pmat.opt.mega <- matrix(c(Pmat, Popt, Pmega, Pobj), nrow = 4, ncol = 1)
  rownames(Pmat.opt.mega) <- c("Pmat", "Popt", "Pmega", "Pobj")
  colnames(Pmat.opt.mega) <- "Value"
  return(Pmat.opt.mega)
}

Froese.plus <- function(catch_prop, midLengths, Linf, K, t0, M, a, b, Lmat, ages){
  Lts <- HCR.VBGF(Linf, K, t0, ages)
  Wts <- HCR.LtWt.fit(a, b, Lts)
  stable.biomass <- (1*exp(-M*ages))*Wts
  Lopt <- Lts[stable.biomass == max(stable.biomass)]
  Lopt0.9 <- Lopt * 0.9
  Lopt1.1 <- round(1.1 * Lopt, 1)
  yrs <- 1
  Pobjs.out <- matrix(NA, nrow = 4, ncol = length(yrs), dimnames = list(c("Pmat", "Popt", "Pmega", "Pobj"), yrs))
  for(i in 1:length(yrs)){
    cab.Pobjs.yr <- Calc_Pobjs(catch_prop, midLengths, Lmat, Lopt0.9, Lopt1.1)
    Pobjs.out[,i] <- as.numeric(cab.Pobjs.yr)
  }
  Pcalc.out <- list()
  Pcalc.out[[1]] <- Pobjs.out
  Pcalc.out[[2]] <- c(Lmat, Lopt, Lopt1.1)
  names(Pcalc.out[[2]]) <- c("Lmat", "Lopt", "Lmega")
  Pcalc.out[[3]]<- Lmat/Lopt
  names(Pcalc.out)[[1]] <- "Pout"
  names(Pcalc.out)[[2]] <- "Lx"
  names(Pcalc.out)[[3]] <- "Lmat/Lopt"
  return(Pcalc.out)
}



modpath <- "D:/DPLBM/Pub1/Longlived/True"
setwd(modpath)
LFQlongmodel <- readRDS("LFQlongmodel.rds")
iters <- 300


#' set up catch prop
catch <- list()
catch_year <- list()
LFQlongmodel1 <- LFQlongmodel
for(i in 1:iters){
  catch[[i]] <- t(LFQlongmodel1[[i]]$catch)
  catch_year[[i]] <- rep(NA, length.out = ncol(catch[[i]]))
  catch_year[[i]] <- colSums(catch[[i]])
}

catch_prop <- list()
for(i in 1:length(catch_year)){
  catch_prop[[i]] <- catch_year[[i]]/sum(catch_year[[i]])
}

ages <- list()
for(i in 1:iters){
  ages[[i]] <- 0:LFQlongmodel1[[i]]$tmax
}

t0 <- -0.01
a = 0.0123
b = 3.035

Linf <- list()
K <- list()
M <- list()
Lmat <- list()
midLengths <- list()
for(i in 1:iters){
  Linf[[i]] <- LFQlongmodel1[[i]]$Linf 
  K[[i]] <- LFQlongmodel1[[i]]$K
  M[[i]] <- LFQlongmodel1[[i]]$M
  Lmat[[i]] <- LFQlongmodel1[[i]]$Lmat
  midLengths[[i]] <- LFQlongmodel1[[i]]$midLengths
}

Longmodel_res_LBRP <- list()
for(i in 1:iters){
  Longmodel_res_LBRP[[i]] <- Froese.plus(catch_prop[[i]], midLengths[[i]], Linf[[i]], K[[i]], t0, M[[i]], a, b, Lmat[[i]], ages[[i]])
}

#' save data
saveRDS(Longmodel_res_LBRP, file = "Longmodel_res_LBRP.rds")

rm(list = ls())








## Sustainability Benchmarks ####
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




modpath <- "D:/DPLBM/Pub1/Longlived/True"
setwd(modpath)
LFQlongmodel <- readRDS("LFQlongmodel.rds")
LFQlongmodel_etc <- readRDS("Longmodel_list.rds")
iters <- 300


#' BH eq - fishmethods
Longmodel_res_bheq <- list()
LFQlongmodel1 <- LFQlongmodel 

max <- list()
for (i in 1:iters){
  max[[i]] <- max(LFQlongmodel_etc$inds[[i]][[1]][,2])
  if (max[[i]] > LFQlongmodel1[[i]]$Linf){
    max[[i]] = LFQlongmodel1[[i]]$Linf
  }
}

res_cc <- readRDS("Longmodel_res_cc.rds")
for(i in 1:iters){
  LFQlongmodel1[[i]]$L50 <- res_cc[[i]]$L50
}


for(i in 1:iters){
  Longmodel_res_bheq[[i]] <- bheq(LFQlongmodel_etc$inds[[i]][[1]][,2], type = 2, K = LFQlongmodel1[[i]]$K, Linf = LFQlongmodel1[[i]]$Linf, Lc = LFQlongmodel1[[i]]$L50, La = max[[i]], nboot = 500)
}

saveRDS(Longmodel_res_bheq, file = "Longmodel_res_bheq.rds")


#' Extract inputs
list <- list()
for(i in 1:iters){
  list[[i]] <- LFQlongmodel_etc$inds[[i]]
}


lwa = 0.0123
lwb = 3.035
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
  linf[[i]] <- LFQlongmodel1[[i]]$Linf
  vbk[[i]] <- LFQlongmodel1[[i]]$K
  ages[[i]] <- 0:LFQlongmodel1[[i]]$tmax
  L_a[[i]] <- linf[[i]]*(1-exp(-vbk[[i]]*(ages[[i]] - t0)))
  W_a[[i]] <- lwa*L_a[[i]]^lwb
  Mat_a[[i]] <- as.numeric(sum(list[[i]][[1]][,4] > LFQlongmodel1[[i]]$Lmat))
  M[[i]] <- LFQlongmodel1[[i]]$M
}

for(i in 1:iters){
  F[[i]] <- Longmodel_res_bheq[[i]]$z - M[[i]]
  if(F[[i]] < 0){
    F[[i]] = 0
  }
}


#' run model
Longmodel_res_SB <- list(SPR = rep(NA, 100), YPR = rep(NA, 100))
for (i in 1:iters){
  Longmodel_res_SB$SPR[[i]] <- SPR_SB(ages[[i]], Mat_a[[i]], W_a[[i]], M[[i]], F[[i]])
  Longmodel_res_SB$YPR[[i]] <- YPR_SB(F[[i]], ages[[i]], M[[i]], R0, W_a[[i]])
}

#' save SB
saveRDS(Longmodel_res_SB, file = "Longmodel_res_SB.rds")

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
}

#' save data
LBSPR_outs <- list(pLCatch = list(rep(NA, 100)), SL50 = NA, SL95 = NA, FM = NA, SPR = NA, SPR_Var = NA, SL50_Var = NA, SL95_Var = NA, FM_Var = NA)
for (i in 1:iters){
  LBSPR_outs$pLCatch[[i]] <- lbspr_res[[i]]@pLCatch
  LBSPR_outs$SL50[[i]] <- lbspr_res[[i]]@Ests[,"SL50"]
  LBSPR_outs$SL95[[i]] <- lbspr_res[[i]]@Ests[,"SL95"]
  LBSPR_outs$FM[[i]] <- lbspr_res[[i]]@Ests[,"FM"]
  LBSPR_outs$SPR[[i]] <- lbspr_res[[i]]@Ests[,"SPR"]
  LBSPR_outs$SPR_Var[[i]] <- lbspr_res[[i]]@Vars[,"SPR"]
  LBSPR_outs$SL50_Var[[i]] <- lbspr_res[[i]]@Vars[,"SL50"]
  LBSPR_outs$SL95_Var[[i]] <- lbspr_res[[i]]@Vars[,"SL95"]
  LBSPR_outs$FM_Var[[i]] <- lbspr_res[[i]]@Vars[,"FM"]
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
                            nseasons = 1)
}


# par(mfrow=c(2,2))
# plot(lh$L_a, type="l", lwd=3, col="forestgreen", xlab="Age", ylab="Length")
# plot(lh$W_a, type="l", lwd=3, col="forestgreen", xlab="Age", ylab="Weight")
# plot(lh$Mat_l, type="l", lwd=3, col="forestgreen", xlab="Length", ylab="Proportion mature")
# plot(lh$S_l, type="l", lwd=3, col="forestgreen", xlab="Length", ylab="Proportion vulnerable to gear")


#' set up data
for(i in 1:iters){
  LFQlongmodel1[[i]] <- lfqModify(LFQlongmodel1[[i]], bin_size = 3)
  # LFQlongmodel1[[i]]$catch <- rowSums(LFQlongmodel1[[i]]$catch)
  LFQlongmodel1[[i]]$catch <- as.matrix(LFQlongmodel1[[i]]$catch)
}

#' LF list
LF <- list()
for(i in 1:iters){
  LF[[i]] <- t(LFQlongmodel1[[i]]$catch)
  rownames(LF[[i]]) <- as.numeric(1:12)
  colnames(LF[[i]]) <- LFQlongmodel1[[i]]$midLengths
}

#' years with length data
data_LF <- list()
for (i in 1:iters){
  data_LF[[i]] <- list("years" = 1:12, "LF" = LF[[i]])
}

#' Run LIME
Longmodel_res_LIME <- list()

# a <- c(which(p %in% "The model is likely not converged"))
for (i in 1:iters){
  Longmodel_res_LIME[[i]] <- run_LIME(modpath = NULL,
                                      lh = lh[[i]],
                                      input_data = data_LF[[i]],
                                      est_sigma = "log_sigma_R",
                                      data_avail = "LC")
}

p <- list()
for (i in 1:iters){
  p[[i]] <- mean(Longmodel_res_LIME[[i]]$Report$SPR_t)
}


saveRDS(Longmodel_res_LIME, file = "Longmodel_res_LIME.rds")
