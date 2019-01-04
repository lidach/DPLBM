## Packages ################################
require(TropFishR)
require(LBSPR)
require(LIME)
require(fishmethods)






## True Parameter Input ################################
# setwd("D:/DPLBM/Pub1/Equilibrium")
# Equili <- readRDS("D:/DPLBM/Pub1/Equilibrium/LFQequilimodel.rds")
Equili <- readRDS("LFQequilimodel.rds")
list <- readRDS("Equilimodel_list.rds")
iters <- 300

Equili1 <- Equili
# max(list$inds[[1]][[1]][,5])
# list$inds[[1]][[1]][,6][max(list$inds[[1]][[1]][,5])]


for(i in 1:iters){
  Equili1[[i]]$Linf <- 64.6
  Equili1[[i]]$K <- 0.21
  Equili1[[i]]$t0 <- -0.01
  Equili1[[i]]$M <- 0.32
  Equili1[[i]]$Lmat <- 34
  Equili1[[i]]$tmax <- 18
}


setwd("D:/DPLBM/Pub1/Equilibrium/True")
saveRDS(Equili1, file = "LFQequilimodel.rds")








## Thompson and Bell ################################
modpath = "D:/DPLBM/Pub1/Equilibrium/True"
setwd(modpath)
LFQequilimodel <- readRDS("LFQequilimodel.rds")
iters <- 300


#' LFQ conversions
LFQequilimodel1 <- LFQequilimodel
for (i in 1:iters){
  LFQequilimodel1[[i]] <- lfqModify(LFQequilimodel1[[i]], bin_size = 2)
  LFQequilimodel1[[i]] <- lfqModify(LFQequilimodel1[[i]], vectorise_catch = T)
}

#' CC analysis
LFQequilimodel2 <- LFQequilimodel1
res_cc <- list()
for (i in 1:iters){
  res_cc[[i]] <- catchCurve(LFQequilimodel2[[i]], reg_int = c(7,28), calc_ogive = T)
}


# res_cc[[i]] <- catchCurve(LFQequilimodel2[[i]], reg_int = NULL, calc_ogive = T)

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
  LFQequilimodel2[[i]]$Z <- res_cc[[i]]$Z
  LFQequilimodel2[[i]]$FM <- as.numeric(res_cc[[i]]$Z - res_cc[[i]]$M)
  LFQequilimodel2[[i]]$E <- res_cc[[i]]$FM[[1]]/res_cc[[i]]$Z
  LFQequilimodel2[[i]]$L50 <- res_cc[[i]]$L50
  LFQequilimodel2[[i]]$L75 <- res_cc[[i]]$L75
}

saveRDS(res_cc, file = "Equili_res_cc.rds")


#' CA
LFQequilimodel3 <- LFQequilimodel2
# a <- c(which(p %in% NA))
source("D:/DPLBM/Pub1/lfqModifydev.R")
for(i in 1:iters){
  LFQequilimodel3[[i]] <- lfqModifydev(LFQequilimodel3[[i]], plus_group = TRUE)
  LFQequilimodel3[[i]]$a <- 0.018
  LFQequilimodel3[[i]]$b <- 2.895
}


dev.off()
vpa_res <- list()
for(i in 1:iters){
  vpa_res[[i]] <- VPA(LFQequilimodel3[[i]], terminalF = LFQequilimodel3[[i]]$FM, analysis_type = "CA", plot = F)
}

for(i in 1:iters){
  LFQequilimodel3[[i]]$FM <- vpa_res[[i]]$FM_calc
}


p <- NA
for (i in 1:iters){
  if(vpa_res[[i]]$annualMeanNr[1] < 0) {
    p[i] <- F }
  else{
    p[i] <- 1}
}



saveRDS(vpa_res, file = "Equili_res_CA.rds")


#' TB
res_TB2 <- list()
for (i in 1:iters){
  res_TB2[[i]] <- predict_mod(LFQequilimodel3[[i]],
                              type = "ThompBell",
                              FM_change = seq(0, 2.5, 0.01),
                              stock_size_1 = 1,
                              curr.E = LFQequilimodel3[[i]]$E,
                              plot = FALSE,
                              hide.progressbar = TRUE)
}



#' save data
saveRDS(res_TB2, file = "Equili_res_TB2.rds")

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



modpath <- "D:/DPLBM/Pub1/Equilibrium/True"
setwd(modpath)
LFQequilimodel <- readRDS("LFQequilimodel.rds")
iters <- 300


#' set up catch prop
catch <- list()
catch_year <- list()
LFQequilimodel1 <- LFQequilimodel
for(i in 1:iters){
  catch[[i]] <- t(LFQequilimodel1[[i]]$catch)
  catch_year[[i]] <- rep(NA, length.out = ncol(catch[[i]]))
  catch_year[[i]] <- colSums(catch[[i]])
}

catch_prop <- list()
for(i in 1:length(catch_year)){
  catch_prop[[i]] <- catch_year[[i]]/sum(catch_year[[i]])
}

ages <- list()
for(i in 1:iters){
  ages[[i]] <- 0:LFQequilimodel1[[i]]$tmax
}

t0 <- -0.01
a <- 0.018
b <- 2.895

Linf <- list()
K <- list()
M <- list()
Lmat <- list()
midLengths <- list()
for(i in 1:iters){
  Linf[[i]] <- LFQequilimodel1[[i]]$Linf 
  K[[i]] <- LFQequilimodel1[[i]]$K
  M[[i]] <- LFQequilimodel1[[i]]$M
  Lmat[[i]] <- LFQequilimodel1[[i]]$Lmat
  midLengths[[i]] <- LFQequilimodel1[[i]]$midLengths
}

Equili_res_LBRP <- list()
for(i in 1:iters){
  Equili_res_LBRP[[i]] <- Froese.plus(catch_prop[[i]], midLengths[[i]], Linf[[i]], K[[i]], t0, M[[i]], a, b, Lmat[[i]], ages[[i]])
}

#' save data
saveRDS(Equili_res_LBRP, file = "Equilimodel_res_LBRP.rds")

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




modpath <- "D:/DPLBM/Pub1/Equilibrium/True"
setwd(modpath)
LFQequilimodel <- readRDS("LFQequilimodel.rds")
LFQequilimodel_etc <- readRDS("Equilimodel_list.rds")
iters <- 300


#' BH eq - fishmethods
Equili_res_bheq <- list()
LFQequilimodel1 <- LFQequilimodel 

max <- list()
for (i in 1:iters){
  max[[i]] <- max(LFQequilimodel_etc$inds[[i]][[1]][,2])
  if (max[[i]] > LFQequilimodel1[[i]]$Linf){
    max[[i]] = LFQequilimodel1[[i]]$Linf
  }
}

res_cc <- readRDS("Equili_res_cc.rds")
for(i in 1:iters){
  LFQequilimodel1[[i]]$L50 <- res_cc[[i]]$L50
}


for(i in 1:iters){
  Equili_res_bheq[[i]] <- bheq(LFQequilimodel_etc$inds[[i]][[1]][,2], type = 2, K = LFQequilimodel1[[i]]$K, Linf = LFQequilimodel1[[i]]$Linf, Lc = LFQequilimodel1[[i]]$L50, La = max[[i]], nboot = 500)
}

saveRDS(Equili_res_bheq, file = "Equili_res_bheq.rds")


#' Extract inputs
list <- list()
for(i in 1:iters){
  list[[i]] <- LFQequilimodel_etc$inds[[i]]
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
  linf[[i]] <- LFQequilimodel1[[i]]$Linf
  vbk[[i]] <- LFQequilimodel1[[i]]$K
  ages[[i]] <- 0:LFQequilimodel1[[i]]$tmax
  L_a[[i]] <- linf[[i]]*(1-exp(-vbk[[i]]*(ages[[i]] - t0)))
  W_a[[i]] <- lwa*L_a[[i]]^lwb
  Mat_a[[i]] <- as.numeric(sum(list[[i]][[1]][,4] > LFQequilimodel1[[i]]$Lmat))
  M[[i]] <- LFQequilimodel1[[i]]$M
}

for(i in 1:iters){
  F[[i]] <- Equili_res_bheq[[i]]$z - M[[i]]
  if(F[[i]] < 0){
    F[[i]] = 0
  }
}


#' run model
Equili_res_SB <- list(SPR = rep(NA, 100), YPR = rep(NA, 100))
for (i in 1:iters){
  Equili_res_SB$SPR[[i]] <- SPR_SB(ages[[i]], Mat_a[[i]], W_a[[i]], M[[i]], F[[i]])
  Equili_res_SB$YPR[[i]] <- YPR_SB(F[[i]], ages[[i]], M[[i]], R0, W_a[[i]])
}

#' save SB
saveRDS(Equili_res_SB, file = "Equili_res_SB.rds")

rm(list = ls())





## LBSPR ####
modpath <- "D:/DPLBM/Pub1/Equilibrium/True"
setwd(modpath)
LFQequilimodel <- readRDS("LFQequilimodel.rds")
iters <- 300

LFQequilimodel1 <- LFQequilimodel
a <- 0.018
b <- 2.895


for(i in 1:iters){
  LFQequilimodel1[[i]] <- lfqModify(LFQequilimodel1[[i]], bin_size = 2)
}


pg <- list()
for(i in 1:iters){
  if(max(LFQequilimodel1[[i]]$midLengths) < LFQequilimodel1[[i]]$Linf){
    closest<-function(xv,sv){
      xv[which(abs(xv-sv)==min(abs(xv-sv)))] 
    }
    
    k <- seq(min(LFQequilimodel1[[i]]$midLengths),120,2)
    l <- closest(k, LFQequilimodel1[[i]]$Linf+2)
    
    pg[[i]] <- c(LFQequilimodel1[[i]]$midLengths[-length(LFQequilimodel1[[i]]$midLengths)], seq(max(LFQequilimodel1[[i]]$midLengths),l,2))
  }else{
    pg[[i]] <- LFQequilimodel1[[i]]$midLengths
  }
}

catch <- list()
for(i in 1:iters){
  LFQequilimodel1[[i]]$catch <- rowSums(LFQequilimodel1[[i]]$catch)
  catch[[i]] <- c(LFQequilimodel1[[i]]$catch,rep(0,length(LFQequilimodel1[[i]]$midLengths[pg[[i]] > max(LFQequilimodel1[[i]]$midLengths)])))
}


res_cc <- readRDS("Equili_res_cc.rds")
for(i in 1:iters){
  LFQequilimodel1[[i]]$SL50 <- res_cc[[i]]$L50
  LFQequilimodel1[[i]]$SL95 <- res_cc[[i]]$L95
}




#' run model
lbspr_res <- list()
for(i in 1:iters){
  PARSequili <- new("LB_pars")
  PARSequili@Linf <- LFQequilimodel1[[i]]$Linf
  PARSequili@MK <- LFQequilimodel1[[i]]$M/LFQequilimodel1[[i]]$K
  PARSequili@L_units <- "cm"
  PARSequili@L50 <- LFQequilimodel[[i]]$Lmat
  PARSequili@L95 <- LFQequilimodel[[i]]$Lmat * 1.3 # definition
  PARSequili@Walpha <- a
  PARSequili@Wbeta <- b
  PARSequili@BinWidth <- 2
  
  SPRequili <- new("LB_lengths", LB_pars = PARSequili)
  SPRequili@LMids <- pg[[i]]
  SPRequili@LData <- as.matrix(catch[[i]])
  SPRequili@L_units <- "cm"
  SPRequili@Years <- 1
  SPRequili@NYears <- 1
  
  lbspr_res[[i]] <- tryCatch(LBSPRfit(LB_pars = PARSequili, LB_lengths = SPRequili, yrs = 1, Control = list(modtype = "GTG")))
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
saveRDS(LBSPR_outs, file = "Equili_res_LBSPR.rds")

rm(list = ls())







## LIME ####
modpath <- "D:/DPLBM/Pub1/Equilibrium/True"
setwd(modpath)
LFQequilimodel <- readRDS("LFQequilimodel.rds")
iters <- 300

#' add SL50 and 95
LFQequilimodel1 <- LFQequilimodel
a = 0.018
b = 2.895

res_cc <- readRDS("Equili_res_cc.rds")
for(i in 1:iters){
  LFQequilimodel1[[i]]$SL50 <- res_cc[[i]]$L50
  LFQequilimodel1[[i]]$SL95 <- res_cc[[i]]$L95
}

#' Input parameters for lh list
lh <- list()
for(i in 1:iters){
  lh[[i]] <- create_lh_list(vbk = LFQequilimodel1[[i]]$K, 
                            linf = LFQequilimodel1[[i]]$Linf,
                            t0 = -0.01,
                            lwa = a,
                            lwb = b,
                            S50 = LFQequilimodel1[[i]]$SL50,
                            S95 = LFQequilimodel1[[i]]$SL95,
                            selex_input = "length",
                            selex_type = "logistic",
                            M50 = LFQequilimodel1[[i]]$Lmat,
                            M95 = NULL,
                            maturity_input = "length",
                            M = LFQequilimodel1[[i]]$M,
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
  LFQequilimodel1[[i]] <- lfqModify(LFQequilimodel1[[i]], bin_size = 2)
  # LFQequilimodel1[[i]]$catch <- rowSums(LFQequilimodel1[[i]]$catch)
  LFQequilimodel1[[i]]$catch <- as.matrix(LFQequilimodel1[[i]]$catch)
}

#' LF list
LF <- list()
for(i in 1:iters){
  LF[[i]] <- t(LFQequilimodel1[[i]]$catch)
  rownames(LF[[i]]) <- as.numeric(1:12)
  colnames(LF[[i]]) <- LFQequilimodel1[[i]]$midLengths
}

#' years with length data
data_LF <- list()
for (i in 1:iters){
  data_LF[[i]] <- list("years" = 1:12, "LF" = LF[[i]])
}

#' Run LIME
Equili_res_LIME <- list()

# a <- c(which(p %in% "The model is likely not converged"))
for (i in 1:300){
  Equili_res_LIME[[i]] <- run_LIME(modpath = NULL,
                                        lh = lh[[i]],
                                        input_data = data_LF[[i]],
                                        est_sigma = "log_sigma_R",
                                        data_avail = "LC")
}

p <- list()
for (i in 1:iters){
  p[[i]] <- Equili_res_LIME[[i]]$Report$SPR_t
}


saveRDS(Equili_res_LIME, file = "Equili_res_LIME.rds")
