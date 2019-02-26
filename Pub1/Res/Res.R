setwd("D:/DPLBM/Pub1/Res")


# Objects to keep ---------------------------------------------------------


source("D:/DPLBM/Pub1/Res/calc_error.R")
iters <- 300
Med_LBRA <- readRDS("D:/DPLBM/Pub1/Medlived/True/Medmodel_res_LBRA.rds")
Med_LBSPR <- readRDS("D:/DPLBM/Pub1/Medlived/True/Medmodel_res_LBSPR.rds")
Med_LIME <- readRDS("D:/DPLBM/Pub1/Medlived/True/Medmodel_res_LIME.rds")
Med_TB <- readRDS("D:/DPLBM/Pub1/Medlived/True/Medmodel_res_TB.rds")



# Medium (base) ---------------------------------------------------------------------
## SPR ####
## 0.44736870
SPRTrue <- rep(0.44736870, iters)
LIMESPR <- sapply(Med_LIME, function(x) mean(x$Report$SPR_t))
TBSPR <- sapply(Med_TB, function(x) x$currents$curr.SPR)

Med_LBRA_SPR <- calc_error(SPRTrue, Med_LBRA$SPR, iters)
Med_LBSPR_SPR <- calc_error(SPRTrue, Med_LBSPR$SPR, iters)
Med_LIME_SPR <- calc_error(SPRTrue, LIMESPR, iters)
Med_TB_SPR <- calc_error(SPRTrue, TBSPR, iters)

# SPRmsy = 0.183
SPRmsyT <- rep(0.183, iters)
LIMESPRmsy <- sapply(Med_LIME, function(x) x$Derived$SBmsy/x$Derived$SB0)

TBSPR <- lapply(Med_TB, function(x) as.matrix(x$SPR))
TBFM <- lapply(Med_TB, function(x) x$FM_change)
TBFmsy <- sapply(Med_TB, function(x) x$df_Es$Fmax)
TBSPRmsy <- NA
for(i in 1:iters){
  rownames(TBSPR[[i]]) <- TBFM[[i]]
  TBSPRmsy[i] <- TBSPR[[i]][which(rownames(TBSPR[[i]]) == TBFmsy[i])]
  }


Med_LBRA_SPRmsy <- calc_error(SPRmsyT, Med_LBRA$SPRmsy, iters)
Med_LBSPR_SPRmsy <- calc_error(SPRmsyT, Med_LBSPR$SPRmsy, iters)
Med_LIME_SPRmsy <- calc_error(SPRmsyT, LIMESPRmsy, iters)
Med_TB_SPRmsy <- calc_error(SPRmsyT, TBSPRmsy, iters)



# F based reference points ####
# F = 0.128
FMTrue <- rep(0.128, iters)
LIMEFM <- sapply(Med_LIME, function(x) mean(x$Report$F_t))
TBFM <- sapply(Med_TB, function(x) x$currents$curr.F)

Med_LBRA_FM <- calc_error(FMTrue, Med_LBRA$FM, iters)
Med_LIME_FM <- calc_error(FMTrue, LIMEFM, iters)
Med_TB_FM <- calc_error(FMTrue, TBFM, iters)
Med_LBSPR_FM <- calc_error(FMTrue, Med_LBSPR$Fmort, iters)

# Fmsy = 0.306
FmsyTrue <- rep(0.306, iters)
LIMEFmsy <- sapply(Med_LIME, function(x) x$Derived$Fmsy)
TBFmsy <- sapply(Med_TB, function(x) x$df_Es$Fmax)

Med_LBRA_Fmsy <- calc_error(FmsyTrue, Med_LBRA$Fmsy, iters)
Med_LIME_Fmsy <- calc_error(FmsyTrue, LIMEFmsy, iters)
Med_TB_Fmsy <- calc_error(FmsyTrue, TBFmsy, iters)
Med_LBSPR_Fmsy <- calc_error(FmsyTrue, Med_LBSPR$Fmsy, iters)

# F01 = 0.204
F01True <- rep(0.204, iters)
TBF01 <- sapply(Med_TB, function(x) x$df_Es$F01)

Med_TB_F01 <- calc_error(F01True, TBF01, iters)

# F/Fmsy = 0.4183007
FFmsyTrue <- rep(0.4183007, iters)
TBFFmsy <- sapply(Med_TB, function(x) x$currents$curr.F/x$df_Es$Fmax)

Med_LBRA_FFmsy <- calc_error(FFmsyTrue, Med_LBRA$FFmsy, iters)
Med_LIME_FFmsy <- calc_error(FFmsyTrue, LIMEFM/LIMEFmsy, iters)
Med_TB_FFmsy <- calc_error(FFmsyTrue, TBFFmsy, iters)
Med_LBSPR_FFmsy <- calc_error(FFmsyTrue, Med_LBSPR$FFmsy, iters)


# F/F01 = 0.627451
FF01True <- rep(0.627451, iters)
TBFF01 <- sapply(Med_TB, function(x) x$currents$curr.F/x$df_Es$F01)

Med_TB_FF01 <- calc_error(FF01True, TBFF01, iters)







# Short -------------------------------------------------------------------
## SPR ####
# 0.413567
Short_LBRA <- readRDS("D:/DPLBM/Pub1/Shortlived/True/Shortmodel_res_LBRA.rds")
Short_LBSPR <- readRDS("D:/DPLBM/Pub1/Shortlived/True/Shortmodel_res_LBSPR.rds")
Short_LIME <- readRDS("D:/DPLBM/Pub1/Shortlived/True/Shortmodel_res_LIME.rds")
Short_TB <- readRDS("D:/DPLBM/Pub1/Shortlived/True/Shortmodel_res_TB.rds")

SPRTrue <- rep(0.413567, iters)
LIMESPR <- sapply(Short_LIME, function(x) mean(x$Report$SPR_t))
TBSPR <- sapply(Short_TB, function(x) x$currents$curr.SPR)

Short_LBRA_SPR <- calc_error(SPRTrue, Short_LBRA$SPR, iters)
Short_LBSPR_SPR <- calc_error(SPRTrue, Short_LBSPR$SPR, iters)
Short_LIME_SPR <- calc_error(SPRTrue, LIMESPR, iters)
Short_TB_SPR <- calc_error(SPRTrue, TBSPR, iters)


# SPRmsy = 0.148
SPRmsyT <- rep(0.148, iters)
LIMESPRmsy <- sapply(Short_LIME, function(x) x$Derived$SBmsy/x$Derived$SB0)

TBSPR <- lapply(Short_TB, function(x) as.matrix(x$SPR))
TBFM <- lapply(Short_TB, function(x) x$FM_change)
TBFmsy <- sapply(Short_TB, function(x) x$df_Es$Fmax)
TBSPRmsy <- NA
for(i in 1:iters){
  rownames(TBSPR[[i]]) <- TBFM[[i]]
  TBSPRmsy[i] <- TBSPR[[i]][which(rownames(TBSPR[[i]]) == TBFmsy[i])]
}

Short_LBRA_SPRmsy <- calc_error(SPRmsyT, Short_LBRA$SPRmsy, iters)
Short_LBSPR_SPRmsy <- calc_error(SPRmsyT, Short_LBSPR$SPRmsy, iters)
Short_LIME_SPRmsy <- calc_error(SPRmsyT, LIMESPRmsy, iters)
Short_TB_SPRmsy <- calc_error(SPRmsyT, TBSPRmsy, iters)



## F based reference points ####
# F = 0.45
FMTrue <- rep(0.45, iters)
LIMEFM <- sapply(Short_LIME, function(x) mean(x$Report$F_t))
TBFM <- sapply(Short_TB, function(x) x$currents$curr.F)

Short_LBRA_FM <- calc_error(FMTrue, Short_LBRA$FM, iters)
Short_LIME_FM <- calc_error(FMTrue, LIMEFM, iters)
Short_TB_FM <- calc_error(FMTrue, TBFM, iters)
Short_LBSPR_FM <- calc_error(FMTrue, Short_LBSPR$Fmort, iters)

# Fmsy = 1.106
FmsyTrue <- rep(1.106, iters)
LIMEFmsy <- sapply(Short_LIME, function(x) x$Derived$Fmsy)
TBFmsy <- sapply(Short_TB, function(x) x$df_Es$Fmax)

Short_LBRA_Fmsy <- calc_error(FmsyTrue, Short_LBRA$Fmsy, iters)
Short_LIME_Fmsy <- calc_error(FmsyTrue, LIMEFmsy, iters)
Short_TB_Fmsy <- calc_error(FmsyTrue, TBFmsy, iters)
Short_LBSPR_Fmsy <- calc_error(FmsyTrue, Short_LBSPR$Fmsy, iters)

# F01 = 0.673
F01True <- rep(0.673, iters)
TBF01 <- sapply(Short_TB, function(x) x$df_Es$F01)

Short_TB_F01 <- calc_error(F01True, TBF01, iters)

# F/Fmsy = 0.4068716
FFmsyTrue <- rep(0.4068716, iters)
TBFFmsy <- sapply(Short_TB, function(x) x$currents$curr.F/x$df_Es$Fmax)

Short_LBRA_FFmsy <- calc_error(FFmsyTrue, Short_LBRA$FFmsy, iters)
Short_LIME_FFmsy <- calc_error(FFmsyTrue, LIMEFM/LIMEFmsy, iters)
Short_TB_FFmsy <- calc_error(FFmsyTrue, TBFFmsy, iters)
Short_LBSPR_FFmsy <- calc_error(FFmsyTrue, Short_LBSPR$FFmsy, iters)


# F/F01 = 0.6686478
FF01True <- rep(0.6686478, iters)
TBFF01 <- sapply(Short_TB, function(x) x$currents$curr.F/x$df_Es$F01)

Short_TB_FF01 <- calc_error(FF01True, TBFF01, iters)





# Long -------------------------------------------------------------------
## SPR ####
# 0.42935
Long_LBRA <- readRDS("D:/DPLBM/Pub1/Longlived/True/Longmodel_res_LBRA.rds")
Long_LBSPR <- readRDS("D:/DPLBM/Pub1/Longlived/True/Longmodel_res_LBSPR.rds")
Long_LIME <- readRDS("D:/DPLBM/Pub1/Longlived/True/Longmodel_res_LIME.rds")
Long_TB <- readRDS("D:/DPLBM/Pub1/Longlived/True/Longmodel_res_TB.rds")

SPRTrue <- rep(0.42935, iters)
LIMESPR <- sapply(Long_LIME, function(x) x$Report$SPR_t)
TBSPR <- sapply(Long_TB, function(x) x$currents$curr.SPR)

Long_LBRA_SPR <- calc_error(SPRTrue, Long_LBRA$SPR, iters)
Long_LBSPR_SPR <- calc_error(SPRTrue, Long_LBSPR$SPR, iters)
Long_LIME_SPR <- calc_error(SPRTrue, LIMESPR, iters)
Long_TB_SPR <- calc_error(SPRTrue, TBSPR, iters)


# SPRmsy = 0.196
SPRmsyT <- rep(0.196, iters)
LIMESPRmsy <- sapply(Long_LIME, function(x) x$Derived$SBmsy/x$Derived$SB0)

TBSPR <- lapply(Long_TB, function(x) as.matrix(x$SPR))
TBFM <- lapply(Long_TB, function(x) x$FM_change)
TBFmsy <- sapply(Long_TB, function(x) x$df_Es$Fmax)
TBSPRmsy <- NA
for(i in 1:iters){
  rownames(TBSPR[[i]]) <- TBFM[[i]]
  TBSPRmsy[i] <- TBSPR[[i]][which(rownames(TBSPR[[i]]) == TBFmsy[i])]
}


Long_LBRA_SPRmsy <- calc_error(SPRmsyT, Long_LBRA$SPRmsy, iters)
Long_LBSPR_SPRmsy <- calc_error(SPRmsyT, Long_LBSPR$SPRmsy, iters)
Long_LIME_SPRmsy <- calc_error(SPRmsyT, LIMESPRmsy, iters)
Long_TB_SPRmsy <- calc_error(SPRmsyT, TBSPRmsy, iters)


## F based reference points ####
# F = 0.075
FMTrue <- rep(0.075, iters)
LIMEFM <- sapply(Long_LIME, function(x) mean(x$Report$F_t))
TBFM <- sapply(Long_TB, function(x) x$currents$curr.F)

Long_LBRA_FM <- calc_error(FMTrue, Long_LBRA$FM, iters)
Long_LIME_FM <- calc_error(FMTrue, LIMEFM, iters)
Long_TB_FM <- calc_error(FMTrue, TBFM, iters)
Long_LBSPR_FM <- calc_error(FMTrue, Long_LBSPR$Fmort, iters)

# Fmsy = 0.159
FmsyTrue <- rep(0.159, iters)
LIMEFmsy <- sapply(Long_LIME, function(x) x$Derived$Fmsy)
TBFmsy <- sapply(Long_TB, function(x) x$df_Es$Fmax)

Long_LBRA_Fmsy <- calc_error(FmsyTrue, Long_LBRA$Fmsy, iters)
Long_LIME_Fmsy <- calc_error(FmsyTrue, LIMEFmsy, iters)
Long_TB_Fmsy <- calc_error(FmsyTrue, TBFmsy, iters)
Long_LBSPR_Fmsy <- calc_error(FmsyTrue, Long_LBSPR$Fmsy, iters)

# F01 = 0.104
F01True <- rep(0.104, iters)
TBF01 <- sapply(Long_TB, function(x) x$df_Es$F01)

Long_TB_F01 <- calc_error(F01True, TBF01, iters)

# F/Fmsy = 0.4716981
FFmsyTrue <- rep(0.4716981, iters)
TBFFmsy <- sapply(Long_TB, function(x) x$currents$curr.F/x$df_Es$Fmax)

Long_LBRA_FFmsy <- calc_error(FFmsyTrue, Long_LBRA$FFmsy, iters)
Long_LIME_FFmsy <- calc_error(FFmsyTrue, LIMEFM/LIMEFmsy, iters)
Long_TB_FFmsy <- calc_error(FFmsyTrue, TBFFmsy, iters)
Long_LBSPR_FFmsy <- calc_error(FFmsyTrue, Long_LBSPR$FFmsy, iters)


# F/F01 = 0.7211538
FF01True <- rep(0.7211538, iters)
TBFF01 <- sapply(Long_TB, function(x) x$currents$curr.F/x$df_Es$F01)

Long_TB_FF01 <- calc_error(FF01True, TBFF01, iters)





# Underexploit ------------------------------------------------------------
Unex_LBRA <- readRDS("D:/DPLBM/Pub1/Underexploit/True/Unexmodel_res_LBRA.rds")
Unex_LBSPR <- readRDS("D:/DPLBM/Pub1/Underexploit/True/Unexmodel_res_LBSPR.rds")
Unex_LIME <- readRDS("D:/DPLBM/Pub1/Underexploit/True/Unexmodel_res_LIME.rds")
Unex_TB <- readRDS("D:/DPLBM/Pub1/Underexploit/True/Unexmodel_res_TB.rds")

## SPR ####
## 0.6342385
SPRTrue <- rep(0.6342385, iters)
LIMESPR <- sapply(Unex_LIME, function(x) x$Report$SPR_t)
TBSPR <- sapply(Unex_TB, function(x) x$currents$curr.SPR)

Unex_LBRA_SPR <- calc_error(SPRTrue, Unex_LBRA$SPR, iters)
Unex_LBSPR_SPR <- calc_error(SPRTrue, Unex_LBSPR$SPR, iters)
Unex_LIME_SPR <- calc_error(SPRTrue, LIMESPR, iters)
Unex_TB_SPR <- calc_error(SPRTrue, TBSPR, iters)


# SPRmsy = 0.183
SPRmsyT <- rep(0.183, iters)
LIMESPRmsy <- sapply(Unex_LIME, function(x) x$Derived$SBmsy/x$Derived$SB0)

TBSPR <- lapply(Unex_TB, function(x) as.matrix(x$SPR))
TBFM <- lapply(Unex_TB, function(x) x$FM_change)
TBFmsy <- sapply(Unex_TB, function(x) x$df_Es$Fmax)
TBSPRmsy <- NA
for(i in 1:iters){
  rownames(TBSPR[[i]]) <- TBFM[[i]]
  TBSPRmsy[i] <- TBSPR[[i]][which(rownames(TBSPR[[i]]) == TBFmsy[i])]
}


Unex_LBRA_SPRmsy <- calc_error(SPRmsyT, Unex_LBRA$SPRmsy, iters)
Unex_LBSPR_SPRmsy <- calc_error(SPRmsyT, Unex_LBSPR$SPRmsy, iters)
Unex_LIME_SPRmsy <- calc_error(SPRmsyT, LIMESPRmsy, iters)
Unex_TB_SPRmsy <- calc_error(SPRmsyT, TBSPRmsy, iters)



## F based reference points ####
# F = 0.064
FMTrue <- rep(0.064, iters)
LIMEFM <- sapply(Unex_LIME, function(x) mean(x$Report$F_t))
TBFM <- sapply(Unex_TB, function(x) x$currents$curr.F)

Unex_LBRA_FM <- calc_error(FMTrue, Unex_LBRA$FM, iters)
Unex_LIME_FM <- calc_error(FMTrue, LIMEFM, iters)
Unex_TB_FM <- calc_error(FMTrue, TBFM, iters)
Unex_LBSPR_FM <- calc_error(FMTrue, Unex_LBSPR$Fmort, iters)

# Fmsy = 0.306
FmsyTrue <- rep(0.306, iters)
LIMEFmsy <- sapply(Unex_LIME, function(x) x$Derived$Fmsy)
TBFmsy <- sapply(Unex_TB, function(x) x$df_Es$Fmax)

Unex_LBRA_Fmsy <- calc_error(FmsyTrue, Unex_LBRA$Fmsy, iters)
Unex_LIME_Fmsy <- calc_error(FmsyTrue, LIMEFmsy, iters)
Unex_TB_Fmsy <- calc_error(FmsyTrue, TBFmsy, iters)
Unex_LBSPR_Fmsy <- calc_error(FmsyTrue, Unex_LBSPR$Fmsy, iters)

# F01 = 0.204
F01True <- rep(0.204, iters)
TBF01 <- sapply(Unex_TB, function(x) x$df_Es$F01)

Unex_TB_F01 <- calc_error(F01True, TBF01, iters)

# F/Fmsy = 0.2091503
FFmsyTrue <- rep(0.2091503, iters)
TBFFmsy <- sapply(Unex_TB, function(x) x$currents$curr.F/x$df_Es$Fmax)

Unex_LBRA_FFmsy <- calc_error(FFmsyTrue, Unex_LBRA$FFmsy, iters)
Unex_LIME_FFmsy <- calc_error(FFmsyTrue, LIMEFM/LIMEFmsy, iters)
Unex_TB_FFmsy <- calc_error(FFmsyTrue, TBFFmsy, iters)
Unex_LBSPR_FFmsy <- calc_error(FFmsyTrue, Unex_LBSPR$FFmsy, iters)


# F/F01 = 0.3137255
FF01True <- rep(0.3137255, iters)
TBFF01 <- sapply(Unex_TB, function(x) x$currents$curr.F/x$df_Es$F01)

Unex_TB_FF01 <- calc_error(FF01True, TBFF01, iters)





# Overexploit ------------------------------------------------------------
Ovex_LBRA <- readRDS("D:/DPLBM/Pub1/Overexploit/True/Ovexmodel_res_LBRA.rds")
Ovex_LBSPR <- readRDS("D:/DPLBM/Pub1/Overexploit/True/Ovexmodel_res_LBSPR.rds")
Ovex_LIME <- readRDS("D:/DPLBM/Pub1/Overexploit/True/Ovexmodel_res_LIME.rds")
Ovex_TB <- readRDS("D:/DPLBM/Pub1/Overexploit/True/Ovexmodel_res_TB.rds")

## SPR ####
## 0.17171820
SPRTrue <- rep(0.17171820, iters)
LIMESPR <- sapply(Ovex_LIME, function(x) x$Report$SPR_t)
TBSPR <- sapply(Ovex_TB, function(x) x$currents$curr.SPR)

Ovex_LBRA_SPR <- calc_error(SPRTrue, Ovex_LBRA$SPR, iters)
Ovex_LBSPR_SPR <- calc_error(SPRTrue, Ovex_LBSPR$SPR, iters)
Ovex_LIME_SPR <- calc_error(SPRTrue, LIMESPR, iters)
Ovex_TB_SPR <- calc_error(SPRTrue, TBSPR, iters)


# SPRmsy = 0.183
SPRmsyT <- rep(0.183, iters)
LIMESPRmsy <- sapply(Ovex_LIME, function(x) x$Derived$SBmsy/x$Derived$SB0)

TBSPR <- lapply(Ovex_TB, function(x) as.matrix(x$SPR))
TBFM <- lapply(Ovex_TB, function(x) x$FM_change)
TBFmsy <- sapply(Ovex_TB, function(x) x$df_Es$Fmax)
TBSPRmsy <- NA
for(i in 1:iters){
  rownames(TBSPR[[i]]) <- TBFM[[i]]
  TBSPRmsy[i] <- TBSPR[[i]][which(rownames(TBSPR[[i]]) == TBFmsy[i])]
}


Ovex_LBRA_SPRmsy <- calc_error(SPRmsyT, Ovex_LBRA$SPRmsy, iters)
Ovex_LBSPR_SPRmsy <- calc_error(SPRmsyT, Ovex_LBSPR$SPRmsy, iters)
Ovex_LIME_SPRmsy <- calc_error(SPRmsyT, LIMESPRmsy, iters)
Ovex_TB_SPRmsy <- calc_error(SPRmsyT, TBSPRmsy, iters)



# F based reference points ####
# F = 0.278
FMTrue <- rep(0.278, iters)
LIMEFM <- sapply(Ovex_LIME, function(x) mean(x$Report$F_t))
TBFM <- sapply(Ovex_TB, function(x) x$currents$curr.F)

Ovex_LBRA_FM <- calc_error(FMTrue, Ovex_LBRA$FM, iters)
Ovex_LIME_FM <- calc_error(FMTrue, LIMEFM, iters)
Ovex_TB_FM <- calc_error(FMTrue, TBFM, iters)
Ovex_LBSPR_FM <- calc_error(FMTrue, Ovex_LBSPR$Fmort, iters)

# Fmsy = 0.306
FmsyTrue <- rep(0.306, iters)
LIMEFmsy <- sapply(Ovex_LIME, function(x) x$Derived$Fmsy)
TBFmsy <- sapply(Ovex_TB, function(x) x$df_Es$Fmax)

Ovex_LBRA_Fmsy <- calc_error(FmsyTrue, Ovex_LBRA$Fmsy, iters)
Ovex_LIME_Fmsy <- calc_error(FmsyTrue, LIMEFmsy, iters)
Ovex_TB_Fmsy <- calc_error(FmsyTrue, TBFmsy, iters)
Ovex_LBSPR_Fmsy <- calc_error(FmsyTrue, Ovex_LBSPR$Fmsy, iters)

# F01 = 0.204
F01True <- rep(0.204, iters)
TBF01 <- sapply(Ovex_TB, function(x) x$df_Es$F01)

Ovex_TB_F01 <- calc_error(F01True, TBF01, iters)

# F/Fmsy = 0.9084967
FFmsyTrue <- rep(0.9084967, iters)
TBFFmsy <- sapply(Ovex_TB, function(x) x$currents$curr.F/x$df_Es$Fmax)

Ovex_LBRA_FFmsy <- calc_error(FFmsyTrue, Ovex_LBRA$FFmsy, iters)
Ovex_LIME_FFmsy <- calc_error(FFmsyTrue, LIMEFM/LIMEFmsy, iters)
Ovex_TB_FFmsy <- calc_error(FFmsyTrue, TBFFmsy, iters)
Ovex_LBSPR_FFmsy <- calc_error(FFmsyTrue, Ovex_LBSPR$FFmsy, iters)


# F/F01 = 1.362745
FF01True <- rep(1.362745, iters)
TBFF01 <- sapply(Ovex_TB, function(x) x$currents$curr.F/x$df_Es$F01)

Ovex_TB_FF01 <- calc_error(FF01True, TBFF01, iters)






# Error -------------------------------------------------------------------
Error_LBRA <- readRDS("D:/DPLBM/Pub1/Recruiterror/True/Errormodel_res_LBRA.rds")
Error_LBSPR <- readRDS("D:/DPLBM/Pub1/Recruiterror/True/Errormodel_res_LBSPR.rds")
Error_LIME <- readRDS("D:/DPLBM/Pub1/Recruiterror/True/Errormodel_res_LIME.rds")
Error_TB <- readRDS("D:/DPLBM/Pub1/Recruiterror/True/Errormodel_res_TB.rds")

## SPR ####
## 0.4473687
SPRTrue <- rep(0.4473687, iters)
LIMESPR <- sapply(Error_LIME, function(x) x$Report$SPR_t)
TBSPR <- sapply(Error_TB, function(x) x$currents$curr.SPR)

Error_LBRA_SPR <- calc_error(SPRTrue, Error_LBRA$SPR, iters)
Error_LBSPR_SPR <- calc_error(SPRTrue, Error_LBSPR$SPR, iters)
Error_LIME_SPR <- calc_error(SPRTrue, LIMESPR, iters)
Error_TB_SPR <- calc_error(SPRTrue, TBSPR, iters)


# SPRmsy = 0.183
SPRmsyT <- rep(0.183, iters)
LIMESPRmsy <- sapply(Error_LIME, function(x) x$Derived$SBmsy/x$Derived$SB0)

TBSPR <- lapply(Error_TB, function(x) as.matrix(x$SPR))
TBFM <- lapply(Error_TB, function(x) x$FM_change)
TBFmsy <- sapply(Error_TB, function(x) x$df_Es$Fmax)
TBSPRmsy <- NA
for(i in 1:iters){
  rownames(TBSPR[[i]]) <- TBFM[[i]]
  TBSPRmsy[i] <- TBSPR[[i]][which(rownames(TBSPR[[i]]) == TBFmsy[i])]
}


Error_LBRA_SPRmsy <- calc_error(SPRmsyT, Error_LBRA$SPRmsy, iters)
Error_LBSPR_SPRmsy <- calc_error(SPRmsyT, Error_LBSPR$SPRmsy, iters)
Error_LIME_SPRmsy <- calc_error(SPRmsyT, LIMESPRmsy, iters)
Error_TB_SPRmsy <- calc_error(SPRmsyT, TBSPRmsy, iters)


# F based reference points ####
# F = 0.128
FMTrue <- rep(0.128, iters)
LIMEFM <- sapply(Error_LIME, function(x) mean(x$Report$F_t))
TBFM <- sapply(Error_TB, function(x) x$currents$curr.F)

Error_LBRA_FM <- calc_error(FMTrue, Error_LBRA$FM, iters)
Error_LIME_FM <- calc_error(FMTrue, LIMEFM, iters)
Error_TB_FM <- calc_error(FMTrue, TBFM, iters)
Error_LBSPR_FM <- calc_error(FMTrue, Error_LBSPR$Fmort, iters)

# Fmsy = 0.306
FmsyTrue <- rep(0.306, iters)
LIMEFmsy <- sapply(Error_LIME, function(x) x$Derived$Fmsy)
TBFmsy <- sapply(Error_TB, function(x) x$df_Es$Fmax)

Error_LBRA_Fmsy <- calc_error(FmsyTrue, Error_LBRA$Fmsy, iters)
Error_LIME_Fmsy <- calc_error(FmsyTrue, LIMEFmsy, iters)
Error_TB_Fmsy <- calc_error(FmsyTrue, TBFmsy, iters)
Error_LBSPR_Fmsy <- calc_error(FmsyTrue, Error_LBSPR$Fmsy, iters)

# F01 = 0.204
F01True <- rep(0.204, iters)
TBF01 <- sapply(Error_TB, function(x) x$df_Es$F01)

Error_TB_F01 <- calc_error(F01True, TBF01, iters)

# F/Fmsy = 0.4183007
FFmsyTrue <- rep(0.4183007, iters)
TBFFmsy <- sapply(Error_TB, function(x) x$currents$curr.F/x$df_Es$Fmax)

Error_LBRA_FFmsy <- calc_error(FFmsyTrue, Error_LBRA$FFmsy, iters)
Error_LIME_FFmsy <- calc_error(FFmsyTrue, LIMEFM/LIMEFmsy, iters)
Error_TB_FFmsy <- calc_error(FFmsyTrue, TBFFmsy, iters)
Error_LBSPR_FFmsy <- calc_error(FFmsyTrue, Error_LBSPR$FFmsy, iters)


# F/F01 = 0.627451
FF01True <- rep(0.627451, iters)
TBFF01 <- sapply(Error_TB, function(x) x$currents$curr.F/x$df_Es$F01)

Error_TB_FF01 <- calc_error(FF01True, TBFF01, iters)












# AR -------------------------------------------------------------------
AR_LBRA <- readRDS("D:/DPLBM/Pub1/AR/True/ARmodel_res_LBRA.rds")
AR_LBSPR <- readRDS("D:/DPLBM/Pub1/AR/True/ARmodel_res_LBSPR.rds")
AR_LIME <- readRDS("D:/DPLBM/Pub1/AR/True/ARmodel_res_LIME.rds")
AR_TB <- readRDS("D:/DPLBM/Pub1/AR/True/ARmodel_res_TB.rds")

## SPR ####
## 0.42063743
SPRTrue <- rep(0.42063743, iters)
LIMESPR <- sapply(AR_LIME, function(x) x$Report$SPR_t)
TBSPR <- sapply(AR_TB, function(x) x$currents$curr.SPR)

AR_LBRA_SPR <- calc_error(SPRTrue, AR_LBRA$SPR, iters)
AR_LBSPR_SPR <- calc_error(SPRTrue, AR_LBSPR$SPR, iters)
AR_LIME_SPR <- calc_error(SPRTrue, LIMESPR, iters)
AR_TB_SPR <- calc_error(SPRTrue, TBSPR, iters)



# SPRmsy = 0.217
SPRmsyT <- rep(0.217, iters)
LIMESPRmsy <- sapply(AR_LIME, function(x) x$Derived$SBmsy/x$Derived$SB0)

TBSPR <- lapply(AR_TB, function(x) as.matrix(x$SPR))
TBFM <- lapply(AR_TB, function(x) x$FM_change)
TBFmsy <- sapply(AR_TB, function(x) x$df_Es$Fmax)
TBSPRmsy <- NA
for(i in 1:iters){
  rownames(TBSPR[[i]]) <- TBFM[[i]]
  TBSPRmsy[i] <- TBSPR[[i]][which(rownames(TBSPR[[i]]) == TBFmsy[i])]
}


AR_LBRA_SPRmsy <- calc_error(SPRmsyT, AR_LBRA$SPRmsy, iters)
AR_LBSPR_SPRmsy <- calc_error(SPRmsyT, AR_LBSPR$SPRmsy, iters)
AR_LIME_SPRmsy <- calc_error(SPRmsyT, LIMESPRmsy, iters)
AR_TB_SPRmsy <- calc_error(SPRmsyT, TBSPRmsy, iters)


# F based reference points ####
# F = 0.128
FMTrue <- rep(0.128, iters)
LIMEFM <- sapply(AR_LIME, function(x) mean(x$Report$F_t))
TBFM <- sapply(AR_TB, function(x) x$currents$curr.F)

AR_LBRA_FM <- calc_error(FMTrue, AR_LBRA$FM, iters)
AR_LIME_FM <- calc_error(FMTrue, LIMEFM, iters)
AR_TB_FM <- calc_error(FMTrue, TBFM, iters)
AR_LBSPR_FM <- calc_error(FMTrue, AR_LBSPR$Fmort, iters)

# Fmsy = 0.247
FmsyTrue <- rep(0.247, iters)
LIMEFmsy <- sapply(AR_LIME, function(x) x$Derived$Fmsy)
TBFmsy <- sapply(AR_TB, function(x) x$df_Es$Fmax)

AR_LBRA_Fmsy <- calc_error(FmsyTrue, AR_LBRA$Fmsy, iters)
AR_LIME_Fmsy <- calc_error(FmsyTrue, LIMEFmsy, iters)
AR_TB_Fmsy <- calc_error(FmsyTrue, TBFmsy, iters)
AR_LBSPR_Fmsy <- calc_error(FmsyTrue, AR_LBSPR$Fmsy, iters)

# F01 = 0.176
F01True <- rep(0.176, iters)
TBF01 <- sapply(AR_TB, function(x) x$df_Es$F01)

AR_TB_F01 <- calc_error(F01True, TBF01, iters)

# F/Fmsy = 0.5182186
FFmsyTrue <- rep(0.5182186, iters)
TBFFmsy <- sapply(AR_TB, function(x) x$currents$curr.F/x$df_Es$Fmax)

AR_LBRA_FFmsy <- calc_error(FFmsyTrue, AR_LBRA$FFmsy, iters)
AR_LIME_FFmsy <- calc_error(FFmsyTrue, LIMEFM/LIMEFmsy, iters)
AR_TB_FFmsy <- calc_error(FFmsyTrue, TBFFmsy, iters)
AR_LBSPR_FFmsy <- calc_error(FFmsyTrue, AR_LBSPR$FFmsy, iters)


# F/F01 = 0.7191011
FF01True <- rep(0.7191011, iters)
TBFF01 <- sapply(AR_TB, function(x) x$currents$curr.F/x$df_Es$F01)

AR_TB_FF01 <- calc_error(FF01True, TBFF01, iters)









# Relative error ----------------------------------------------------------

# Life history
LH_SPR <- data.frame(c(Med_TB_SPR$rel_error, Med_LBSPR_SPR$rel_error, Med_LIME_SPR$rel_error, Med_LBRA_SPR$rel_error, 
                       Short_TB_SPR$rel_error, Short_LBSPR_SPR$rel_error, Short_LIME_SPR$rel_error, Short_LBRA_SPR$rel_error,
                       Long_TB_SPR$rel_error, Long_LBSPR_SPR$rel_error, Long_LIME_SPR$rel_error, Long_LBRA_SPR$rel_error),
                       c(rep("Med", 1200), rep("Short", 1200), rep("Long",1200)),
                       c(rep(rep(c("TB", "LBSPR", "LIME", "LBRA"), each = 300), 3)))
colnames(LH_SPR) <- c("Relative_error", "Scenario", "Method")
LH_SPR$Scenario <- factor(LH_SPR$Scenario, levels = c("Med", "Short", "Long"))
LH_SPR$Method <- factor(LH_SPR$Method, levels = c("TB", "LBSPR", "LIME", "LBRA"))
saveRDS(LH_SPR, file = "files/LH_SPR.rds")




LH_SPRmsy <- data.frame(c(Med_TB_SPRmsy$rel_error, Med_LBSPR_SPRmsy$rel_error, Med_LIME_SPRmsy$rel_error, Med_LBRA_SPRmsy$rel_error, 
                         Short_TB_SPRmsy$rel_error, Short_LBSPR_SPRmsy$rel_error, Short_LIME_SPRmsy$rel_error, Short_LBRA_SPRmsy$rel_error,
                         Long_TB_SPRmsy$rel_error, Long_LBSPR_SPRmsy$rel_error, Long_LIME_SPRmsy$rel_error, Long_LBRA_SPRmsy$rel_error),
                       c(rep("Med", 1200), rep("Short", 1200), rep("Long",1200)),
                       c(rep(rep(c("TB", "LBSPR", "LIME", "LBRA"), each = 300), 3)))
colnames(LH_SPRmsy) <- c("Relative_error", "Scenario", "Method")
LH_SPRmsy$Scenario <- factor(LH_SPRmsy$Scenario, levels = c("Med", "Short", "Long"))
LH_SPRmsy$Method <- factor(LH_SPRmsy$Method, levels = c("TB", "LBSPR", "LIME", "LBRA"))
saveRDS(LH_SPRmsy, file = "files/LH_SPRmsy.rds")




LH_FFmsy <- data.frame(c(Med_TB_FFmsy$rel_error, Med_LBSPR_FFmsy$rel_error, Med_LIME_FFmsy$rel_error, Med_LBRA_FFmsy$rel_error, 
                       Short_TB_FFmsy$rel_error, Short_LBSPR_FFmsy$rel_error, Short_LIME_FFmsy$rel_error, Short_LBRA_FFmsy$rel_error,
                       Long_TB_FFmsy$rel_error, Long_LBSPR_FFmsy$rel_error, Long_LIME_FFmsy$rel_error, Long_LBRA_FFmsy$rel_error),
                     c(rep("Med", 1200), rep("Short", 1200), rep("Long",1200)),
                     c(rep(rep(c("TB", "LBSPR", "LIME", "LBRA"), each = 300), 3)))
colnames(LH_FFmsy) <- c("Relative_error", "Scenario", "Method")
LH_FFmsy$Scenario <- factor(LH_FFmsy$Scenario, levels = c("Med", "Short", "Long"))
LH_FFmsy$Method <- factor(LH_FFmsy$Method, levels = c("TB", "LBSPR", "LIME", "LBRA"))
saveRDS(LH_FFmsy, file = "files/LH_FFmsy.rds")




LH_FM <- data.frame(c(Med_TB_FM$rel_error, Med_LBSPR_FM$rel_error, Med_LIME_FM$rel_error, Med_LBRA_FM$rel_error, 
                         Short_TB_FM$rel_error, Short_LBSPR_FM$rel_error, Short_LIME_FM$rel_error, Short_LBRA_FM$rel_error,
                         Long_TB_FM$rel_error, Long_LBSPR_FM$rel_error, Long_LIME_FM$rel_error, Long_LBRA_FM$rel_error),
                       c(rep("Med", 1200), rep("Short", 1200), rep("Long",1200)),
                       c(rep(rep(c("TB", "LBSPR", "LIME", "LBRA"), each = 300), 3)))
colnames(LH_FM) <- c("Relative_error", "Scenario", "Method")
LH_FM$Scenario <- factor(LH_FM$Scenario, levels = c("Med", "Short", "Long"))
LH_FM$Method <- factor(LH_FM$Method, levels = c("TB", "LBSPR", "LIME", "LBRA"))
saveRDS(LH_FM, file = "files/LH_FM.rds")




LH_Fmsy <- data.frame(c(Med_TB_Fmsy$rel_error, Med_LBSPR_Fmsy$rel_error, Med_LIME_Fmsy$rel_error, Med_LBRA_Fmsy$rel_error, 
                      Short_TB_Fmsy$rel_error, Short_LBSPR_Fmsy$rel_error, Short_LIME_Fmsy$rel_error, Short_LBRA_Fmsy$rel_error,
                      Long_TB_Fmsy$rel_error, Long_LBSPR_Fmsy$rel_error, Long_LIME_Fmsy$rel_error, Long_LBRA_Fmsy$rel_error),
                    c(rep("Med", 1200), rep("Short", 1200), rep("Long",1200)),
                    c(rep(rep(c("TB", "LBSPR", "LIME", "LBRA"), each = 300), 3)))
colnames(LH_Fmsy) <- c("Relative_error", "Scenario", "Method")
LH_Fmsy$Scenario <- factor(LH_Fmsy$Scenario, levels = c("Med", "Short", "Long"))
LH_Fmsy$Method <- factor(LH_Fmsy$Method, levels = c("TB", "LBSPR", "LIME", "LBRA"))
saveRDS(LH_Fmsy, file = "files/LH_Fmsy.rds")







# Exploitation
Exp_SPR <- data.frame(c(Med_TB_SPR$rel_error, Med_LBSPR_SPR$rel_error, Med_LIME_SPR$rel_error, Med_LBRA_SPR$rel_error, 
                       Unex_TB_SPR$rel_error, Unex_LBSPR_SPR$rel_error, Unex_LIME_SPR$rel_error, Unex_LBRA_SPR$rel_error,
                       Ovex_TB_SPR$rel_error, Ovex_LBSPR_SPR$rel_error, Ovex_LIME_SPR$rel_error, Ovex_LBRA_SPR$rel_error),
                     c(rep("Targ", 1200), rep("Unex", 1200), rep("Ovex",1200)),
                     c(rep(rep(c("TB", "LBSPR", "LIME", "LBRA"), each = 300), 3)))
colnames(Exp_SPR) <- c("Relative_error", "Scenario", "Method")
Exp_SPR$Scenario <- factor(Exp_SPR$Scenario, levels = c("Targ", "Unex", "Ovex"))
Exp_SPR$Method <- factor(Exp_SPR$Method, levels = c("TB", "LBSPR", "LIME", "LBRA"))
saveRDS(Exp_SPR, file = "files/Exp_SPR.rds")




Exp_SPRmsy <- data.frame(c(Med_TB_SPRmsy$rel_error, Med_LBSPR_SPRmsy$rel_error, Med_LIME_SPRmsy$rel_error, Med_LBRA_SPRmsy$rel_error, 
                          Unex_TB_SPRmsy$rel_error, Unex_LBSPR_SPRmsy$rel_error, Unex_LIME_SPRmsy$rel_error, Unex_LBRA_SPRmsy$rel_error,
                          Ovex_TB_SPRmsy$rel_error, Ovex_LBSPR_SPRmsy$rel_error, Ovex_LIME_SPRmsy$rel_error, Ovex_LBRA_SPRmsy$rel_error),
                        c(rep("Targ", 1200), rep("Unex", 1200), rep("Ovex",1200)),
                        c(rep(rep(c("TB", "LBSPR", "LIME", "LBRA"), each = 300), 3)))
colnames(Exp_SPRmsy) <- c("Relative_error", "Scenario", "Method")
Exp_SPRmsy$Scenario <- factor(Exp_SPRmsy$Scenario, levels = c("Targ", "Unex", "Ovex"))
Exp_SPRmsy$Method <- factor(Exp_SPRmsy$Method, levels = c("TB", "LBSPR", "LIME", "LBRA"))
saveRDS(Exp_SPRmsy, file = "files/Exp_SPRmsy.rds")




Exp_FFmsy <- data.frame(c(Med_TB_FFmsy$rel_error, Med_LBSPR_FFmsy$rel_error, Med_LIME_FFmsy$rel_error, Med_LBRA_FFmsy$rel_error, 
                         Unex_TB_FFmsy$rel_error, Unex_LBSPR_FFmsy$rel_error, Unex_LIME_FFmsy$rel_error, Unex_LBRA_FFmsy$rel_error,
                         Ovex_TB_FFmsy$rel_error, Ovex_LBSPR_FFmsy$rel_error, Ovex_LIME_FFmsy$rel_error, Ovex_LBRA_FFmsy$rel_error),
                       c(rep("Targ", 1200), rep("Unex", 1200), rep("Ovex",1200)),
                       c(rep(rep(c("TB", "LBSPR", "LIME", "LBRA"), each = 300), 3)))
colnames(Exp_FFmsy) <- c("Relative_error", "Scenario", "Method")
Exp_FFmsy$Scenario <- factor(Exp_FFmsy$Scenario, levels = c("Targ", "Unex", "Ovex"))
Exp_FFmsy$Method <- factor(Exp_FFmsy$Method, levels = c("TB", "LBSPR", "LIME", "LBRA"))
saveRDS(Exp_FFmsy, file = "files/Exp_FFmsy.rds")




Exp_FM <- data.frame(c(Med_TB_FM$rel_error, Med_LBSPR_FM$rel_error, Med_LIME_FM$rel_error, Med_LBRA_FM$rel_error, 
                      Unex_TB_FM$rel_error, Unex_LBSPR_FM$rel_error, Unex_LIME_FM$rel_error, Unex_LBRA_FM$rel_error,
                      Ovex_TB_FM$rel_error, Ovex_LBSPR_FM$rel_error, Ovex_LIME_FM$rel_error, Ovex_LBRA_FM$rel_error),
                    c(rep("Targ", 1200), rep("Unex", 1200), rep("Ovex",1200)),
                    c(rep(rep(c("TB", "LBSPR", "LIME", "LBRA"), each = 300), 3)))
colnames(Exp_FM) <- c("Relative_error", "Scenario", "Method")
Exp_FM$Scenario <- factor(Exp_FM$Scenario, levels = c("Targ", "Unex", "Ovex"))
Exp_FM$Method <- factor(Exp_FM$Method, levels = c("TB", "LBSPR", "LIME", "LBRA"))
saveRDS(Exp_FM, file = "files/Exp_FM.rds")




Exp_Fmsy <- data.frame(c(Med_TB_Fmsy$rel_error, Med_LBSPR_Fmsy$rel_error, Med_LIME_Fmsy$rel_error, Med_LBRA_Fmsy$rel_error, 
                        Unex_TB_Fmsy$rel_error, Unex_LBSPR_Fmsy$rel_error, Unex_LIME_Fmsy$rel_error, Unex_LBRA_Fmsy$rel_error,
                        Ovex_TB_Fmsy$rel_error, Ovex_LBSPR_Fmsy$rel_error, Ovex_LIME_Fmsy$rel_error, Ovex_LBRA_Fmsy$rel_error),
                      c(rep("Targ", 1200), rep("Unex", 1200), rep("Ovex",1200)),
                      c(rep(rep(c("TB", "LBSPR", "LIME", "LBRA"), each = 300), 3)))
colnames(Exp_Fmsy) <- c("Relative_error", "Scenario", "Method")
Exp_Fmsy$Scenario <- factor(Exp_Fmsy$Scenario, levels = c("Targ", "Unex", "Ovex"))
Exp_Fmsy$Method <- factor(Exp_Fmsy$Method, levels = c("TB", "LBSPR", "LIME", "LBRA"))
saveRDS(Exp_Fmsy, file = "files/Exp_Fmsy.rds")






# Recruitment
Rec_SPR <- data.frame(c(Med_TB_SPR$rel_error, Med_LBSPR_SPR$rel_error, Med_LIME_SPR$rel_error, Med_LBRA_SPR$rel_error, 
                       Error_TB_SPR$rel_error, Error_LBSPR_SPR$rel_error, Error_LIME_SPR$rel_error, Error_LBRA_SPR$rel_error,
                       AR_TB_SPR$rel_error, AR_LBSPR_SPR$rel_error, AR_LIME_SPR$rel_error, AR_LBRA_SPR$rel_error),
                     c(rep("None", 1200), rep("Error", 1200), rep("AR",1200)),
                     c(rep(rep(c("TB", "LBSPR", "LIME", "LBRA"), each = 300), 3)))
colnames(Rec_SPR) <- c("Relative_error", "Scenario", "Method")
Rec_SPR$Scenario <- factor(Rec_SPR$Scenario, levels = c("None", "Error", "AR"))
Rec_SPR$Method <- factor(Rec_SPR$Method, levels = c("TB", "LBSPR", "LIME", "LBRA"))
saveRDS(Rec_SPR, file = "files/Rec_SPR.rds")




Rec_SPRmsy <- data.frame(c(Med_TB_SPRmsy$rel_error, Med_LBSPR_SPRmsy$rel_error, Med_LIME_SPRmsy$rel_error, Med_LBRA_SPRmsy$rel_error, 
                          Error_TB_SPRmsy$rel_error, Error_LBSPR_SPRmsy$rel_error, Error_LIME_SPRmsy$rel_error, Error_LBRA_SPRmsy$rel_error,
                          AR_TB_SPRmsy$rel_error, AR_LBSPR_SPRmsy$rel_error, AR_LIME_SPRmsy$rel_error, AR_LBRA_SPRmsy$rel_error),
                        c(rep("None", 1200), rep("Error", 1200), rep("AR",1200)),
                        c(rep(rep(c("TB", "LBSPR", "LIME", "LBRA"), each = 300), 3)))
colnames(Rec_SPRmsy) <- c("Relative_error", "Scenario", "Method")
Rec_SPRmsy$Scenario <- factor(Rec_SPRmsy$Scenario, levels = c("None", "Error", "AR"))
Rec_SPRmsy$Method <- factor(Rec_SPRmsy$Method, levels = c("TB", "LBSPR", "LIME", "LBRA"))
saveRDS(Rec_SPRmsy, file = "files/Rec_SPRmsy.rds")




Rec_FFmsy <- data.frame(c(Med_TB_FFmsy$rel_error, Med_LBSPR_FFmsy$rel_error, Med_LIME_FFmsy$rel_error, Med_LBRA_FFmsy$rel_error, 
                         Error_TB_FFmsy$rel_error, Error_LBSPR_FFmsy$rel_error, Error_LIME_FFmsy$rel_error, Error_LBRA_FFmsy$rel_error,
                         AR_TB_FFmsy$rel_error, AR_LBSPR_FFmsy$rel_error, AR_LIME_FFmsy$rel_error, AR_LBRA_FFmsy$rel_error),
                       c(rep("None", 1200), rep("Error", 1200), rep("AR",1200)),
                       c(rep(rep(c("TB", "LBSPR", "LIME", "LBRA"), each = 300), 3)))
colnames(Rec_FFmsy) <- c("Relative_error", "Scenario", "Method")
Rec_FFmsy$Scenario <- factor(Rec_FFmsy$Scenario, levels = c("None", "Error", "AR"))
Rec_FFmsy$Method <- factor(Rec_FFmsy$Method, levels = c("TB", "LBSPR", "LIME", "LBRA"))
saveRDS(Rec_FFmsy, file = "files/Rec_FFmsy.rds")




Rec_FM <- data.frame(c(Med_TB_FM$rel_error, Med_LBSPR_FM$rel_error, Med_LIME_FM$rel_error, Med_LBRA_FM$rel_error, 
                      Error_TB_FM$rel_error, Error_LBSPR_FM$rel_error, Error_LIME_FM$rel_error, Error_LBRA_FM$rel_error,
                      AR_TB_FM$rel_error, AR_LBSPR_FM$rel_error, AR_LIME_FM$rel_error, AR_LBRA_FM$rel_error),
                    c(rep("None", 1200), rep("Error", 1200), rep("AR",1200)),
                    c(rep(rep(c("TB", "LBSPR", "LIME", "LBRA"), each = 300), 3)))
colnames(Rec_FM) <- c("Relative_error", "Scenario", "Method")
Rec_FM$Scenario <- factor(Rec_FM$Scenario, levels = c("None", "Error", "AR"))
Rec_FM$Method <- factor(Rec_FM$Method, levels = c("TB", "LBSPR", "LIME", "LBRA"))
saveRDS(Rec_FM, file = "files/Rec_FM.rds")




Rec_Fmsy <- data.frame(c(Med_TB_Fmsy$rel_error, Med_LBSPR_Fmsy$rel_error, Med_LIME_Fmsy$rel_error, Med_LBRA_Fmsy$rel_error, 
                        Error_TB_Fmsy$rel_error, Error_LBSPR_Fmsy$rel_error, Error_LIME_Fmsy$rel_error, Error_LBRA_Fmsy$rel_error,
                        AR_TB_Fmsy$rel_error, AR_LBSPR_Fmsy$rel_error, AR_LIME_Fmsy$rel_error, AR_LBRA_Fmsy$rel_error),
                      c(rep("None", 1200), rep("Error", 1200), rep("AR",1200)),
                      c(rep(rep(c("TB", "LBSPR", "LIME", "LBRA"), each = 300), 3)))
colnames(Rec_Fmsy) <- c("Relative_error", "Scenario", "Method")
Rec_Fmsy$Scenario <- factor(Rec_Fmsy$Scenario, levels = c("None", "Error", "AR"))
Rec_Fmsy$Method <- factor(Rec_Fmsy$Method, levels = c("TB", "LBSPR", "LIME", "LBRA"))
saveRDS(Rec_Fmsy, file = "files/Rec_Fmsy.rds")




# MRE and MARE ------------------------------------------------------------

# LH
LH_SPR.err <- list()
LH_SPR.err$TB <- cbind(c(Med_TB_SPR$MRE, Short_TB_SPR$MRE, Long_TB_SPR$MRE), c(Med_TB_SPR$MARE, Short_TB_SPR$MARE, Long_TB_SPR$MARE))
LH_SPR.err$LBSPR <- cbind(c(Med_LBSPR_SPR$MRE, Short_LBSPR_SPR$MRE, Long_LBSPR_SPR$MRE), c(Med_LBSPR_SPR$MARE, Short_LBSPR_SPR$MARE, Long_LBSPR_SPR$MARE))
LH_SPR.err$LIME <- cbind(c(Med_LIME_SPR$MRE, Short_LIME_SPR$MRE, Long_LIME_SPR$MRE), c(Med_LIME_SPR$MARE, Short_LIME_SPR$MARE, Long_LIME_SPR$MARE))
LH_SPR.err$LBRA <- cbind(c(Med_LBRA_SPR$MRE, Short_LBRA_SPR$MRE, Long_LBRA_SPR$MRE), c(Med_LBRA_SPR$MARE, Short_LBRA_SPR$MARE, Long_LBRA_SPR$MARE))

for(i in 1:4){
  colnames(LH_SPR.err[[i]]) <- c("MRE", "MARE")
  rownames(LH_SPR.err[[i]]) <- c("Med", "Short", "Long")
}

saveRDS(LH_SPR.err, file = "files/LH_SPR_err.rds")




LH_SPRmsy.err <- list()
LH_SPRmsy.err$TB <- cbind(c(Med_TB_SPRmsy$MRE, Short_TB_SPRmsy$MRE, Long_TB_SPRmsy$MRE), c(Med_TB_SPRmsy$MARE, Short_TB_SPRmsy$MARE, Long_TB_SPRmsy$MARE))
LH_SPRmsy.err$LBSPR <- cbind(c(Med_LBSPR_SPRmsy$MRE, Short_LBSPR_SPRmsy$MRE, Long_LBSPR_SPRmsy$MRE), c(Med_LBSPR_SPRmsy$MARE, Short_LBSPR_SPRmsy$MARE, Long_LBSPR_SPRmsy$MARE))
LH_SPRmsy.err$LIME <- cbind(c(Med_LIME_SPRmsy$MRE, Short_LIME_SPRmsy$MRE, Long_LIME_SPRmsy$MRE), c(Med_LIME_SPRmsy$MARE, Short_LIME_SPRmsy$MARE, Long_LIME_SPRmsy$MARE))
LH_SPRmsy.err$LBRA <- cbind(c(Med_LBRA_SPRmsy$MRE, Short_LBRA_SPRmsy$MRE, Long_LBRA_SPRmsy$MRE), c(Med_LBRA_SPRmsy$MARE, Short_LBRA_SPRmsy$MARE, Long_LBRA_SPRmsy$MARE))

for(i in 1:4){
  colnames(LH_SPRmsy.err[[i]]) <- c("MRE", "MARE")
  rownames(LH_SPRmsy.err[[i]]) <- c("Med", "Short", "Long")
}

saveRDS(LH_SPRmsy.err, file = "files/LH_SPRmsy_err.rds")




LH_FFmsy.err <- list()
LH_FFmsy.err$TB <- cbind(c(Med_TB_FFmsy$MRE, Short_TB_FFmsy$MRE, Long_TB_FFmsy$MRE), c(Med_TB_FFmsy$MARE, Short_TB_FFmsy$MARE, Long_TB_FFmsy$MARE))
LH_FFmsy.err$LBSPR <- cbind(c(Med_LBSPR_FFmsy$MRE, Short_LBSPR_FFmsy$MRE, Long_LBSPR_FFmsy$MRE), c(Med_LBSPR_FFmsy$MARE, Short_LBSPR_FFmsy$MARE, Long_LBSPR_FFmsy$MARE))
LH_FFmsy.err$LIME <- cbind(c(Med_LIME_FFmsy$MRE, Short_LIME_FFmsy$MRE, Long_LIME_FFmsy$MRE), c(Med_LIME_FFmsy$MARE, Short_LIME_FFmsy$MARE, Long_LIME_FFmsy$MARE))
LH_FFmsy.err$LBRA <- cbind(c(Med_LBRA_FFmsy$MRE, Short_LBRA_FFmsy$MRE, Long_LBRA_FFmsy$MRE), c(Med_LBRA_FFmsy$MARE, Short_LBRA_FFmsy$MARE, Long_LBRA_FFmsy$MARE))

for(i in 1:4){
  colnames(LH_FFmsy.err[[i]]) <- c("MRE", "MARE")
  rownames(LH_FFmsy.err[[i]]) <- c("Med", "Short", "Long")
}

saveRDS(LH_FFmsy.err, file = "files/LH_FFmsy_err.rds")




LH_FM.err <- list()
LH_FM.err$TB <- cbind(c(Med_TB_FM$MRE, Short_TB_FM$MRE, Long_TB_FM$MRE), c(Med_TB_FM$MARE, Short_TB_FM$MARE, Long_TB_FM$MARE))
LH_FM.err$LBSPR <- cbind(c(Med_LBSPR_FM$MRE, Short_LBSPR_FM$MRE, Long_LBSPR_FM$MRE), c(Med_LBSPR_FM$MARE, Short_LBSPR_FM$MARE, Long_LBSPR_FM$MARE))
LH_FM.err$LIME <- cbind(c(Med_LIME_FM$MRE, Short_LIME_FM$MRE, Long_LIME_FM$MRE), c(Med_LIME_FM$MARE, Short_LIME_FM$MARE, Long_LIME_FM$MARE))
LH_FM.err$LBRA <- cbind(c(Med_LBRA_FM$MRE, Short_LBRA_FM$MRE, Long_LBRA_FM$MRE), c(Med_LBRA_FM$MARE, Short_LBRA_FM$MARE, Long_LBRA_FM$MARE))

for(i in 1:4){
  colnames(LH_FM.err[[i]]) <- c("MRE", "MARE")
  rownames(LH_FM.err[[i]]) <- c("Med", "Short", "Long")
}

saveRDS(LH_FM.err, file = "files/LH_FM_err.rds")





LH_Fmsy.err <- list()
LH_Fmsy.err$TB <- cbind(c(Med_TB_Fmsy$MRE, Short_TB_Fmsy$MRE, Long_TB_Fmsy$MRE), c(Med_TB_Fmsy$MARE, Short_TB_Fmsy$MARE, Long_TB_Fmsy$MARE))
LH_Fmsy.err$LBSPR <- cbind(c(Med_LBSPR_Fmsy$MRE, Short_LBSPR_Fmsy$MRE, Long_LBSPR_Fmsy$MRE), c(Med_LBSPR_Fmsy$MARE, Short_LBSPR_Fmsy$MARE, Long_LBSPR_Fmsy$MARE))
LH_Fmsy.err$LIME <- cbind(c(Med_LIME_Fmsy$MRE, Short_LIME_Fmsy$MRE, Long_LIME_Fmsy$MRE), c(Med_LIME_Fmsy$MARE, Short_LIME_Fmsy$MARE, Long_LIME_Fmsy$MARE))
LH_Fmsy.err$LBRA <- cbind(c(Med_LBRA_Fmsy$MRE, Short_LBRA_Fmsy$MRE, Long_LBRA_Fmsy$MRE), c(Med_LBRA_Fmsy$MARE, Short_LBRA_Fmsy$MARE, Long_LBRA_Fmsy$MARE))

for(i in 1:4){
  colnames(LH_Fmsy.err[[i]]) <- c("MRE", "MARE")
  rownames(LH_Fmsy.err[[i]]) <- c("Med", "Short", "Long")
}

saveRDS(LH_Fmsy.err, file = "files/LH_Fmsy_err.rds")




# Exp
Exp_SPR.err <- list()
Exp_SPR.err$TB <- cbind(c(Med_TB_SPR$MRE, Unex_TB_SPR$MRE, Ovex_TB_SPR$MRE), c(Med_TB_SPR$MARE, Unex_TB_SPR$MARE, Ovex_TB_SPR$MARE))
Exp_SPR.err$LBSPR <- cbind(c(Med_LBSPR_SPR$MRE, Unex_LBSPR_SPR$MRE, Ovex_LBSPR_SPR$MRE), c(Med_LBSPR_SPR$MARE, Unex_LBSPR_SPR$MARE, Ovex_LBSPR_SPR$MARE))
Exp_SPR.err$LIME <- cbind(c(Med_LIME_SPR$MRE, Unex_LIME_SPR$MRE, Ovex_LIME_SPR$MRE), c(Med_LIME_SPR$MARE, Unex_LIME_SPR$MARE, Ovex_LIME_SPR$MARE))
Exp_SPR.err$LBRA <- cbind(c(Med_LBRA_SPR$MRE, Unex_LBRA_SPR$MRE, Ovex_LBRA_SPR$MRE), c(Med_LBRA_SPR$MARE, Unex_LBRA_SPR$MARE, Ovex_LBRA_SPR$MARE))

for(i in 1:4){
  colnames(Exp_SPR.err[[i]]) <- c("MRE", "MARE")
  rownames(Exp_SPR.err[[i]]) <- c("Targ", "Unex", "Ovex")
}

saveRDS(Exp_SPR.err, file = "files/Exp_SPR_err.rds")




Exp_SPRmsy.err <- list()
Exp_SPRmsy.err$TB <- cbind(c(Med_TB_SPRmsy$MRE, Unex_TB_SPRmsy$MRE, Ovex_TB_SPRmsy$MRE), c(Med_TB_SPRmsy$MARE, Unex_TB_SPRmsy$MARE, Ovex_TB_SPRmsy$MARE))
Exp_SPRmsy.err$LBSPR <- cbind(c(Med_LBSPR_SPRmsy$MRE, Unex_LBSPR_SPRmsy$MRE, Ovex_LBSPR_SPRmsy$MRE), c(Med_LBSPR_SPRmsy$MARE, Unex_LBSPR_SPRmsy$MARE, Ovex_LBSPR_SPRmsy$MARE))
Exp_SPRmsy.err$LIME <- cbind(c(Med_LIME_SPRmsy$MRE, Unex_LIME_SPRmsy$MRE, Ovex_LIME_SPRmsy$MRE), c(Med_LIME_SPRmsy$MARE, Unex_LIME_SPRmsy$MARE, Ovex_LIME_SPRmsy$MARE))
Exp_SPRmsy.err$LBRA <- cbind(c(Med_LBRA_SPRmsy$MRE, Unex_LBRA_SPRmsy$MRE, Ovex_LBRA_SPRmsy$MRE), c(Med_LBRA_SPRmsy$MARE, Unex_LBRA_SPRmsy$MARE, Ovex_LBRA_SPRmsy$MARE))

for(i in 1:4){
  colnames(Exp_SPRmsy.err[[i]]) <- c("MRE", "MARE")
  rownames(Exp_SPRmsy.err[[i]]) <- c("Targ", "Unex", "Ovex")
}

saveRDS(Exp_SPRmsy.err, file = "files/Exp_SPRmsy_err.rds")




Exp_FFmsy.err <- list()
Exp_FFmsy.err$TB <- cbind(c(Med_TB_FFmsy$MRE, Unex_TB_FFmsy$MRE, Ovex_TB_FFmsy$MRE), c(Med_TB_FFmsy$MARE, Unex_TB_FFmsy$MARE, Ovex_TB_FFmsy$MARE))
Exp_FFmsy.err$LBSPR <- cbind(c(Med_LBSPR_FFmsy$MRE, Unex_LBSPR_FFmsy$MRE, Ovex_LBSPR_FFmsy$MRE), c(Med_LBSPR_FFmsy$MARE, Unex_LBSPR_FFmsy$MARE, Ovex_LBSPR_FFmsy$MARE))
Exp_FFmsy.err$LIME <- cbind(c(Med_LIME_FFmsy$MRE, Unex_LIME_FFmsy$MRE, Ovex_LIME_FFmsy$MRE), c(Med_LIME_FFmsy$MARE, Unex_LIME_FFmsy$MARE, Ovex_LIME_FFmsy$MARE))
Exp_FFmsy.err$LBRA <- cbind(c(Med_LBRA_FFmsy$MRE, Unex_LBRA_FFmsy$MRE, Ovex_LBRA_FFmsy$MRE), c(Med_LBRA_FFmsy$MARE, Unex_LBRA_FFmsy$MARE, Ovex_LBRA_FFmsy$MARE))

for(i in 1:4){
  colnames(Exp_FFmsy.err[[i]]) <- c("MRE", "MARE")
  rownames(Exp_FFmsy.err[[i]]) <- c("Targ", "Unex", "Ovex")
}

saveRDS(Exp_FFmsy.err, file = "files/Exp_FFmsy_err.rds")




Exp_FM.err <- list()
Exp_FM.err$TB <- cbind(c(Med_TB_FM$MRE, Unex_TB_FM$MRE, Ovex_TB_FM$MRE), c(Med_TB_FM$MARE, Unex_TB_FM$MARE, Ovex_TB_FM$MARE))
Exp_FM.err$LBSPR <- cbind(c(Med_LBSPR_FM$MRE, Unex_LBSPR_FM$MRE, Ovex_LBSPR_FM$MRE), c(Med_LBSPR_FM$MARE, Unex_LBSPR_FM$MARE, Ovex_LBSPR_FM$MARE))
Exp_FM.err$LIME <- cbind(c(Med_LIME_FM$MRE, Unex_LIME_FM$MRE, Ovex_LIME_FM$MRE), c(Med_LIME_FM$MARE, Unex_LIME_FM$MARE, Ovex_LIME_FM$MARE))
Exp_FM.err$LBRA <- cbind(c(Med_LBRA_FM$MRE, Unex_LBRA_FM$MRE, Ovex_LBRA_FM$MRE), c(Med_LBRA_FM$MARE, Unex_LBRA_FM$MARE, Ovex_LBRA_FM$MARE))

for(i in 1:4){
  colnames(Exp_FM.err[[i]]) <- c("MRE", "MARE")
  rownames(Exp_FM.err[[i]]) <- c("Targ", "Unex", "Ovex")
}

saveRDS(Exp_FM.err, file = "files/Exp_FM_err.rds")





Exp_Fmsy.err <- list()
Exp_Fmsy.err$TB <- cbind(c(Med_TB_Fmsy$MRE, Unex_TB_Fmsy$MRE, Ovex_TB_Fmsy$MRE), c(Med_TB_Fmsy$MARE, Unex_TB_Fmsy$MARE, Ovex_TB_Fmsy$MARE))
Exp_Fmsy.err$LBSPR <- cbind(c(Med_LBSPR_Fmsy$MRE, Unex_LBSPR_Fmsy$MRE, Ovex_LBSPR_Fmsy$MRE), c(Med_LBSPR_Fmsy$MARE, Unex_LBSPR_Fmsy$MARE, Ovex_LBSPR_Fmsy$MARE))
Exp_Fmsy.err$LIME <- cbind(c(Med_LIME_Fmsy$MRE, Unex_LIME_Fmsy$MRE, Ovex_LIME_Fmsy$MRE), c(Med_LIME_Fmsy$MARE, Unex_LIME_Fmsy$MARE, Ovex_LIME_Fmsy$MARE))
Exp_Fmsy.err$LBRA <- cbind(c(Med_LBRA_Fmsy$MRE, Unex_LBRA_Fmsy$MRE, Ovex_LBRA_Fmsy$MRE), c(Med_LBRA_Fmsy$MARE, Unex_LBRA_Fmsy$MARE, Ovex_LBRA_Fmsy$MARE))

for(i in 1:4){
  colnames(Exp_Fmsy.err[[i]]) <- c("MRE", "MARE")
  rownames(Exp_Fmsy.err[[i]]) <- c("Targ", "Unex", "Ovex")
}

saveRDS(Exp_Fmsy.err, file = "files/Exp_Fmsy_err.rds")




# Rec
Rec_SPR.err <- list()
Rec_SPR.err$TB <- cbind(c(Med_TB_SPR$MRE, Error_TB_SPR$MRE, AR_TB_SPR$MRE), c(Med_TB_SPR$MARE, Error_TB_SPR$MARE, AR_TB_SPR$MARE))
Rec_SPR.err$LBSPR <- cbind(c(Med_LBSPR_SPR$MRE, Error_LBSPR_SPR$MRE, AR_LBSPR_SPR$MRE), c(Med_LBSPR_SPR$MARE, Error_LBSPR_SPR$MARE, AR_LBSPR_SPR$MARE))
Rec_SPR.err$LIME <- cbind(c(Med_LIME_SPR$MRE, Error_LIME_SPR$MRE, AR_LIME_SPR$MRE), c(Med_LIME_SPR$MARE, Error_LIME_SPR$MARE, AR_LIME_SPR$MARE))
Rec_SPR.err$LBRA <- cbind(c(Med_LBRA_SPR$MRE, Error_LBRA_SPR$MRE, AR_LBRA_SPR$MRE), c(Med_LBRA_SPR$MARE, Error_LBRA_SPR$MARE, AR_LBRA_SPR$MARE))

for(i in 1:4){
  colnames(Rec_SPR.err[[i]]) <- c("MRE", "MARE")
  rownames(Rec_SPR.err[[i]]) <- c("None", "Error", "AR")
}

saveRDS(Rec_SPR.err, file = "files/Rec_SPR_err.rds")




Rec_SPRmsy.err <- list()
Rec_SPRmsy.err$TB <- cbind(c(Med_TB_SPRmsy$MRE, Error_TB_SPRmsy$MRE, AR_TB_SPRmsy$MRE), c(Med_TB_SPRmsy$MARE, Error_TB_SPRmsy$MARE, AR_TB_SPRmsy$MARE))
Rec_SPRmsy.err$LBSPR <- cbind(c(Med_LBSPR_SPRmsy$MRE, Error_LBSPR_SPRmsy$MRE, AR_LBSPR_SPRmsy$MRE), c(Med_LBSPR_SPRmsy$MARE, Error_LBSPR_SPRmsy$MARE, AR_LBSPR_SPRmsy$MARE))
Rec_SPRmsy.err$LIME <- cbind(c(Med_LIME_SPRmsy$MRE, Error_LIME_SPRmsy$MRE, AR_LIME_SPRmsy$MRE), c(Med_LIME_SPRmsy$MARE, Error_LIME_SPRmsy$MARE, AR_LIME_SPRmsy$MARE))
Rec_SPRmsy.err$LBRA <- cbind(c(Med_LBRA_SPRmsy$MRE, Error_LBRA_SPRmsy$MRE, AR_LBRA_SPRmsy$MRE), c(Med_LBRA_SPRmsy$MARE, Error_LBRA_SPRmsy$MARE, AR_LBRA_SPRmsy$MARE))

for(i in 1:4){
  colnames(Rec_SPRmsy.err[[i]]) <- c("MRE", "MARE")
  rownames(Rec_SPRmsy.err[[i]]) <- c("None", "Error", "AR")
}

saveRDS(Rec_SPRmsy.err, file = "files/Rec_SPRmsy_err.rds")




Rec_FFmsy.err <- list()
Rec_FFmsy.err$TB <- cbind(c(Med_TB_FFmsy$MRE, Error_TB_FFmsy$MRE, AR_TB_FFmsy$MRE), c(Med_TB_FFmsy$MARE, Error_TB_FFmsy$MARE, AR_TB_FFmsy$MARE))
Rec_FFmsy.err$LBSPR <- cbind(c(Med_LBSPR_FFmsy$MRE, Error_LBSPR_FFmsy$MRE, AR_LBSPR_FFmsy$MRE), c(Med_LBSPR_FFmsy$MARE, Error_LBSPR_FFmsy$MARE, AR_LBSPR_FFmsy$MARE))
Rec_FFmsy.err$LIME <- cbind(c(Med_LIME_FFmsy$MRE, Error_LIME_FFmsy$MRE, AR_LIME_FFmsy$MRE), c(Med_LIME_FFmsy$MARE, Error_LIME_FFmsy$MARE, AR_LIME_FFmsy$MARE))
Rec_FFmsy.err$LBRA <- cbind(c(Med_LBRA_FFmsy$MRE, Error_LBRA_FFmsy$MRE, AR_LBRA_FFmsy$MRE), c(Med_LBRA_FFmsy$MARE, Error_LBRA_FFmsy$MARE, AR_LBRA_FFmsy$MARE))

for(i in 1:4){
  colnames(Rec_FFmsy.err[[i]]) <- c("MRE", "MARE")
  rownames(Rec_FFmsy.err[[i]]) <- c("None", "Error", "AR")
}

saveRDS(Rec_FFmsy.err, file = "files/Rec_FFmsy_err.rds")




Rec_FM.err <- list()
Rec_FM.err$TB <- cbind(c(Med_TB_FM$MRE, Error_TB_FM$MRE, AR_TB_FM$MRE), c(Med_TB_FM$MARE, Error_TB_FM$MARE, AR_TB_FM$MARE))
Rec_FM.err$LBSPR <- cbind(c(Med_LBSPR_FM$MRE, Error_LBSPR_FM$MRE, AR_LBSPR_FM$MRE), c(Med_LBSPR_FM$MARE, Error_LBSPR_FM$MARE, AR_LBSPR_FM$MARE))
Rec_FM.err$LIME <- cbind(c(Med_LIME_FM$MRE, Error_LIME_FM$MRE, AR_LIME_FM$MRE), c(Med_LIME_FM$MARE, Error_LIME_FM$MARE, AR_LIME_FM$MARE))
Rec_FM.err$LBRA <- cbind(c(Med_LBRA_FM$MRE, Error_LBRA_FM$MRE, AR_LBRA_FM$MRE), c(Med_LBRA_FM$MARE, Error_LBRA_FM$MARE, AR_LBRA_FM$MARE))

for(i in 1:4){
  colnames(Rec_FM.err[[i]]) <- c("MRE", "MARE")
  rownames(Rec_FM.err[[i]]) <- c("None", "Error", "AR")
}

saveRDS(Rec_FM.err, file = "files/Rec_FM_err.rds")





Rec_Fmsy.err <- list()
Rec_Fmsy.err$TB <- cbind(c(Med_TB_Fmsy$MRE, Error_TB_Fmsy$MRE, AR_TB_Fmsy$MRE), c(Med_TB_Fmsy$MARE, Error_TB_Fmsy$MARE, AR_TB_Fmsy$MARE))
Rec_Fmsy.err$LBSPR <- cbind(c(Med_LBSPR_Fmsy$MRE, Error_LBSPR_Fmsy$MRE, AR_LBSPR_Fmsy$MRE), c(Med_LBSPR_Fmsy$MARE, Error_LBSPR_Fmsy$MARE, AR_LBSPR_Fmsy$MARE))
Rec_Fmsy.err$LIME <- cbind(c(Med_LIME_Fmsy$MRE, Error_LIME_Fmsy$MRE, AR_LIME_Fmsy$MRE), c(Med_LIME_Fmsy$MARE, Error_LIME_Fmsy$MARE, AR_LIME_Fmsy$MARE))
Rec_Fmsy.err$LBRA <- cbind(c(Med_LBRA_Fmsy$MRE, Error_LBRA_Fmsy$MRE, AR_LBRA_Fmsy$MRE), c(Med_LBRA_Fmsy$MARE, Error_LBRA_Fmsy$MARE, AR_LBRA_Fmsy$MARE))

for(i in 1:4){
  colnames(Rec_Fmsy.err[[i]]) <- c("MRE", "MARE")
  rownames(Rec_Fmsy.err[[i]]) <- c("None", "Error", "AR")
}

saveRDS(Rec_Fmsy.err, file = "files/Rec_Fmsy_err.rds")



# rm(list = setdiff(ls(), c("calc_error", "iters", "Med_LIME", "Med_LBRA", "Med_TB", "Med_LBSPR")))




