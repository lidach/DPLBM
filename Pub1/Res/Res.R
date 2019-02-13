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
LIMESPR <- sapply(Med_LIME, function(x) x$Report$SPR_t)
TBSPR <- sapply(Med_TB, function(x) x$currents$curr.SPR)

Med_LBRA_SPR <- calc_error(SPRTrue, Med_LBRA$SPR, iters)
Med_LBSPR_SPR <- calc_error(SPRTrue, Med_LBSPR$SPR, iters)
Med_LIME_SPR <- calc_error(SPRTrue, LIMESPR, iters)
Med_TB_SPR <- calc_error(SPRTrue, TBSPR, iters)

# SPRmsy = 0.183
SPRmsyT <- rep(0.183, iters)
# LIMESPR *****
# TBSPR ****

Med_LBRA_SPRmsy <- calc_error(SPRmsyT, Med_LBRA$SPRmsy, iters)
Med_LBSPR_SPRmsy <- calc_error(SPRmsyT, Med_LBSPR$SPRmsy, iters)
# Med_LIME_SPRmsy <- calc_error(SPRmsyT, LIMESPR, iters)
# Med_TB_SPRmsy <- calc_error(SPRmsyT, TBSPR, iters)



# F based reference points ####
# F = 0.128
FMTrue <- rep(0.128, iters)
LIMEFM <- sapply(Med_LIME, function(x) x$Report$F_t)
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
LIMEFFmsy <- sapply(Med_LIME, function(x) x$Derived$FFmsy)
TBFFmsy <- sapply(Med_TB, function(x) x$currents$curr.F/x$df_Es$Fmax)

Med_LBRA_FFmsy <- calc_error(FFmsyTrue, Med_LBRA$FFmsy, iters)
Med_LIME_FFmsy <- calc_error(FFmsyTrue, LIMEFFmsy, iters)
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
LIMESPR <- sapply(Short_LIME, function(x) x$Report$SPR_t)
TBSPR <- sapply(Short_TB, function(x) x$currents$curr.SPR)

Short_LBRA_SPR <- calc_error(SPRTrue, Short_LBRA$SPR, iters)
Short_LBSPR_SPR <- calc_error(SPRTrue, Short_LBSPR$SPR, iters)
Short_LIME_SPR <- calc_error(SPRTrue, LIMESPR, iters)
Short_TB_SPR <- calc_error(SPRTrue, TBSPR, iters)


# F = 0.45






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


# F = 0.075






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






# Error -------------------------------------------------------------------
Error_LBRA <- readRDS("D:/DPLBM/Pub1/Recruiterror/True/Errormodel_res_LBRA.rds")
Error_LBSPR <- readRDS("D:/DPLBM/Pub1/Recruiterror/True/Errormodel_res_LBSPR.rds")
Error_LIME <- readRDS("D:/DPLBM/Pub1/Recruiterror/True/Errormodel_res_LIME.rds")
Error_TB <- readRDS("D:/DPLBM/Pub1/Recruiterror/True/Errormodel_res_TB.rds")

## SPR ####
## 0.17171820
SPRTrue <- rep(0.17171820, iters)
LIMESPR <- sapply(Error_LIME, function(x) x$Report$SPR_t)
TBSPR <- sapply(Error_TB, function(x) x$currents$curr.SPR)

Error_LBRA_SPR <- calc_error(SPRTrue, Error_LBRA$SPR, iters)
Error_LBSPR_SPR <- calc_error(SPRTrue, Error_LBSPR$SPR, iters)
Error_LIME_SPR <- calc_error(SPRTrue, LIMESPR, iters)
Ovex_TB_SPR <- calc_error(SPRTrue, TBSPR, iters)







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




# MRE and MARE ------------------------------------------------------------

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


# rm(list = setdiff(ls(), c("calc_error", "iters", "Med_LIME", "Med_LBRA", "Med_TB", "Med_LBSPR")))




