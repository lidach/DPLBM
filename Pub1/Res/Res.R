setwd("D:/Msc_Thesis/R/Res")



## Objects to keep
source("calc_error.R")



# Medlived Scenario (base model) ------------------------------------------
Med_LBRA <- readRDS("~/DPLBM/Pub1/Medlived/True/Medmodel_res_LBRA.rds")
Med_LBSPR <- readRDS("~/DPLBM/Pub1/Medlived/True/Medmodel_res_LBSPR.rds")
Med_LIME <- readRDS("~/DPLBM/Pub1/Medlived/True/Medmodel_res_LIME.rds")
Med_TB <- readRDS("~/DPLBM/Pub1/Medlived/True/Medmodel_res_TB.rds")
Med_list <- readRDS("~/DPLBM/Pub1/Medlived/Medlist.rds")
Med_list <- Med_list$ref
iters <- 300


## SPR 
SPRTrue <- sapply(Med_list, function(x) x$SPR)
LIMESPR <- sapply(Med_LIME, function(x) x$Report$SPR_t)
TBSPR <- sapply(Med_TB, function(x) x$currents$curr.SPR)

Med_LBRA_FM <- calc_error(SPRTrue, Med_LBRA$SPR, iters)
Med_LBSPR_FM <- calc_error(SPRTrue, Med_LBSPR$SPR, iters)
Med_LIME_FM <- calc_error(SPRTrue, LIMESPR, iters)
Med_TB_FM <- calc_error(SPRTrue, TBSPR, iters)


## F based reference points
## F 
FMTrue <- rep(0.128, iters)
LIMEFM <- sapply(Med_LIME, function(x) x$Report$F_t)
TBFM <- sapply(Med_TB, function(x) x$currents$curr.F)

Med_LBRA_FM <- calc_error(FMTrue, Med_LBRA$FM, iters)
Med_LIME_FM <- calc_error(FMTrue, LIMEFM, iters)
Med_TB_FM <- calc_error(FMTrue, TBFM, iters)

## Fmsy
FmsyTrue <- sapply(Med_list, function(x) x$Fdmsy)
LIMEFmsy <- sapply(Med_LIME, function(x) x$Derived$Fmsy)
TBFmsy <- sapply(Med_TB, function(x) x$df_Es$Fmax) # Fmax vs Fmsy??

Med_LBRA_Fmsy <- calc_error(FmsyTrue, Med_LBRA$Fmsy, iters)
Med_LIME_Fmsy <- calc_error(FmsyTrue, LIMEFmsy, iters)
Med_TB_Fmsy <- calc_error(FmsyTrue, TBFmsy, iters)

## F01
F01True <- sapply(Med_list, function(x) x$YPR$F01[25])
TBF01 <- sapply(Med_TB, function(x) x$df_Es$F01)

Med_TB_F01 <- calc_error(F01True, TBF01, iters)

## F/Fmsy
FFmsyTrue <- sapply(Med_list, function(x) x$states$F.Fmsy[25])
LIMEFFmsy <- sapply(Med_LIME, function(x) x$Derived$FFmsy)
TBFFmsy <- sapply(Med_TB, function(x) x$currents$curr.F/x$df_Es$Fmax)

Med_LBRA_FFmsy <- calc_error(FFmsyTrue, Med_LBRA$FFmsy, iters)
Med_LIME_FFmsy <- calc_error(FFmsyTrue, LIMEFFmsy, iters)
Med_TB_FFmsy <- calc_error(FFmsyTrue, TBFFmsy, iters)

## F/F01
FF01True <- sapply(Med_list, function(x) x$states$F.F01[25])
TBFF01 <- sapply(Med_TB, function(x) x$currents$curr.F/x$df_Es$F01)

Med_TB_FF01 <- calc_error(FF01True, TBFF01, iters)



## Res tab
Med_res <- matrix(nrow = 8, ncol = 6)
colnames(Med_res) <- c("SPR", "FM", "Fmsy", "FFmsy" ,"F01", "FF01")
rownames(Med_res) <- c("TB MRE", "LBSPR MRE", "LIME MRE", "LBRA MRE", 
                       "TB MARE", "LBSPR MARE", "LIME MARE", "LBRA MARE")

#
Med_res["TB MRE", "SPR"] <- Med_TB_SPR$MRE
Med_res["LBSPR MRE", "SPR"] <- Med_LBSPR_SPR$MRE
Med_res["LIME MRE", "SPR"] <- Med_LIME_SPR$MRE
Med_res["LBRA MRE", "SPR"] <- Med_LBRA_SPR$MRE

Med_res["TB MARE", "SPR"] <- Med_TB_SPR$MARE
Med_res["LBSPR MARE", "SPR"] <- Med_LBSPR_SPR$MARE
Med_res["LIME MARE", "SPR"] <- Med_LIME_SPR$MARE
Med_res["LBRA MARE", "SPR"] <- Med_LBRA_SPR$MARE

#
Med_res["TB MRE", "FM"] <- Med_TB_FM$MRE
Med_res["LIME MRE", "FM"] <- Med_LIME_FM$MRE
Med_res["LBRA MRE", "FM"] <- Med_LBRA_FM$MRE

Med_res["TB MARE", "FM"] <- Med_TB_FM$MARE
Med_res["LIME MARE", "FM"] <- Med_LIME_FM$MARE
Med_res["LBRA MARE", "FM"] <- Med_LBRA_FM$MARE

#
Med_res["TB MRE", "Fmsy"] <- Med_TB_Fmsy$MRE
Med_res["LIME MRE", "Fmsy"] <- Med_LIME_Fmsy$MRE
Med_res["LBRA MRE", "Fmsy"] <- Med_LBRA_Fmsy$MRE

Med_res["TB MARE", "Fmsy"] <- Med_TB_Fmsy$MARE
Med_res["LIME MARE", "Fmsy"] <- Med_LIME_Fmsy$MARE
Med_res["LBRA MARE", "Fmsy"] <- Med_LBRA_Fmsy$MARE

#
Med_res["TB MRE", "FFmsy"] <- Med_TB_FFmsy$MRE
Med_res["LIME MRE", "FFmsy"] <- Med_LIME_FFmsy$MRE
Med_res["LBRA MRE", "FFmsy"] <- Med_LBRA_FFmsy$MRE

Med_res["TB MARE", "FFmsy"] <- Med_TB_FFmsy$MARE
Med_res["LIME MARE", "FFmsy"] <- Med_LIME_FFmsy$MARE
Med_res["LBRA MARE", "FFmsy"] <- Med_LBRA_FFmsy$MARE

#
Med_res["TB MRE", "F01"] <- Med_TB_F01$MRE

Med_res["TB MARE", "F01"] <- Med_TB_F01$MARE

#
Med_res["TB MRE", "FF01"] <- Med_TB_FF01$MRE

Med_res["TB MARE", "FF01"] <- Med_TB_FF01$MARE


# 
# 
# par(mfrow = c(2,2), mar = c(1.5,3,1.5,3))
# 
# boxplot(Med_TB_SPR$rel_error, main = "TB SPR")
# abline(h = 0, lty = 2, col = "red")
# 
# boxplot(Med_LBSPR_SPR$rel_error, main = "LBSPR SPR")
# abline(h = 0, lty = 2, col = "red")
# 
# boxplot(Med_LIME_SPR$rel_error, main = "LIME SPR")
# abline(h = 0, lty = 2, col = "red")
# 
# boxplot(Med_LBRA_SPR$rel_error, main = "LBRA SPR")
# abline(h = 0, lty = 2, col = "red")
# 
# 
# #
# boxplot(Med_TB_FM$rel_error, main = "TB FM")
# abline(h = 0, lty = 2, col = "red")
# 
# boxplot(Med_LIME_FM$rel_error, main = "LIME FM")
# abline(h = 0, lty = 2, col = "red")
# 
# boxplot(Med_LBRA_FM$rel_error, main = "LBRA FM")
# abline(h = 0, lty = 2, col = "red")
# 
# # 
# boxplot(Med_TB_Fmsy$rel_error, main = "TB Fmsy")
# abline(h = 0, lty = 2, col = "red")
# 
# boxplot(Med_LIME_Fmsy$rel_error, main = "LIME Fmsy")
# abline(h = 0, lty = 2, col = "red")
# 
# boxplot(Med_LBRA_Fmsy$rel_error, main = "LBRA Fmsy")
# abline(h = 0, lty = 2, col = "red")
# 
# # 
# boxplot(Med_TB_FFmsy$rel_error, main = "TB FFmsy")
# abline(h = 0, lty = 2, col = "red")
# 
# boxplot(Med_LIME_FFmsy$rel_error, main = "LIME FFmsy")
# abline(h = 0, lty = 2, col = "red")
# 
# boxplot(Med_LBRA_FFmsy$rel_error, main = "LBRA FFmsy")
# abline(h = 0, lty = 2, col = "red")
# 
# 
# # 
# boxplot(Med_TB_FM$rel_error, main = "TB FM")
# abline(h = 0, lty = 2, col = "red")
# 
# boxplot(Med_TB_F01$rel_error, main = "TB F01")
# abline(h = 0, lty = 2, col = "red")
# 
# boxplot(Med_TB_FF01$rel_error, main = "TB FF01")
# abline(h = 0, lty = 2, col = "red")


# rm(list = setdiff(ls(), c("calc_error", "list")))