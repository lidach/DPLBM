setwd("D:/DPLBM/Pub1/Res")

# Packages and Objects to save --------------------------------------------


if(!require(ggplot2)) install.packages("ggplot"); require(ggplot2)
if(!require(gridExtra)) install.packages("gridExtra"); require(gridExtra)



# LH comparison -----------------------------------------------------------
## SPR ####
LH_SPR <- readRDS("files/LH_SPR.rds")

Plot.LH_SPR <- ggplot(LH_SPR, aes(x = Method, y = Relative_error, fill = Method)) +
  geom_violin(trim = F) +
  labs(title = "Spawning Potential Ratio", y = "Relative Error") +
  stat_summary(fun.y = median, geom = "point", size = 2, color = "black") +
  geom_hline(yintercept = 0, color = 'darkgrey', size = 1) +
  facet_wrap(~ Scenario, ncol = 3, scale = "fixed") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 20, colour = "black")) +
  theme(axis.text.y = element_text(size = 20))
Plot.LH_SPR


LH_SPR.err <- readRDS("files/LH_SPR_err.rds")

tab <- data.frame(
  TB = round(LH_SPR.err$TB[1, c(1:2)],4),
  LBSPR = round(LH_SPR.err$LBSPR[1, c(1:2)],4),
  LIME = round(LH_SPR.err$LIME[1, c(1:2)],4),
  LBRA = round(LH_SPR.err$LBRA[1, c(1:2)],4),
  TB = round(LH_SPR.err$TB[2, c(1:2)],4),
  LBSPR = round(LH_SPR.err$LBSPR[2, c(1:2)],4),
  LIME = round(LH_SPR.err$LIME[2, c(1:2)],4),
  LBRA = round(LH_SPR.err$LBRA[2, c(1:2)],4),
  TB = round(LH_SPR.err$TB[3, c(1:2)],4),
  LBSPR = round(LH_SPR.err$LBSPR[3, c(1:2)],4),
  LIME = round(LH_SPR.err$LIME[3, c(1:2)],4),
  LBRA = round(LH_SPR.err$LBRA[3, c(1:2)],4)
)

colnames(tab) <- c(expression(TB[Medium]),
                   expression(LBSPR[Medium]),
                   expression(LIME[Medium]),
                   expression(LBRA[Medium]),
                   expression(TB[Short]),
                   expression(LBSPR[Short]),
                   expression(LIME[Short]),
                   expression(LBRA[Short]),
                   expression(TB[Long]),
                   expression(LBSPR[Long]),
                   expression(LIME[Long]),
                   expression(LBRA[Long]))
rownames(tab) <- c("MRE", "MARE")


tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(tab, theme=tt)


Plot.LH_SPR1 <- grid.arrange(Plot.LH_SPR, tbl,
                                nrow = 2,
                                as.table = TRUE,
                                heights = c(3,1))

# ggsave("Plot_SPRTRUE_LH.png", plot = Plot_SPRTrue_LH,
#        scale = 1, width = 14, height = 7, units = "in")




## FFmsy ####
LH_FFmsy <- readRDS("files/LH_FFmsy.rds")

Plot.LH_FFmsy <- ggplot(LH_FFmsy, aes(x = Method, y = Relative_error, fill = Method)) +
  geom_violin(trim = F) +
  labs(title = "F/Fmsy", y = "Relative Error") +
  stat_summary(fun.y = median, geom = "point", size = 2, color = "black") +
  geom_hline(yintercept = 0, color = 'darkgrey', size = 1) +
  facet_wrap(~ Scenario, ncol = 3, scale = "fixed") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 20, colour = "black")) +
  theme(axis.text.y = element_text(size = 20))
Plot.LH_FFmsy


LH_FFmsy.err <- readRDS("files/LH_FFmsy_err.rds")

tab <- data.frame(
  TB = round(LH_FFmsy.err$TB[1, c(1:2)],4),
  LBSPR = round(LH_FFmsy.err$LBSPR[1, c(1:2)],4),
  LIME = round(LH_FFmsy.err$LIME[1, c(1:2)],4),
  LBRA = round(LH_FFmsy.err$LBRA[1, c(1:2)],4),
  TB = round(LH_FFmsy.err$TB[2, c(1:2)],4),
  LBSPR = round(LH_FFmsy.err$LBSPR[2, c(1:2)],4),
  LIME = round(LH_FFmsy.err$LIME[2, c(1:2)],4),
  LBRA = round(LH_FFmsy.err$LBRA[2, c(1:2)],4),
  TB = round(LH_FFmsy.err$TB[3, c(1:2)],4),
  LBSPR = round(LH_FFmsy.err$LBSPR[3, c(1:2)],4),
  LIME = round(LH_FFmsy.err$LIME[3, c(1:2)],4),
  LBRA = round(LH_FFmsy.err$LBRA[3, c(1:2)],4)
)

colnames(tab) <- c(expression(TB[Medium]),
                   expression(LBSPR[Medium]),
                   expression(LIME[Medium]),
                   expression(LBRA[Medium]),
                   expression(TB[Short]),
                   expression(LBSPR[Short]),
                   expression(LIME[Short]),
                   expression(LBRA[Short]),
                   expression(TB[Long]),
                   expression(LBSPR[Long]),
                   expression(LIME[Long]),
                   expression(LBRA[Long]))
rownames(tab) <- c("MRE", "MARE")


tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(tab, theme=tt)


Plot.LH_FFmsy1 <- grid.arrange(Plot.LH_FFmsy, tbl,
                             nrow = 2,
                             as.table = TRUE,
                             heights = c(3,1))

# ggsave("Plot_FFmsyTRUE_LH.png", plot = Plot_FFmsyTrue_LH,
#        scale = 1, width = 14, height = 7, units = "in")




## FM ####
LH_FM <- readRDS("files/LH_FM.rds")

Plot.LH_FM <- ggplot(LH_FM, aes(x = Method, y = Relative_error, fill = Method)) +
  geom_violin(trim = F) +
  labs(title = "Fishing mortality", y = "Relative Error") +
  stat_summary(fun.y = median, geom = "point", size = 2, color = "black") +
  geom_hline(yintercept = 0, color = 'darkgrey', size = 1) +
  facet_wrap(~ Scenario, ncol = 3, scale = "fixed") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 20, colour = "black")) +
  theme(axis.text.y = element_text(size = 20))
Plot.LH_FM


LH_FM.err <- readRDS("files/LH_FM_err.rds")

tab <- data.frame(
  TB = round(LH_FM.err$TB[1, c(1:2)],4),
  LBSPR = round(LH_FM.err$LBSPR[1, c(1:2)],4),
  LIME = round(LH_FM.err$LIME[1, c(1:2)],4),
  LBRA = round(LH_FM.err$LBRA[1, c(1:2)],4),
  TB = round(LH_FM.err$TB[2, c(1:2)],4),
  LBSPR = round(LH_FM.err$LBSPR[2, c(1:2)],4),
  LIME = round(LH_FM.err$LIME[2, c(1:2)],4),
  LBRA = round(LH_FM.err$LBRA[2, c(1:2)],4),
  TB = round(LH_FM.err$TB[3, c(1:2)],4),
  LBSPR = round(LH_FM.err$LBSPR[3, c(1:2)],4),
  LIME = round(LH_FM.err$LIME[3, c(1:2)],4),
  LBRA = round(LH_FM.err$LBRA[3, c(1:2)],4)
)

colnames(tab) <- c(expression(TB[Medium]),
                   expression(LBSPR[Medium]),
                   expression(LIME[Medium]),
                   expression(LBRA[Medium]),
                   expression(TB[Short]),
                   expression(LBSPR[Short]),
                   expression(LIME[Short]),
                   expression(LBRA[Short]),
                   expression(TB[Long]),
                   expression(LBSPR[Long]),
                   expression(LIME[Long]),
                   expression(LBRA[Long]))
rownames(tab) <- c("MRE", "MARE")


tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(tab, theme=tt)


Plot.LH_FM1 <- grid.arrange(Plot.LH_FM, tbl,
                             nrow = 2,
                             as.table = TRUE,
                             heights = c(3,1))

# ggsave("Plot_FMTRUE_LH.png", plot = Plot_FMTrue_LH,
#        scale = 1, width = 14, height = 7, units = "in")






## Fmsy ####
LH_Fmsy <- readRDS("files/LH_Fmsy.rds")

Plot.LH_Fmsy <- ggplot(LH_Fmsy, aes(x = Method, y = Relative_error, fill = Method)) +
  geom_violin(trim = F) +
  labs(title = "Fmsy", y = "Relative Error") +
  stat_summary(fun.y = median, geom = "point", size = 2, color = "black") +
  geom_hline(yintercept = 0, color = 'darkgrey', size = 1) +
  facet_wrap(~ Scenario, ncol = 3, scale = "fixed") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 20, colour = "black")) +
  theme(axis.text.y = element_text(size = 20))
Plot.LH_Fmsy


LH_Fmsy.err <- readRDS("files/LH_Fmsy_err.rds")

tab <- data.frame(
  TB = round(LH_Fmsy.err$TB[1, c(1:2)],4),
  LBSPR = round(LH_Fmsy.err$LBSPR[1, c(1:2)],4),
  LIME = round(LH_Fmsy.err$LIME[1, c(1:2)],4),
  LBRA = round(LH_Fmsy.err$LBRA[1, c(1:2)],4),
  TB = round(LH_Fmsy.err$TB[2, c(1:2)],4),
  LBSPR = round(LH_Fmsy.err$LBSPR[2, c(1:2)],4),
  LIME = round(LH_Fmsy.err$LIME[2, c(1:2)],4),
  LBRA = round(LH_Fmsy.err$LBRA[2, c(1:2)],4),
  TB = round(LH_Fmsy.err$TB[3, c(1:2)],4),
  LBSPR = round(LH_Fmsy.err$LBSPR[3, c(1:2)],4),
  LIME = round(LH_Fmsy.err$LIME[3, c(1:2)],4),
  LBRA = round(LH_Fmsy.err$LBRA[3, c(1:2)],4)
)

colnames(tab) <- c(expression(TB[Medium]),
                   expression(LBSPR[Medium]),
                   expression(LIME[Medium]),
                   expression(LBRA[Medium]),
                   expression(TB[Short]),
                   expression(LBSPR[Short]),
                   expression(LIME[Short]),
                   expression(LBRA[Short]),
                   expression(TB[Long]),
                   expression(LBSPR[Long]),
                   expression(LIME[Long]),
                   expression(LBRA[Long]))
rownames(tab) <- c("MRE", "MARE")


tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(tab, theme=tt)


Plot.LH_Fmsy1 <- grid.arrange(Plot.LH_Fmsy, tbl,
                            nrow = 2,
                            as.table = TRUE,
                            heights = c(3,1))

# ggsave("Plot_FmsyTRUE_LH.png", plot = Plot_FmsyTrue_LH,
#        scale = 1, width = 14, height = 7, units = "in")





## SPRmsy ####
LH_SPRmsy <- readRDS("files/LH_SPRmsy.rds")

Plot.LH_SPRmsy <- ggplot(LH_SPRmsy, aes(x = Method, y = Relative_error, fill = Method)) +
  geom_violin(trim = F) +
  labs(title = "SPRmsy", y = "Relative Error") +
  stat_summary(fun.y = median, geom = "point", size = 2, color = "black") +
  geom_hline(yintercept = 0, color = 'darkgrey', size = 1) +
  facet_wrap(~ Scenario, ncol = 3, scale = "fixed") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 20, colour = "black")) +
  theme(axis.text.y = element_text(size = 20))
Plot.LH_SPRmsy


LH_SPRmsy.err <- readRDS("files/LH_SPRmsy_err.rds")

tab <- data.frame(
  TB = round(LH_SPRmsy.err$TB[1, c(1:2)],4),
  LBSPR = round(LH_SPRmsy.err$LBSPR[1, c(1:2)],4),
  LIME = round(LH_SPRmsy.err$LIME[1, c(1:2)],4),
  LBRA = round(LH_SPRmsy.err$LBRA[1, c(1:2)],4),
  TB = round(LH_SPRmsy.err$TB[2, c(1:2)],4),
  LBSPR = round(LH_SPRmsy.err$LBSPR[2, c(1:2)],4),
  LIME = round(LH_SPRmsy.err$LIME[2, c(1:2)],4),
  LBRA = round(LH_SPRmsy.err$LBRA[2, c(1:2)],4),
  TB = round(LH_SPRmsy.err$TB[3, c(1:2)],4),
  LBSPR = round(LH_SPRmsy.err$LBSPR[3, c(1:2)],4),
  LIME = round(LH_SPRmsy.err$LIME[3, c(1:2)],4),
  LBRA = round(LH_SPRmsy.err$LBRA[3, c(1:2)],4)
)

colnames(tab) <- c(expression(TB[Medium]),
                   expression(LBSPR[Medium]),
                   expression(LIME[Medium]),
                   expression(LBRA[Medium]),
                   expression(TB[Short]),
                   expression(LBSPR[Short]),
                   expression(LIME[Short]),
                   expression(LBRA[Short]),
                   expression(TB[Long]),
                   expression(LBSPR[Long]),
                   expression(LIME[Long]),
                   expression(LBRA[Long]))
rownames(tab) <- c("MRE", "MARE")


tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(tab, theme=tt)


Plot.LH_SPRmsy1 <- grid.arrange(Plot.LH_SPRmsy, tbl,
                              nrow = 2,
                              as.table = TRUE,
                              heights = c(3,1))

# ggsave("Plot_SPRmsyTRUE_LH.png", plot = Plot_SPRmsyTrue_LH,
#        scale = 1, width = 14, height = 7, units = "in")





# Exp comparison ----------------------------------------------------------


Exp_SPR <- readRDS("files/Exp_SPR.rds")

Plot.Exp_SPR <- ggplot(Exp_SPR, aes(x = Method, y = Relative_error, fill = Method)) +
  geom_violin(trim = F) +
  labs(title = "Spawning Potential Ratio", y = "Relative Error") +
  stat_summary(fun.y = median, geom = "point", size = 2, color = "black") +
  geom_hline(yintercept = 0, color = 'darkgrey', size = 1) +
  facet_wrap(~ Scenario, ncol = 3, scale = "fixed") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 20, colour = "black")) +
  theme(axis.text.y = element_text(size = 20))
Plot.Exp_SPR


Exp_SPR.err <- readRDS("files/Exp_SPR_err.rds")

tab <- data.frame(
  TB = round(Exp_SPR.err$TB[1, c(1:2)],4),
  LBSPR = round(Exp_SPR.err$LBSPR[1, c(1:2)],4),
  LIME = round(Exp_SPR.err$LIME[1, c(1:2)],4),
  LBRA = round(Exp_SPR.err$LBRA[1, c(1:2)],4),
  TB = round(Exp_SPR.err$TB[2, c(1:2)],4),
  LBSPR = round(Exp_SPR.err$LBSPR[2, c(1:2)],4),
  LIME = round(Exp_SPR.err$LIME[2, c(1:2)],4),
  LBRA = round(Exp_SPR.err$LBRA[2, c(1:2)],4),
  TB = round(Exp_SPR.err$TB[3, c(1:2)],4),
  LBSPR = round(Exp_SPR.err$LBSPR[3, c(1:2)],4),
  LIME = round(Exp_SPR.err$LIME[3, c(1:2)],4),
  LBRA = round(Exp_SPR.err$LBRA[3, c(1:2)],4)
)

colnames(tab) <- c(expression(TB[Target]),
                   expression(LBSPR[Target]),
                   expression(LIME[Target]),
                   expression(LBRA[Target]),
                   expression(TB[Unex]),
                   expression(LBSPR[Unex]),
                   expression(LIME[Unex]),
                   expression(LBRA[Unex]),
                   expression(TB[Ovex]),
                   expression(LBSPR[Ovex]),
                   expression(LIME[Ovex]),
                   expression(LBRA[Ovex]))
rownames(tab) <- c("MRE", "MARE")


tt <- ttheme_default(coExpead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(tab, theme=tt)


Plot.Exp_SPR1 <- grid.arrange(Plot.Exp_SPR, tbl,
                             nrow = 2,
                             as.table = TRUE,
                             heights = c(3,1))

# ggsave("Plot_SPRTRUE_Exp.png", plot = Plot_SPRTrue_Exp,
#        scale = 1, width = 14, height = 7, units = "in")







## SPR ####
Exp_SPR <- readRDS("files/Exp_SPR.rds")

Plot.Exp_SPR <- ggplot(Exp_SPR, aes(x = Method, y = Relative_error, fill = Method)) +
  geom_violin(trim = F) +
  labs(title = "Spawning Potential Ratio", y = "Relative Error") +
  stat_summary(fun.y = median, geom = "point", size = 2, color = "black") +
  geom_hline(yintercept = 0, color = 'darkgrey', size = 1) +
  facet_wrap(~ Scenario, ncol = 3, scale = "fixed") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 20, colour = "black")) +
  theme(axis.text.y = element_text(size = 20))
Plot.Exp_SPR


Exp_SPR.err <- readRDS("files/Exp_SPR_err.rds")

tab <- data.frame(
  TB = round(Exp_SPR.err$TB[1, c(1:2)],4),
  LBSPR = round(Exp_SPR.err$LBSPR[1, c(1:2)],4),
  LIME = round(Exp_SPR.err$LIME[1, c(1:2)],4),
  LBRA = round(Exp_SPR.err$LBRA[1, c(1:2)],4),
  TB = round(Exp_SPR.err$TB[2, c(1:2)],4),
  LBSPR = round(Exp_SPR.err$LBSPR[2, c(1:2)],4),
  LIME = round(Exp_SPR.err$LIME[2, c(1:2)],4),
  LBRA = round(Exp_SPR.err$LBRA[2, c(1:2)],4),
  TB = round(Exp_SPR.err$TB[3, c(1:2)],4),
  LBSPR = round(Exp_SPR.err$LBSPR[3, c(1:2)],4),
  LIME = round(Exp_SPR.err$LIME[3, c(1:2)],4),
  LBRA = round(Exp_SPR.err$LBRA[3, c(1:2)],4)
)

colnames(tab) <- c(expression(TB[Target[),
                   expression(LBSPR[Target[),
                   expression(LIME[Target[),
                   expression(LBRA[Target[),
                   expression(TB[Unex]),
                   expression(LBSPR[Unex]),
                   expression(LIME[Unex]),
                   expression(LBRA[Unex]),
                   expression(TB[Ovex]),
                   expression(LBSPR[Ovex]),
                   expression(LIME[Ovex]),
                   expression(LBRA[Ovex]))
rownames(tab) <- c("MRE", "MARE")


tt <- ttheme_default(coExpead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(tab, theme=tt)


Plot.Exp_SPR1 <- grid.arrange(Plot.Exp_SPR, tbl,
                             nrow = 2,
                             as.table = TRUE,
                             heights = c(3,1))

# ggsave("Plot_SPRTRUE_Exp.png", plot = Plot_SPRTrue_Exp,
#        scale = 1, width = 14, height = 7, units = "in")




## FFmsy ####
Exp_FFmsy <- readRDS("files/Exp_FFmsy.rds")

Plot.Exp_FFmsy <- ggplot(Exp_FFmsy, aes(x = Method, y = Relative_error, fill = Method)) +
  geom_violin(trim = F) +
  labs(title = "F/Fmsy", y = "Relative Error") +
  stat_summary(fun.y = median, geom = "point", size = 2, color = "black") +
  geom_hline(yintercept = 0, color = 'darkgrey', size = 1) +
  facet_wrap(~ Scenario, ncol = 3, scale = "fixed") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 20, colour = "black")) +
  theme(axis.text.y = element_text(size = 20))
Plot.Exp_FFmsy


Exp_FFmsy.err <- readRDS("files/Exp_FFmsy_err.rds")

tab <- data.frame(
  TB = round(Exp_FFmsy.err$TB[1, c(1:2)],4),
  LBSPR = round(Exp_FFmsy.err$LBSPR[1, c(1:2)],4),
  LIME = round(Exp_FFmsy.err$LIME[1, c(1:2)],4),
  LBRA = round(Exp_FFmsy.err$LBRA[1, c(1:2)],4),
  TB = round(Exp_FFmsy.err$TB[2, c(1:2)],4),
  LBSPR = round(Exp_FFmsy.err$LBSPR[2, c(1:2)],4),
  LIME = round(Exp_FFmsy.err$LIME[2, c(1:2)],4),
  LBRA = round(Exp_FFmsy.err$LBRA[2, c(1:2)],4),
  TB = round(Exp_FFmsy.err$TB[3, c(1:2)],4),
  LBSPR = round(Exp_FFmsy.err$LBSPR[3, c(1:2)],4),
  LIME = round(Exp_FFmsy.err$LIME[3, c(1:2)],4),
  LBRA = round(Exp_FFmsy.err$LBRA[3, c(1:2)],4)
)

colnames(tab) <- c(expression(TB[Target[),
                   expression(LBSPR[Target[),
                   expression(LIME[Target[),
                   expression(LBRA[Target[),
                   expression(TB[Unex]),
                   expression(LBSPR[Unex]),
                   expression(LIME[Unex]),
                   expression(LBRA[Unex]),
                   expression(TB[Ovex]),
                   expression(LBSPR[Ovex]),
                   expression(LIME[Ovex]),
                   expression(LBRA[Ovex]))
rownames(tab) <- c("MRE", "MARE")


tt <- ttheme_default(coExpead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(tab, theme=tt)


Plot.Exp_FFmsy1 <- grid.arrange(Plot.Exp_FFmsy, tbl,
                               nrow = 2,
                               as.table = TRUE,
                               heights = c(3,1))

# ggsave("Plot_FFmsyTRUE_Exp.png", plot = Plot_FFmsyTrue_Exp,
#        scale = 1, width = 14, height = 7, units = "in")




## FM ####
Exp_FM <- readRDS("files/Exp_FM.rds")

Plot.Exp_FM <- ggplot(Exp_FM, aes(x = Method, y = Relative_error, fill = Method)) +
  geom_violin(trim = F) +
  labs(title = "Fishing mortality", y = "Relative Error") +
  stat_summary(fun.y = median, geom = "point", size = 2, color = "black") +
  geom_hline(yintercept = 0, color = 'darkgrey', size = 1) +
  facet_wrap(~ Scenario, ncol = 3, scale = "fixed") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 20, colour = "black")) +
  theme(axis.text.y = element_text(size = 20))
Plot.Exp_FM


Exp_FM.err <- readRDS("files/Exp_FM_err.rds")

tab <- data.frame(
  TB = round(Exp_FM.err$TB[1, c(1:2)],4),
  LBSPR = round(Exp_FM.err$LBSPR[1, c(1:2)],4),
  LIME = round(Exp_FM.err$LIME[1, c(1:2)],4),
  LBRA = round(Exp_FM.err$LBRA[1, c(1:2)],4),
  TB = round(Exp_FM.err$TB[2, c(1:2)],4),
  LBSPR = round(Exp_FM.err$LBSPR[2, c(1:2)],4),
  LIME = round(Exp_FM.err$LIME[2, c(1:2)],4),
  LBRA = round(Exp_FM.err$LBRA[2, c(1:2)],4),
  TB = round(Exp_FM.err$TB[3, c(1:2)],4),
  LBSPR = round(Exp_FM.err$LBSPR[3, c(1:2)],4),
  LIME = round(Exp_FM.err$LIME[3, c(1:2)],4),
  LBRA = round(Exp_FM.err$LBRA[3, c(1:2)],4)
)

colnames(tab) <- c(expression(TB[Target[),
                   expression(LBSPR[Target[),
                   expression(LIME[Target[),
                   expression(LBRA[Target[),
                   expression(TB[Unex]),
                   expression(LBSPR[Unex]),
                   expression(LIME[Unex]),
                   expression(LBRA[Unex]),
                   expression(TB[Ovex]),
                   expression(LBSPR[Ovex]),
                   expression(LIME[Ovex]),
                   expression(LBRA[Ovex]))
rownames(tab) <- c("MRE", "MARE")


tt <- ttheme_default(coExpead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(tab, theme=tt)


Plot.Exp_FM1 <- grid.arrange(Plot.Exp_FM, tbl,
                            nrow = 2,
                            as.table = TRUE,
                            heights = c(3,1))

# ggsave("Plot_FMTRUE_Exp.png", plot = Plot_FMTrue_Exp,
#        scale = 1, width = 14, height = 7, units = "in")






## Fmsy ####
Exp_Fmsy <- readRDS("files/Exp_Fmsy.rds")

Plot.Exp_Fmsy <- ggplot(Exp_Fmsy, aes(x = Method, y = Relative_error, fill = Method)) +
  geom_violin(trim = F) +
  labs(title = "Fmsy", y = "Relative Error") +
  stat_summary(fun.y = median, geom = "point", size = 2, color = "black") +
  geom_hline(yintercept = 0, color = 'darkgrey', size = 1) +
  facet_wrap(~ Scenario, ncol = 3, scale = "fixed") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 20, colour = "black")) +
  theme(axis.text.y = element_text(size = 20))
Plot.Exp_Fmsy


Exp_Fmsy.err <- readRDS("files/Exp_Fmsy_err.rds")

tab <- data.frame(
  TB = round(Exp_Fmsy.err$TB[1, c(1:2)],4),
  LBSPR = round(Exp_Fmsy.err$LBSPR[1, c(1:2)],4),
  LIME = round(Exp_Fmsy.err$LIME[1, c(1:2)],4),
  LBRA = round(Exp_Fmsy.err$LBRA[1, c(1:2)],4),
  TB = round(Exp_Fmsy.err$TB[2, c(1:2)],4),
  LBSPR = round(Exp_Fmsy.err$LBSPR[2, c(1:2)],4),
  LIME = round(Exp_Fmsy.err$LIME[2, c(1:2)],4),
  LBRA = round(Exp_Fmsy.err$LBRA[2, c(1:2)],4),
  TB = round(Exp_Fmsy.err$TB[3, c(1:2)],4),
  LBSPR = round(Exp_Fmsy.err$LBSPR[3, c(1:2)],4),
  LIME = round(Exp_Fmsy.err$LIME[3, c(1:2)],4),
  LBRA = round(Exp_Fmsy.err$LBRA[3, c(1:2)],4)
)

colnames(tab) <- c(expression(TB[Target[),
                   expression(LBSPR[Target[),
                   expression(LIME[Target[),
                   expression(LBRA[Target[),
                   expression(TB[Unex]),
                   expression(LBSPR[Unex]),
                   expression(LIME[Unex]),
                   expression(LBRA[Unex]),
                   expression(TB[Ovex]),
                   expression(LBSPR[Ovex]),
                   expression(LIME[Ovex]),
                   expression(LBRA[Ovex]))
rownames(tab) <- c("MRE", "MARE")


tt <- ttheme_default(coExpead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(tab, theme=tt)


Plot.Exp_Fmsy1 <- grid.arrange(Plot.Exp_Fmsy, tbl,
                              nrow = 2,
                              as.table = TRUE,
                              heights = c(3,1))

# ggsave("Plot_FmsyTRUE_Exp.png", plot = Plot_FmsyTrue_Exp,
#        scale = 1, width = 14, height = 7, units = "in")





## SPRmsy ####
Exp_SPRmsy <- readRDS("files/Exp_SPRmsy.rds")

Plot.Exp_SPRmsy <- ggplot(Exp_SPRmsy, aes(x = Method, y = Relative_error, fill = Method)) +
  geom_violin(trim = F) +
  labs(title = "SPRmsy", y = "Relative Error") +
  stat_summary(fun.y = median, geom = "point", size = 2, color = "black") +
  geom_hline(yintercept = 0, color = 'darkgrey', size = 1) +
  facet_wrap(~ Scenario, ncol = 3, scale = "fixed") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 20, colour = "black")) +
  theme(axis.text.y = element_text(size = 20))
Plot.Exp_SPRmsy


Exp_SPRmsy.err <- readRDS("files/Exp_SPRmsy_err.rds")

tab <- data.frame(
  TB = round(Exp_SPRmsy.err$TB[1, c(1:2)],4),
  LBSPR = round(Exp_SPRmsy.err$LBSPR[1, c(1:2)],4),
  LIME = round(Exp_SPRmsy.err$LIME[1, c(1:2)],4),
  LBRA = round(Exp_SPRmsy.err$LBRA[1, c(1:2)],4),
  TB = round(Exp_SPRmsy.err$TB[2, c(1:2)],4),
  LBSPR = round(Exp_SPRmsy.err$LBSPR[2, c(1:2)],4),
  LIME = round(Exp_SPRmsy.err$LIME[2, c(1:2)],4),
  LBRA = round(Exp_SPRmsy.err$LBRA[2, c(1:2)],4),
  TB = round(Exp_SPRmsy.err$TB[3, c(1:2)],4),
  LBSPR = round(Exp_SPRmsy.err$LBSPR[3, c(1:2)],4),
  LIME = round(Exp_SPRmsy.err$LIME[3, c(1:2)],4),
  LBRA = round(Exp_SPRmsy.err$LBRA[3, c(1:2)],4)
)

colnames(tab) <- c(expression(TB[Target[),
                   expression(LBSPR[Target[),
                   expression(LIME[Target[),
                   expression(LBRA[Target[),
                   expression(TB[Unex]),
                   expression(LBSPR[Unex]),
                   expression(LIME[Unex]),
                   expression(LBRA[Unex]),
                   expression(TB[Ovex]),
                   expression(LBSPR[Ovex]),
                   expression(LIME[Ovex]),
                   expression(LBRA[Ovex]))
rownames(tab) <- c("MRE", "MARE")


tt <- ttheme_default(coExpead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(tab, theme=tt)


Plot.Exp_SPRmsy1 <- grid.arrange(Plot.Exp_SPRmsy, tbl,
                               nrow = 2,
                               as.table = TRUE,
                               heights = c(3,1))

# ggsave("Plot_SPRmsyTRUE_Exp.png", plot = Plot_SPRmsyTrue_Exp,
#        scale = 1, width = 14, height = 7, units = "in")





# Rec comparison ----------------------------------------------------------



## SPR ####
Rec_SPR <- readRDS("files/Rec_SPR.rds")

Plot.Rec_SPR <- ggplot(Rec_SPR, aes(x = Method, y = Relative_error, fill = Method)) +
  geom_violin(trim = F) +
  labs(title = "Spawning Potential Ratio", y = "Relative Error") +
  stat_summary(fun.y = median, geom = "point", size = 2, color = "black") +
  geom_hline(yintercept = 0, color = 'darkgrey', size = 1) +
  facet_wrap(~ Scenario, ncol = 3, scale = "fixed") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 20, colour = "black")) +
  theme(axis.text.y = element_text(size = 20))
Plot.Rec_SPR


Rec_SPR.err <- readRDS("files/Rec_SPR_err.rds")

tab <- data.frame(
  TB = round(Rec_SPR.err$TB[1, c(1:2)],4),
  LBSPR = round(Rec_SPR.err$LBSPR[1, c(1:2)],4),
  LIME = round(Rec_SPR.err$LIME[1, c(1:2)],4),
  LBRA = round(Rec_SPR.err$LBRA[1, c(1:2)],4),
  TB = round(Rec_SPR.err$TB[2, c(1:2)],4),
  LBSPR = round(Rec_SPR.err$LBSPR[2, c(1:2)],4),
  LIME = round(Rec_SPR.err$LIME[2, c(1:2)],4),
  LBRA = round(Rec_SPR.err$LBRA[2, c(1:2)],4),
  TB = round(Rec_SPR.err$TB[3, c(1:2)],4),
  LBSPR = round(Rec_SPR.err$LBSPR[3, c(1:2)],4),
  LIME = round(Rec_SPR.err$LIME[3, c(1:2)],4),
  LBRA = round(Rec_SPR.err$LBRA[3, c(1:2)],4)
)

colnames(tab) <- c(expression(TB[None]),
                   expression(LBSPR[None]),
                   expression(LIME[None]),
                   expression(LBRA[None]),
                   expression(TB[Error]),
                   expression(LBSPR[Error]),
                   expression(LIME[Error]),
                   expression(LBRA[Error]),
                   expression(TB[AR]),
                   expression(LBSPR[AR]),
                   expression(LIME[AR]),
                   expression(LBRA[AR]))
rownames(tab) <- c("MRE", "MARE")


tt <- ttheme_default(coRecead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(tab, theme=tt)


Plot.Rec_SPR1 <- grid.arrange(Plot.Rec_SPR, tbl,
                             nrow = 2,
                             as.table = TRUE,
                             heights = c(3,1))

# ggsave("Plot_SPRTRUE_Rec.png", plot = Plot_SPRTrue_Rec,
#        scale = 1, width = 14, height = 7, units = "in")




## FFmsy ####
Rec_FFmsy <- readRDS("files/Rec_FFmsy.rds")

Plot.Rec_FFmsy <- ggplot(Rec_FFmsy, aes(x = Method, y = Relative_error, fill = Method)) +
  geom_violin(trim = F) +
  labs(title = "F/Fmsy", y = "Relative Error") +
  stat_summary(fun.y = median, geom = "point", size = 2, color = "black") +
  geom_hline(yintercept = 0, color = 'darkgrey', size = 1) +
  facet_wrap(~ Scenario, ncol = 3, scale = "fixed") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 20, colour = "black")) +
  theme(axis.text.y = element_text(size = 20))
Plot.Rec_FFmsy


Rec_FFmsy.err <- readRDS("files/Rec_FFmsy_err.rds")

tab <- data.frame(
  TB = round(Rec_FFmsy.err$TB[1, c(1:2)],4),
  LBSPR = round(Rec_FFmsy.err$LBSPR[1, c(1:2)],4),
  LIME = round(Rec_FFmsy.err$LIME[1, c(1:2)],4),
  LBRA = round(Rec_FFmsy.err$LBRA[1, c(1:2)],4),
  TB = round(Rec_FFmsy.err$TB[2, c(1:2)],4),
  LBSPR = round(Rec_FFmsy.err$LBSPR[2, c(1:2)],4),
  LIME = round(Rec_FFmsy.err$LIME[2, c(1:2)],4),
  LBRA = round(Rec_FFmsy.err$LBRA[2, c(1:2)],4),
  TB = round(Rec_FFmsy.err$TB[3, c(1:2)],4),
  LBSPR = round(Rec_FFmsy.err$LBSPR[3, c(1:2)],4),
  LIME = round(Rec_FFmsy.err$LIME[3, c(1:2)],4),
  LBRA = round(Rec_FFmsy.err$LBRA[3, c(1:2)],4)
)

colnames(tab) <- c(expression(TB[None]),
                   expression(LBSPR[None]),
                   expression(LIME[None]),
                   expression(LBRA[None]),
                   expression(TB[Error]),
                   expression(LBSPR[Error]),
                   expression(LIME[Error]),
                   expression(LBRA[Error]),
                   expression(TB[AR]),
                   expression(LBSPR[AR]),
                   expression(LIME[AR]),
                   expression(LBRA[AR]))
rownames(tab) <- c("MRE", "MARE")


tt <- ttheme_default(coRecead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(tab, theme=tt)


Plot.Rec_FFmsy1 <- grid.arrange(Plot.Rec_FFmsy, tbl,
                               nrow = 2,
                               as.table = TRUE,
                               heights = c(3,1))

# ggsave("Plot_FFmsyTRUE_Rec.png", plot = Plot_FFmsyTrue_Rec,
#        scale = 1, width = 14, height = 7, units = "in")




## FM ####
Rec_FM <- readRDS("files/Rec_FM.rds")

Plot.Rec_FM <- ggplot(Rec_FM, aes(x = Method, y = Relative_error, fill = Method)) +
  geom_violin(trim = F) +
  labs(title = "Fishing mortality", y = "Relative Error") +
  stat_summary(fun.y = median, geom = "point", size = 2, color = "black") +
  geom_hline(yintercept = 0, color = 'darkgrey', size = 1) +
  facet_wrap(~ Scenario, ncol = 3, scale = "fixed") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 20, colour = "black")) +
  theme(axis.text.y = element_text(size = 20))
Plot.Rec_FM


Rec_FM.err <- readRDS("files/Rec_FM_err.rds")

tab <- data.frame(
  TB = round(Rec_FM.err$TB[1, c(1:2)],4),
  LBSPR = round(Rec_FM.err$LBSPR[1, c(1:2)],4),
  LIME = round(Rec_FM.err$LIME[1, c(1:2)],4),
  LBRA = round(Rec_FM.err$LBRA[1, c(1:2)],4),
  TB = round(Rec_FM.err$TB[2, c(1:2)],4),
  LBSPR = round(Rec_FM.err$LBSPR[2, c(1:2)],4),
  LIME = round(Rec_FM.err$LIME[2, c(1:2)],4),
  LBRA = round(Rec_FM.err$LBRA[2, c(1:2)],4),
  TB = round(Rec_FM.err$TB[3, c(1:2)],4),
  LBSPR = round(Rec_FM.err$LBSPR[3, c(1:2)],4),
  LIME = round(Rec_FM.err$LIME[3, c(1:2)],4),
  LBRA = round(Rec_FM.err$LBRA[3, c(1:2)],4)
)

colnames(tab) <- c(expression(TB[None]),
                   expression(LBSPR[None]),
                   expression(LIME[None]),
                   expression(LBRA[None]),
                   expression(TB[Error]),
                   expression(LBSPR[Error]),
                   expression(LIME[Error]),
                   expression(LBRA[Error]),
                   expression(TB[AR]),
                   expression(LBSPR[AR]),
                   expression(LIME[AR]),
                   expression(LBRA[AR]))
rownames(tab) <- c("MRE", "MARE")


tt <- ttheme_default(coRecead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(tab, theme=tt)


Plot.Rec_FM1 <- grid.arrange(Plot.Rec_FM, tbl,
                            nrow = 2,
                            as.table = TRUE,
                            heights = c(3,1))

# ggsave("Plot_FMTRUE_Rec.png", plot = Plot_FMTrue_Rec,
#        scale = 1, width = 14, height = 7, units = "in")






## Fmsy ####
Rec_Fmsy <- readRDS("files/Rec_Fmsy.rds")

Plot.Rec_Fmsy <- ggplot(Rec_Fmsy, aes(x = Method, y = Relative_error, fill = Method)) +
  geom_violin(trim = F) +
  labs(title = "Fmsy", y = "Relative Error") +
  stat_summary(fun.y = median, geom = "point", size = 2, color = "black") +
  geom_hline(yintercept = 0, color = 'darkgrey', size = 1) +
  facet_wrap(~ Scenario, ncol = 3, scale = "fixed") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 20, colour = "black")) +
  theme(axis.text.y = element_text(size = 20))
Plot.Rec_Fmsy


Rec_Fmsy.err <- readRDS("files/Rec_Fmsy_err.rds")

tab <- data.frame(
  TB = round(Rec_Fmsy.err$TB[1, c(1:2)],4),
  LBSPR = round(Rec_Fmsy.err$LBSPR[1, c(1:2)],4),
  LIME = round(Rec_Fmsy.err$LIME[1, c(1:2)],4),
  LBRA = round(Rec_Fmsy.err$LBRA[1, c(1:2)],4),
  TB = round(Rec_Fmsy.err$TB[2, c(1:2)],4),
  LBSPR = round(Rec_Fmsy.err$LBSPR[2, c(1:2)],4),
  LIME = round(Rec_Fmsy.err$LIME[2, c(1:2)],4),
  LBRA = round(Rec_Fmsy.err$LBRA[2, c(1:2)],4),
  TB = round(Rec_Fmsy.err$TB[3, c(1:2)],4),
  LBSPR = round(Rec_Fmsy.err$LBSPR[3, c(1:2)],4),
  LIME = round(Rec_Fmsy.err$LIME[3, c(1:2)],4),
  LBRA = round(Rec_Fmsy.err$LBRA[3, c(1:2)],4)
)

colnames(tab) <- c(expression(TB[None]),
                   expression(LBSPR[None]),
                   expression(LIME[None]),
                   expression(LBRA[None]),
                   expression(TB[Error]),
                   expression(LBSPR[Error]),
                   expression(LIME[Error]),
                   expression(LBRA[Error]),
                   expression(TB[AR]),
                   expression(LBSPR[AR]),
                   expression(LIME[AR]),
                   expression(LBRA[AR]))
rownames(tab) <- c("MRE", "MARE")


tt <- ttheme_default(coRecead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(tab, theme=tt)


Plot.Rec_Fmsy1 <- grid.arrange(Plot.Rec_Fmsy, tbl,
                              nrow = 2,
                              as.table = TRUE,
                              heights = c(3,1))

# ggsave("Plot_FmsyTRUE_Rec.png", plot = Plot_FmsyTrue_Rec,
#        scale = 1, width = 14, height = 7, units = "in")





## SPRmsy ####
Rec_SPRmsy <- readRDS("files/Rec_SPRmsy.rds")

Plot.Rec_SPRmsy <- ggplot(Rec_SPRmsy, aes(x = Method, y = Relative_error, fill = Method)) +
  geom_violin(trim = F) +
  labs(title = "SPRmsy", y = "Relative Error") +
  stat_summary(fun.y = median, geom = "point", size = 2, color = "black") +
  geom_hline(yintercept = 0, color = 'darkgrey', size = 1) +
  facet_wrap(~ Scenario, ncol = 3, scale = "fixed") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 20, colour = "black")) +
  theme(axis.text.y = element_text(size = 20))
Plot.Rec_SPRmsy


Rec_SPRmsy.err <- readRDS("files/Rec_SPRmsy_err.rds")

tab <- data.frame(
  TB = round(Rec_SPRmsy.err$TB[1, c(1:2)],4),
  LBSPR = round(Rec_SPRmsy.err$LBSPR[1, c(1:2)],4),
  LIME = round(Rec_SPRmsy.err$LIME[1, c(1:2)],4),
  LBRA = round(Rec_SPRmsy.err$LBRA[1, c(1:2)],4),
  TB = round(Rec_SPRmsy.err$TB[2, c(1:2)],4),
  LBSPR = round(Rec_SPRmsy.err$LBSPR[2, c(1:2)],4),
  LIME = round(Rec_SPRmsy.err$LIME[2, c(1:2)],4),
  LBRA = round(Rec_SPRmsy.err$LBRA[2, c(1:2)],4),
  TB = round(Rec_SPRmsy.err$TB[3, c(1:2)],4),
  LBSPR = round(Rec_SPRmsy.err$LBSPR[3, c(1:2)],4),
  LIME = round(Rec_SPRmsy.err$LIME[3, c(1:2)],4),
  LBRA = round(Rec_SPRmsy.err$LBRA[3, c(1:2)],4)
)

colnames(tab) <- c(expression(TB[None]),
                   expression(LBSPR[None]),
                   expression(LIME[None]),
                   expression(LBRA[None]),
                   expression(TB[Error]),
                   expression(LBSPR[Error]),
                   expression(LIME[Error]),
                   expression(LBRA[Error]),
                   expression(TB[AR]),
                   expression(LBSPR[AR]),
                   expression(LIME[AR]),
                   expression(LBRA[AR]))
rownames(tab) <- c("MRE", "MARE")


tt <- ttheme_default(coRecead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(tab, theme=tt)


Plot.Rec_SPRmsy1 <- grid.arrange(Plot.Rec_SPRmsy, tbl,
                               nrow = 2,
                               as.table = TRUE,
                               heights = c(3,1))

# ggsave("Plot_SPRmsyTRUE_Rec.png", plot = Plot_SPRmsyTrue_Rec,
#        scale = 1, width = 14, height = 7, units = "in")




