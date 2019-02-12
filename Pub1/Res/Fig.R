

# Packages and Objects to save --------------------------------------------


if(!require(ggplot2)) install.packages("ggplot"); require(ggplot2)
if(!require(gridExtra)) install.packages("gridExtra"); require(gridExtra)



# LH comparison -----------------------------------------------------------

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







# Rec comparison ----------------------------------------------------------


