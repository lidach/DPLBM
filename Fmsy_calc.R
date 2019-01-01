## calculate expected abundance given constant fishing mortality
calc_abun <- function(ages, S_a, M, F, R0){

  N_a <- R0

  Fmat <- NA
  for(i in 1:length(S_a)){
    Fmat[i] <- S_a[i] * F
  }
  Ftotal <- Fmat

  dt <- diff(ages) # ******
  for(i in 2:length(ages)){
    if(i < length(ages)) N_a[i] <- N_a[i-1] * exp(-(M + Ftotal[i-1])*dt[i-1]) # ******
    if(i == length(ages)) N_a[i] <- (N_a[i-1] * exp(-(M + Ftotal[i-1])*dt[i-1]))/
        (1 - exp(-(M + Ftotal[i-1])*dt[i-1])) # ******
  }
  return(N_a)
}


calc_msy <- function(ages, S_a, M, F, R0, W_a){

  dt <- diff(ages)[1] # ******
  Nages <- calc_abun(ages = ages, S_a = S_a, M = M, F = F, R0 = R0)
  YPR <- sum(Nages * W_a * (1-exp(-M - F)) * (F) / (M + F)) * dt # ******

  return(YPR)
}

calc_spr <- function(ages, S_a, M, F, R0, SW_a){ #**************

  dt <- diff(ages)[1] # ******
  Nages <- calc_abun(ages = ages, S_a = S_a, M = M, F = F, R0 = R0)
  SPR <- sum(Nages * SW_a * (1-exp(-M - F)))* dt # ******

  return(SPR)
}


linf <- 64.6
vbk <- 0.21
tincr <- 1/12 # ******
ages <- seq(0,18,tincr) # ******
lwa = 0.018
lwb = 2.895
t0 <- -0.01
L_a <- linf*(1-exp(-vbk*(ages - t0)))
W_a <- lwa*L_a^lwb
M <- 0.32
SL50 <- 11
ML50 <- 34
R0 <- 1
CV <- 0.1
wqs <- SL50 * 0.2
wmat <- ML50 * 0.2
binwidth <- 1


mids <- seq((binwidth/2), linf*1.5, binwidth)
highs <- mids + (binwidth/2)
lows <- mids - (binwidth/2)

L_a <- linf*(1-exp(-vbk*(ages - t0)))
W_a <- lwa*L_a^lwb
W_l <- lwa*mids^lwb


M <- 0.32
R0 <- 1


#' probability
lbprobs <- function(mnl, sdl) return(pnorm(highs, mnl, sdl) - pnorm(lows, mnl, sdl))
vlprobs <- Vectorize(lbprobs, vectorize.args = c("mnl", "sdl"))
plba_a <- t(vlprobs(L_a, L_a * CV))
plba_a <- plba_a / rowSums(plba_a)


#' selectivity and maturity at age
S_l <- 1 / (1 + exp(-(mids - SL50) /
                      (wqs / ( log(0.75/(1-0.75)) - log(0.25/(1-0.25)) ))) )
S_a <- apply(t(plba_a)*S_l, 2, sum)
Mat_l <- 1 / (1 + exp(-(mids - ML50) /
                        (wmat / ( log(0.75/(1-0.75)) - log(0.25/(1-0.25)) ))) )
Mat_a <- apply(t(plba_a)*Mat_l, 2, sum)

SW_a <- W_a*Mat_a #******************



#
n <- 300
MSY_plot2 <- data.frame(val = NA, F =  seq(0,3,length.out = n))
SPR_plot2 <- MSY_plot2
for(i in 1:300){
  MSY_plot2$val[i] <- calc_msy(ages = ages, S_a = S_a, M = M , F = MSY_plot2$F[i], R0 = R0, W_a = W_a)
  SPR_plot2$val[i] <- calc_spr(ages = ages, S_a = S_a, M = M , F = SPR_plot2$F[i], R0 = R0, SW_a = SW_a)
}

plot(val ~ F, MSY_plot2, type = "l", lwd = 2, lty = 1, ylab = "YPR")
Fmsy <- MSY_plot2$F[which.max(MSY_plot2$val)]
Fmsy

plot(val/max(val) ~ F, SPR_plot2, type = "l", lwd = 2, lty = 1, ylab = "SPR")
(SPR_plot2$val/max(SPR_plot2$val))[which.max(MSY_plot2$val)]
