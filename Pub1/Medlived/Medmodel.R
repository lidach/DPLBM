## Packages ############################
# devtools::install_github("tokami/TropFishR")
# devtools::install_github("tokami/fishdynr")

require(TropFishR)
require(fishdynr)
require(parallel)

#' run virtualPop2 externally





#'---
## Med Model Simulation ############################
source("~/DPLBM/virtualPop_par.R")

#' OM inputs
Linf.mu = 64.6
Linf.cv = 0.1
K.mu = 0.21
K.cv = 0.1
ts = 0
C = 0
t0 = -0.01
M = 0.32
burnin = 10
tincr = 1/12
fished_t = seq(0,35,1/12)
timemax = 35
Fishscen = "Constant"
F_scen = M/2.5
harvest_rate = NaN
gear_types = "trawl"
L50 = 0.25*Linf.mu
wqs = L50*0.2 # definition
Lmat.f = 34
Lmat.m = 34
wmat.f = Lmat.f*0.2 # definition
wmat.m = Lmat.m*0.2
LWa = 0.018
LWb = 2.895
rec_dyn = "BH"
rmaxBH = 1e4
betaBH = 1
SigmaR = 0.01
rho = 0.43
numSamp = 200
N0 = 1000
repro_wt = c(0,0,1,1,1,0,0,0,0,0,0,0)
# modpath = "D:/DPLBM/Pub1/Medlived"
modpath = getwd()
iteration = 1:300

#' Run Model
set.seed(1)
lutGu <- virtualPop_par(Linf.mu = Linf.mu,
                    Linf.cv = Linf.cv,
                    K.mu = K.mu,
                    K.cv = K.cv,
                    ts = ts,
                    C = C,
                    t0 = t0,
                    M = M,
                    tincr = tincr,
                    Fishscen = "Constant",
                    F_scen = F_scen,
                    harvest_rate = harvest_rate,
                    gear_types = "trawl",
                    L50 = L50,
                    wqs = wqs,
                    fished_t = fished_t,
                    timemax = timemax,
                    Lmat.f = Lmat.f,
                    wmat.f = wmat.f,
                    Lmat.m = Lmat.m,
                    wmat.m = wmat.m,
                    LWa= LWa,
                    LWb= LWb,
                    rec_dyn = "BH",
                    rmaxBH = rmaxBH,
                    betaBH = betaBH,
                    SigmaR = SigmaR,
                    N0= N0,
                    repro_wt= repro_wt,
                    addSurvey = FALSE,
                    numSamp = 200,
                    progressBar = FALSE,
                    plot = FALSE,
                    modpath = modpath,
                    iteration = iteration,
                    seed = NULL)


setwd(modpath)
iters <- length(iteration)

Medmodel_ls <- list.files(path = modpath, recursive = T, pattern = "True.rds")
Medmodel <- lapply(Medmodel_ls, function(x) {
  dir = readRDS(file = x)
  Medmodel = dir$lfqbin
})


diri <- list()
for (i in 1:iters){
  diri[[i]] <- readRDS(file = Medmodel_ls[i])
  Medmodel[[i]]$catch <- Medmodel[[i]]$catch[,which(format(Medmodel[[i]]$dates, "%Y") %in% "2014")]
  Medmodel[[i]]$dates <- Medmodel[[i]]$dates[format(Medmodel[[i]]$dates, "%Y") %in% "2014"]
  class(Medmodel[[i]]) <- "lfq"
  Medmodel[[i]]$par <- list(Linf = Linf.mu,
                            K = K.mu,
                            t_anchor = 0.2840183)
}



mean(unlist(lapply(diri,
                   function(x)
                     x$refLev$states$SPR[which(rownames(diri[[i]]$refLev$states) == "2014")])))


saveRDS(Medmodel, file = "LFQmedmodel.rds")



#'---
#' Extracting parameters
ref <- list.files(path = modpath, recursive = T, pattern = "True.rds")
ref <- lapply(ref, function(x){
  dir = readRDS(file = x)
  ref = dir$refLev
})

for(i in 1:length(ref)){
  ref[[i]][["years"]] <- ref[[i]][["years"]][25]
  ref[[i]][["SPR"]] <- ref[[i]][["SPR"]][25]
  ref[[i]][["Pmega"]] <- ref[[i]][["Pmega"]][25]
  ref[[i]][["Popt"]] <- ref[[i]][["Popt"]][25]
  ref[[i]][["Pmat"]] <- ref[[i]][["Pmat"]][25]
  ref[[i]][["Pobj"]] <- ref[[i]][["Pobj"]][25]
  ref[[i]][["F01"]] <- ref[[i]][["F01"]][25]
  ref[[i]][["Fmax"]] <- ref[[i]][["Fmax"]][25]
  ref[[i]][["F05"]] <- ref[[i]][["F05"]][25]
}


#' for SB
inds <- list.files(path = modpath, recursive = T, pattern = "True.rds")
inds <- lapply(inds, function(x){
  dir = readRDS(file = x)
  inds = dir$inds
})

listL <- list()
listA <- list()
listW <- list()
listLmat <- list()
listLinf <- list()
listK <- list()
for(i in 1:length(inds)){
  inds[[i]] <-  c(inds[[i]]["34"], inds[[i]]["34.0833333333333"], inds[[i]]["34.1666666666667"], inds[[i]]["34.25"], inds[[i]]["34.3333333333333"], inds[[i]]["34.4166666666667"], inds[[i]]["34.5"], inds[[i]]["34.5833333333333"], inds[[i]]["34.6666666666667"], inds[[i]]["34.75"], inds[[i]]["34.8333333333333"], inds[[i]]["34.9166666666667"])
  listL[[i]] <- sapply(inds[[i]], '[', 'L')
  listL[[i]] <- unname(unlist(listL[[i]], recursive = F))
  listA[[i]] <- sapply(inds[[i]], '[', 'A')
  listA[[i]] <- unname(unlist(listA[[i]], recursive = F))
  listW[[i]] <- sapply(inds[[i]], '[', 'W')
  listW[[i]] <- unname(unlist(listW[[i]], recursive = F))
  listLmat[[i]] <- sapply(inds[[i]], '[', 'Lmat')
  listLmat[[i]] <- unname(unlist(listLmat[[i]], recursive = F))
  listLinf[[i]] <- sapply(inds[[i]], '[', 'Linf')
  listLinf[[i]] <- unname(unlist(listLinf[[i]], recursive = F))
  listK[[i]] <- sapply(inds[[i]], '[', 'K')
  listK[[i]] <- unname(unlist(listK[[i]], recursive = F))
  inds[[i]] <- list(cbind(A = listA[[i]], L = listL[[i]], W = listW[[i]], Lmat = listLmat[[i]], Linf = listLinf[[i]], K = listK[[i]]))
}


Medmodel_list <- list(ref = ref, inds = inds)
saveRDS(Medmodel_list, file = "Medlist.rds")
saveRDS(inds, file = "inds.rds")


unlink(Medmodel_ls)
rm(list = ls())





