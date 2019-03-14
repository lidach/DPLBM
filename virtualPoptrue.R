
#' @title Virtual fish population true values
#'
#' @param K.mu mean K (growth parameter from von Bertalanffy growth function)
#' @param K.cv coefficient of variation on K
#' @param Linf.mu mean Linf (infinite length parameter from von Bertalanffy growth function)
#' @param Linf.cv coefficient of variation on Linf
#' @param t0 theoretical age at length 0
#' @param ts summer point (range 0 to 1) (parameter from seasonally oscillating von Bertalanffy growth function)
#' @param C strength of seasonal oscillation (range 0 to 1) (parameter from seasonally oscillating von Bertalanffy growth function)
#' @param LWa length-weight relationship constant 'a' (W = a*L^b). Model assumed length in cm and weight in kg.
#' @param LWb length-weight relationship constant 'b' (W = a*L^b). Model assumed length in cm and weight in kg.
#' @param Lmat.f length at maturity for females (where 50\% of individuals are mature)
#' @param wmat.f width between 25\% and 75\% quantiles for Lmat for females
#' @param Lmat.m length at maturity for males(where 50\% of individuals are mature)
#' @param wmat.m width between 25\% and 75\% quantiles for Lmat for males
#' @param rmaxBH parameter for Beverton-Holt stock recruitment relationship (see \code{\link[fishdynr]{srrBH}})
#' @param betaBH parameter for Beverton-Holt stock recruitment relationship (see \code{\link[fishdynr]{srrBH}})
#' @param srr.cv coefficient of variation stock recruitment relationship
#' @param repro_wt weight of reproduction (vector of monthly reproduction weight)
#' @param M natural mortality
#' @param Etf  effort (E = F / q); single numeric, numeric vector for effort per year, or matrix for different fleets (columns) and different years (rows)
#' @param qtf catchability (default 0.005); single numeric, numeric vector for effort per year, or matrix for different fleets (columns)  and different years (rows)
#' @param harvest_rate Fishing mortality (i.e. 'F' = C/B); if NaN Etf and qtf are used to estimate the harvest_rate
#' @param gear_types Character(s) defining the gear of the fishing fleet(s) (so far: either "trawl" or "gillnet")
#' @param L50 minimum length of capture (in cm). Where selectivity equals 0.5. Assumes logistic ogive typical of trawl net selectivity.
#' @param wqs width of selectivity ogive (in cm)
#' @param sel_list list with selectivities parameters for gillnet selectivity
#' @param bin.size resulting bin size for length frequencies (in cm)
#' @param timemin time at start of simulation (in years). Typically set to zero.
#' @param timemax time at end of simulation (in years).
#' @param timemin.date date corresponding to timemin (of "Date" class)
#' @param tincr time increment for simulation (default = 1/12; i.e. 1 month)
#' @param N0 starting number of individuals
#' @param fished_t times when stock is fished; when NA no exploitation simulated
#' @param lfqFrac fraction of fished stock that are sampled for length frequency data (default = 0.1).
#' @param spmYears number of years which are used for fitting the surplus production model (Default 10)
#' @param binSizeVPA bin size for the application of VPA (to derive reference points in ypr)
#' @param progressBar Logical. Should progress bar be shown in console (Default=TRUE)
#' @param plot Logical. Should the standard plots be printed (Default=TRUE)
#' @param seed integer; indicating the seed for set.seed()
#' @param SSBset set SSB level
#'
#' @description See \code{\link[fishdynr]{dt_growth_soVB}} for information on growth function.
#' The model creates variation in growth based on a mean phi prime value for the population,
#' which describes relationship between individual Linf and K values. See Vakily (1992)
#' for more details.
#'
#' @details The model takes around 5 to 10 years to reach equilibrium, i.e. no biomass changes independent from fishing activity, the actual time is dependent on N0, K.mu, Lmat, repro_wt  and rmax.BH. For the estimation of carrying capacity the first 10 years of the simulation are disregarded and only subsequent years where no fishing took place are used to estimate the annual mean carrying capacity (K). If fishing is simulated for all years or fishing activities start before ten years after simulation start no carrying capacity is estimated.
#'
#' @return a list containing ypr, pbr, and spr
#'
#' @references
#' Vakily, J.M., 1992. Determination and comparison of bivalve growth,
#' with emphasis on Thailand and other tropical areas. WorldFish.
#'
#' Munro, J.L., Pauly, D., 1983. A simple method for comparing the growth
#' of fishes and invertebrates. Fishbyte 1, 5-6.
#'
#' Pauly, D., Munro, J., 1984. Once more on the comparison of growth
#' in fish and invertebrates. Fishbyte (Philippines).
#'
#' @importFrom graphics hist
#' @importFrom stats rlnorm runif weighted.mean
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats qnorm rnorm
#' @importFrom TropFishR VBGF
#'
#' @export
#'
virtualPoptrue <- function(tincr = 1/12,
                            K.mu = 0.5, K.cv = 0.1,
                            Linf.mu = 80, Linf.cv = 0.1,
                            t0 = -0.03,
                            ts = 0, C = 0.85,
                            LWa = 0.01, LWb = 3,   ## for cm in g? then final B results also in g
                            Lmat.f = 0.5*Linf.mu, wmat.f = Lmat.f*0.2,
                            Lmat.m = 0.45*Linf.mu, wmat.m = Lmat.m*0.15,
                            rmaxBH = 1e6,
                            betaBH = 1, srr.cv = 0.1,
                            repro_wt = c(0,0,0,1,0,0,0,0,0,0,0,0),
                            M = 0.7,
                            Etf = 500,
                            qtf = 0.001,
                            harvest_rate = NaN,
                            gear_types = "trawl",   # alternative: "gillnet", or "dome"
                            L50 = 0.25*Linf.mu,
                            wqs = L50*0.2,
                            sel_list = list(mesh_size=100, mesh_size1=60,select_dist="lognormal",select_p1=3, select_p2=0.5),  # parameters adapted from the tilapia data set (increased spread)
                            bin.size = 1,
                            timemin = 0, timemax = 30, timemin.date = as.Date("1980-01-01"),
                            N0 = 0, ##1000,
                            fished_t = seq(0,100,tincr),
                            lfqFrac = 1,
                            numSamp = NA,
                            spmYears = 10,
                            addSurvey = FALSE,
                            survey_t = seq(17,25,1/6),
                            survey_L50 = 0.1 * Linf.mu,
                            survey_wqs = survey_L50 * 0.2,
                            numSampSurvey = 500,
                            binSizeVPA = bin.size,
                            progressBar = TRUE,
                            plot = TRUE,
                            seed = NULL,
                            ## new pars
                            SSBset = 1e10){

  ## Fishing mortality - effort - catchability
  ## if E == single value, assuming one fleet and same effort for all fished years
  if(length(as.numeric(Etf))==1){
    Emat <- as.matrix(rep(Etf,length(fished_t)))
  }
  ## if E == matrix, rows = years and columns = fleets


  if(class(Etf) == "matrix"){
    Emat <- Etf
  }else if(length(Etf)>1){
    Emat <- as.matrix(Etf)
  }

  ## adapt q to Emat
  if(length(qtf)==1){
    qmat <- matrix(qtf, ncol=dim(Emat)[2], nrow=dim(Emat)[1])
  }
  if(class(qtf) == "matrix"){
    if(dim(qtf)[1] != dim(Emat)[1]){
      qmat <- matrix(rep(qft[1,], dim(Emat)[1]),ncol=dim(qft)[2],byrow=TRUE)
    }else{
      qmat <- qtf
    }
  }else if(length(qtf)>1){
    qmat <- as.matrix(qtf)
  }

  ## If no harvest_rate provided assuming that effort * catchability = fishing mortality
  if(!is.na(harvest_rate[1]) & !is.nan(harvest_rate[1])){
    if(length(as.numeric(harvest_rate))==1){
      harvest_rate <- rep(harvest_rate, length(fished_t))
    }else{
      harvest_rate <- matrix(rep(harvest_rate, each = length(fished_t)),
                             ncol = length(harvest_rate), nrow = length(fished_t))
    }
  }else{
    harvest_rate <- Emat * qmat
  }

  selfunc <- function(Lt, fleetNo){
    if(is.na(fleetNo)){
      gear_typesX <- gear_types
      L50X <- L50
      wqsX <- wqs
      sel_listX <- sel_list
    }else{
      gear_typesX <- gear_types[fleetNo]
    }
    switch(gear_typesX,
           trawl ={
             if(!is.na(fleetNo)){
               L50X <- L50[fleetNo]
               wqsX <- wqs[fleetNo]
             }
             pSel <- logisticSelect(Lt=Lt, L50=L50X, wqs=wqsX)
           },
           gillnet={
             if(!is.na(fleetNo)){
               sel_listX <- sel_list[[fleetNo]]
             }
             pSel <- do.call(fishdynr::gillnet, c(list(Lt=Lt),sel_listX))
           },
           dome={
             if(!is.na(fleetNo)){
               sel_listX <- sel_list[[fleetNo]]
             }
             pSel <- do.call(fishdynr::domeSelect, c(list(Lt=Lt),sel_listX))
           },
           stop(paste("\n",gear_typesX,"not recognized, possible options are: \n","trawl \n","gillnet \n")))
    return(pSel)
  }


  ## ## if multiple fleets target the same stock, the harvest rate of each fleet is scaled according to the combined harvest rate - this only works if all fleets would have the same gear!
  ## if(class(harvest_rate) == "matrix"){
  ##     multimat <- harvest_rate / rowSums(harvest_rate)
  ##     harvest_rate <- rowSums(harvest_rate * multimat)
  ##    }
  ## times
  timeseq = seq(from=timemin, to=timemax, by=tincr)
  if(!zapsmall(1/tincr) == length(repro_wt)) stop("length of repro_wt must equal the number of tincr in one year")
  repro_wt <- repro_wt/sum(repro_wt)
  repro_t <- rep(repro_wt, length=length(timeseq))
  # repro_t <- seq(timemin+repro_toy, timemax+repro_toy, by=1)

  # make empty lfq object
  lfq <- vector(mode="list", length(timeseq))
  names(lfq) <- timeseq

  indsSamp <- vector(mode="list", length(timeseq))
  names(indsSamp) <- timeseq

  ## make empty lfq object
  lfqSurvey <- vector(mode="list", length(timeseq))
  names(lfqSurvey) <- timeseq

  indsSampSurvey <- vector(mode="list", length(timeseq))
  names(indsSampSurvey) <- timeseq


  # Estimate tmaxrecr
  tmaxrecr <- (which.max(repro_wt)-1)*tincr

  # mean phiprime
  phiprime.mu = log10(K.mu) + 2*log10(Linf.mu)



  # required functions ------------------------------------------------------
  date2yeardec <- function(date){as.POSIXlt(date)$year+1900 + (as.POSIXlt(date)$yday)/365}
  yeardec2date <- function(yeardec){as.Date(strptime(paste(yeardec%/%1, ceiling(yeardec%%1*365+1), sep="-"), format="%Y-%j"))}

  make.inds <- function(
    id=NaN, A = 0, L = 0, W=NaN, sex = NaN, mat=0,
    K = K.mu, Winf=NaN, Linf=NaN, phiprime=NaN,
    F=NaN, Z=NaN, Fd=0, alive=1, FSurvey = 0
  ){
    inds <- data.frame(
      id = id,
      A = A,
      L = L,
      W = W,
      sex = sex,
      Lmat=NaN,
      mat = mat,
      K = K,
      Linf = Linf,
      Winf = Winf,
      phiprime = phiprime,
      F = F,
      Z = Z,
      Fd = Fd,
      alive = alive,
      FSurvey = FSurvey
    )
    lastID <<- max(inds$id)
    return(inds)
  }

  express.inds <- function(inds){
    inds$Linf <- exp(rnorm(nrow(inds), log(Linf.mu), Linf.cv)) ## old:  Linf.mu * rlnorm(nrow(inds), 0, Linf.cv)
    inds$Winf <- LWa*inds$Linf^LWb
    # inds$K <- 10^(phiprime.mu - 2*log10(inds$Linf)) * rlnorm(nrow(inds), 0, K.cv)
    inds$K <- exp(rnorm(nrow(inds), log(K.mu), K.cv)) ## old: K.mu * rlnorm(nrow(inds), 0, K.cv)
    inds$W <- LWa*inds$L^LWb
    inds$phiprime <- log10(inds$K) + 2*log10(inds$Linf)

    inds$sex <- rbinom(nrow(inds), size = 1, prob = .5)  # 0 = F, 1 = M
    inds$Lmat[which(inds$sex == 0)] <- rnorm(length(which(inds$sex == 0)), mean=Lmat.f, sd=wmat.f/diff(qnorm(c(0.25, 0.75))))
    inds$Lmat[which(inds$sex == 1)] <- rnorm(length(which(inds$sex == 1)), mean=Lmat.m, sd=wmat.m/diff(qnorm(c(0.25, 0.75))))
    inds$L <- dt_growth_soVB(Linf = inds$Linf, K = inds$K, ts = ts, C = C, L1 = 0, t1 = t0, t2 = 0)
    return(inds)
  }

  grow.inds <- function(inds){
    # grow
    L2 <- dt_growth_soVB(Linf = inds$Linf, K = inds$K, ts = ts, C = C, L1 = inds$L, t1 = tj-tincr, t2 = tj)
    # update length and weight
    inds$L <- L2
    inds$W <- LWa*inds$L^LWb
    # age inds
    inds$A <- inds$A + tincr
    return(inds)
  }

  mature.inds <- function(inds){
    # p <- pmat_w(inds$L, Lmat, wmat) # probability of being mature at length L
    # p1t <- 1-((1-p)^tincr)
    # inds$mat <- ifelse(runif(nrow(inds)) < p1t | inds$mat == 1, 1, 0)
    inds$mat <- ifelse((inds$L > inds$Lmat | inds$mat == 1), 1, 0)
    return(inds)
  }

  reproduce.inds <- function(inds,  save = FALSE){
    ## reproduction can only occur of population contains >1 mature individual
    if(repro > 0){
      ## calc. SSB
      SSB <- SSBset ##sum(inds$W*inds$mat, na.rm = TRUE)
      n.recruits <- ceiling(srrBH(rmaxBH, betaBH, SSB) * repro)
      ## add noise to recruitment process
      n.recruits <- exp(rnorm(1, log(n.recruits), srr.cv)) ## old: n.recruits * rlnorm(1, 0, sdlog = srr.cv)
      ## save SSB + n.recruits for stock recruitment plot
      if(save) stockRec <<- rbind(stockRec, data.frame(SSB = SSB, recruits = n.recruits, time = tj))
      ## make recruits
      offspring <- make.inds(
        id = seq(lastID+1, length.out=n.recruits)
      )
      ## express genes in recruits
      offspring <- express.inds(offspring)
      ##combine all individuals
      inds <- rbind(inds, offspring)
    }
    return(inds)
  }

  death.inds <- function(inds, f0 = FALSE){
    ## multiple fleets
    if(class(harvest_rate)=="matrix"){
      if(dim(harvest_rate)[2]>=2){
        pSel <- matrix(NaN, ncol=dim(harvest_rate)[2],nrow=dim(inds)[1])
        for(seli in 1:(dim(harvest_rate)[2])){
          pSel[,seli] <- selfunc(Lt = inds$L, fleetNo = seli)
        }
        ## effective fishing mortality (in relation to selectivity) - per fleet with mutliple fleets
        Feff <- pSel * Fmax
        ## single fishing mortality value (per year) scaled according to F of each fleet
        ## this calculation only works if there is fishery (Fmax in denominator not allowed to be 0, otherwise F = NaN and then Z = NaN), thus:
        if(all(Fmax == 0)){
          inds$F <- 0
        }else{
          inds$F <- as.numeric(rowSums(Feff * Fmax) / sum(Fmax))
        }
      }else{
        ## single fleet
        pSel <- selfunc(Lt = inds$L, fleetNo = NA)
        inds$F <- as.numeric(pSel * Fmax)
      }
    }else{
      ## single fleet
      pSel <- selfunc(Lt = inds$L, fleetNo = NA)
      inds$F <- as.numeric(pSel * Fmax)
    }

    ## vulnerability to survey
    pSelSurvey <- logisticSelect(Lt=inds$L, L50=survey_L50, wqs=survey_wqs)
    inds$FSurvey <- as.numeric(pSelSurvey * 1)

    inds$Z <- M + inds$F
    if(f0) inds$Z <- M
    pDeath <- 1 - exp(-inds$Z*tincr)
    dead <- which(runif(nrow(inds)) < pDeath)
    # determine if natural or fished
    if(length(dead) > 0){
      inds$alive[dead] <- 0
      tmp <- cbind(inds$F[dead], inds$Z[dead])
      # Fd=1 for fished individuals; Fd=0, for those that died naturally
      Fd <- apply(tmp, 1, FUN=function(x){sample(c(0,1), size=1, prob=c(M/x[2], x[1]/x[2]) )})
      inds$Fd[dead] <- Fd
      rm(tmp)
    }
    return(inds)
  }

  remove.inds <- function(inds){
    dead <- which(inds$alive == 0)
    if(length(dead)>0) {inds <- inds[-dead,]}
    return(inds)
  }

  record.inds <- function(inds, ids=1:10, rec=NULL){
    if(is.null(rec)) {
      rec <- vector(mode="list", length(ids))
      names(rec) <- ids
      inds <- inds
    } else {
      ids <- as.numeric(names(rec))
    }
    if(length(rec) > 0) {
      inds.rows.rec <- which(!is.na(match(inds$id, ids)))
      if(length(inds.rows.rec) > 0){
        for(ii in inds.rows.rec){
          match.id <- match(inds$id[ii], ids)
          if(is.null(rec[[match.id]])) {
            rec[[match.id]] <- inds[ii,]
          } else {
            rec[[match.id]] <- rbind(rec[[match.id]], inds[ii,])
          }
        }
      }
    }
    return(rec)
  }

  spmOpt <- function(x, B){
    K = x[1]
    r = x[2]
    n = x[3]
    Bthat <- rep(NA, length(B))
    Bthat[1] <- B[1]
    for(i in 2:length(Bthat)){
      Bthat[i] <- Bthat[i-1] + ((r / (n - 1)) * Bthat[i-1] * (1 - (Bthat[i-1] / K)^(n-1)))*tincr
    }
    sum((B - Bthat)^2)
  }


  spmPlot <- function(pars){
    B0 <- pars[1]
    K <- pars[2]
    r <- pars[3]
    n <- pars[4]
    Bthat <- rep(NA, length(resf0$pop$B))
    Bthat[1] <- B0
    for(i in 2:length(Bthat)){
      Bthat[i] <- Bthat[i-1] + ((r / (n - 1)) * Bthat[i-1] * (1 - (Bthat[i-1] / K)^(n-1)))*tincr
    }
    Bthat
  }


  ## run model ---------------------------------------------------------------

  ## seed values
  if(is.null(seed) || is.null(seed)){
    seed <- floor(runif(1,1,1000))
  }
  set.seed(seed)

  ## Initial population
  lastID <- 0
  inds <- make.inds(
    id=0 ##seq(N0)
  )
  inds <- inds[-1,]
  inds <- express.inds(inds = inds)

  ## results object
  res <- list()
  res$pop <- list(
    dates = yeardec2date( date2yeardec(timemin.date) + (timeseq - timemin) ),
    N = NaN*timeseq,
    B = NaN*timeseq,
    SSB = NaN*timeseq
  )

  stockRec <- as.data.frame(matrix(ncol=3,nrow=0))
  colnames(stockRec) <- c("SSB", "recruits","time")

  ## for total catches per time step (e.g. for correction factor in VPA)
  catches <- vector("numeric",length(timeseq))
  ypr <- vector("numeric",length(timeseq))
  bpr <- vector("numeric",length(timeseq))
  spr <- vector("numeric",length(timeseq))

  j <- 1

  while(j < max(seq(timeseq)) && (j < 1/tincr | nrow(inds) > 0)){

    tj <- timeseq[j]


    ## harvest rate applied? lfq sampled?
    if(is.na(fished_t[1]) | is.nan(fished_t[1])){ ## before: length(fished_t) == 0 : as I see it fished_t never has length 0, even if set ot NA or NaN, it woudl have length 1
      Fmax <- 0
      lfqSamp <- 0
    } else if(min(sqrt((tj-fished_t)^2)) < 1e-8){
      ## time index for fished_t
      tfish <- which.min(abs(fished_t - tj))
      ## provide yearly Fmax value (per fleet if multiple fleets simulated)
      if(class(harvest_rate) == "matrix"){
        Fmax <- harvest_rate[tfish,]
      }else if(length(harvest_rate)>1){
        Fmax <- harvest_rate[tfish]
      }else{
        Fmax <- harvest_rate
      }
      lfqSamp <- 1
    } else {
      Fmax <- 0
      lfqSamp <- 0
    }
    if(addSurvey & min(sqrt((tj-survey_t)^2)) < 1e-8){
      lfqSampSurvey <- 1
    }else lfqSampSurvey <- 0

    repro <- repro_t[j]

    ## population processes
    inds <- grow.inds(inds)
    inds <- mature.inds(inds)
    if(tj < 1) inds <- reproduce.inds(inds = inds, save = TRUE)
    inds <- death.inds(inds)

    ypr[j] <- sum(inds$W[inds$Fd == 1])

    inds <- remove.inds(inds)

    bpr[j] <- sum(inds$W)
    spr[j] <- sum(inds$W[inds$mat == 1])

    j <- j + 1

  }

  return(list(ypr = sum(ypr)/rmaxBH, bpr = sum(bpr)/rmaxBH, spr = sum(spr)))
}
