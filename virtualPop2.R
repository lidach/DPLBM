#' @title Virtual fish population
#'
#' @param burnin operating model burn-in period (default- 10 years)
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
#' @param rec_dyn patterns for recruitment dynamics. Options include "BH" (Beverton-Holt), "Constant", and "AR" (autocorrelated recruitment)
#' @param rmaxBH parameter for Beverton-Holt stock recruitment relationship (see \code{\link[fishdynr]{srrBH}})
#' @param betaBH parameter for Beverton-Holt stock recruitment relationship (see \code{\link[fishdynr]{srrBH}})
#' @param srr.cv coefficient of variation stock recruitment relationship
#' @param SigmaR recruitment standard deviation- default: 0.737, the median across all fish species (Thorson et al. 2014). Used for generating recruitment deviates in the simulation
#' @param rho changing the value of rho will add autocorrelation to the generated recruitment time series for a simulation study
#' @param repro_wt weight of reproduction (vector of monthly reproduction weight)
#' @param M natural mortality
#' @param Etf  effort (E = F / q); single numeric, numeric vector for effort per year, or matrix for different fleets (columns) and different years (rows)
#' @param qtf catchability (default 0.005); single numeric, numeric vector for effort per year, or matrix for different fleets (columns)  and different years (rows)
#' @param Fishscen patterns for fishing dynamics. Options include "None" (then uses harvest_rate, Etf, and qtf), "Constant", "Fish_down" (fishing down), "Two-way", and "Multi_peak". Depends on variables F_scen, F_high, F_low, and SigmaF 
#' @param SigmaF fishing mortality standard deviation. Used for generating fishing mortality deviates in the simulation
#' @param F_scen Some reference point (Fmsy, F40, etc) used to create different fishing scenarios
#' @param F_high Upper bound of fishing mortality with F_scen used to create different fishing scenarios
#' @param F_low Lower bound of fishing mortality with F_scen used to create different fishing scenarios
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
#' @param modpath model path for saving simulated populations; can be set as NULL to run within R environment only without saving locally
#' @param iteration vector of iterations of simulated data (default = 1, i.e. run only 1 iteration)
#'
#' @description See \code{\link[fishdynr]{dt_growth_soVB}} for information on growth function.
#' The model creates variation in growth based on a mean phi prime value for the population,
#' which describes relationship between individual Linf and K values. See Vakily (1992)
#' for more details.
#'
#' @details The model takes around 5 to 10 years to reach equilibrium, i.e. no biomass changes independent from fishing activity, the actual time is dependent on N0, K.mu, Lmat, repro_wt  and rmax.BH. For the estimation of carrying capacity the first 10 years of the simulation are disregarded and only subsequent years where no fishing took place are used to estimate the annual mean carrying capacity (K). If fishing is simulated for all years or fishing activities start before ten years after simulation start no carrying capacity is estimated.
#'
#' @return a list containing growth parameters and length frequency object
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
#' Thorson, J.T., Jensen, O.P, and Zipkin, E.F. 2014. How variable is recruitment 
#' for exploited marine fishes? A hierarchical model for testing life history theory. 
#' Canadian Journal of Fisheries and Aquatic Sciences 71(7):973-983.
#'
#' @importFrom graphics hist
#' @importFrom stats rlnorm runif weighted.mean
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats qnorm rnorm
#' @importFrom TropFishR VBGF
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' res <- virtualPop()
#' names(res)
#'
#' op <- par(mfcol=c(2,1), mar=c(4,4,1,1))
#' plot(N ~ dates, data=res$pop, t="l")
#' plot(B ~ dates, data=res$pop, t="l", ylab="B, SSB")
#' lines(SSB ~ dates, data=res$pop, t="l", lty=2)
#' par(op)
#'
#' pal <- colorRampPalette(c("grey30",5,7,2), bias=2)
#' with(res$lfqbin, image(x=dates, y=midLengths, z=t(catch), col=pal(100)))
#'
#' ### biased results with single monthly sample
#' inds <- res$inds[[1]]
#' plot(mat ~ L, data = inds, pch=".", cex=5, col=rgb(0,0,0,0.05))
#' fit <- glm(mat ~ L, data = inds, family = binomial(link = "logit"))
#' summary(fit)
#'
#' newdat <- data.frame(L = seq(min(inds$L), max(inds$L), length.out=100))
#' newdat$pmat <- pmat_w(newdat$L, Lmat = 40, wmat=40*0.2)
#' pred <- predict(fit, newdata=newdat, se.fit=TRUE)
#' # Combine the hypothetical data and predicted values
#' newdat <- cbind(newdat, pred)
#' # Calculate confidence intervals
#' std <- qnorm(0.95 / 2 + 0.5)
#' newdat$ymin <- fit$family$linkinv(newdat$fit - std * newdat$se.fit)
#' newdat$ymax <- fit$family$linkinv(newdat$fit + std * newdat$se.fit)
#' newdat$fit <- fit$family$linkinv(newdat$fit)  # Rescale to 0-1
#'
#' plot(mat ~ L, data = inds, pch=".", cex=5, col=rgb(0,0,0,0.05))
#' lines(pmat ~ L, newdat, col=8, lty=3)
#' polygon(
#'   x = c(newdat$L, rev(newdat$L)),
#'   y = c(newdat$ymax, rev(newdat$ymin)),
#'   col = adjustcolor(2, alpha.f = 0.3),
#'   border = adjustcolor(2, alpha.f = 0.3)
#' )
#' lines(fit ~ L, newdat, col=2)
#'
#' lrPerc <- function(alpha, beta, p) (log(p/(1-p))-alpha)/beta
#' ( L50 <- lrPerc(alpha=coef(fit)[1], beta=coef(fit)[2], p=0.5) )
#' lines(x=c(L50,L50,0), y=c(-100,0.5,0.5), lty=2, col=2)
#' text(x=L50, y=0.5, labels = paste0("L50 = ", round(L50,2)), pos=4, col=2 )
#'
#'
#'
#' ### all samples combined
#' inds <- do.call("rbind", res$inds)
#' plot(mat ~ L, data = inds, pch=".", cex=5, col=rgb(0,0,0,0.05))
#' fit <- glm(mat ~ L, data = inds, family = binomial(link = "logit"))
#' summary(fit)
#'
#' newdat <- data.frame(L = seq(min(inds$L), max(inds$L), length.out=100))
#' newdat$pmat <- pmat_w(newdat$L, Lmat = 40, wmat=40*0.2)
#' pred <- predict(fit, newdata=newdat, se.fit=TRUE)
#' # Combine the hypothetical data and predicted values
#' newdat <- cbind(newdat, pred)
#' # Calculate confidence intervals
#' std <- qnorm(0.95 / 2 + 0.5)
#' newdat$ymin <- fit$family$linkinv(newdat$fit - std * newdat$se.fit)
#' newdat$ymax <- fit$family$linkinv(newdat$fit + std * newdat$se.fit)
#' newdat$fit <- fit$family$linkinv(newdat$fit)  # Rescale to 0-1
#'
#' plot(mat ~ L, data = inds, pch=".", cex=5, col=rgb(0,0,0,0.05))
#' lines(pmat ~ L, newdat, col=8, lty=3)
#' polygon(
#'   x = c(newdat$L, rev(newdat$L)),
#'   y = c(newdat$ymax, rev(newdat$ymin)),
#'   col = adjustcolor(2, alpha.f = 0.3),
#'   border = adjustcolor(2, alpha.f = 0.3)
#' )
#' lines(fit ~ L, newdat, col=2)
#'
#' lrPerc <- function(alpha, beta, p) (log(p/(1-p))-alpha)/beta
#' ( L50 <- lrPerc(alpha=coef(fit)[1], beta=coef(fit)[2], p=0.5) )
#' lines(x=c(L50,L50,0), y=c(-100,0.5,0.5), lty=2, col=2)
#' text(x=L50, y=0.5, labels = paste0("L50 = ", round(L50,2)), pos=4, col=2 )
#'
#'
#' }
#'
#'



virtualPop2 <- function(tincr = 1/12,
                        K.mu = 0.5, K.cv = 0.1,
                        Linf.mu = 80, Linf.cv = 0.1,
                        t0 = -0.03,
                        ts = 0, C = 0.85,
                        LWa = 0.01, LWb = 3,   ## for cm in g? then final B results also in g
                        Lmat.f = 0.5*Linf.mu, wmat.f = Lmat.f*0.2,
                        Lmat.m = 0.45*Linf.mu, wmat.m = Lmat.m*0.15,
                        rmaxBH = 1000,
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
                        timemin = 0, timemax = 25, timemin.date = as.Date("1980-01-01"),
                        N0 = 1000,
                        fished_t = seq(17,25,tincr),
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
                        seed = NULL){
  if(is.null(modpath) & length(iteration) > 1)
    stop("must specify path (modpath) to save simulation iterations")
  if(is.null(modpath)){
    iteration <- 1
    seed <- NULL}
  for (iter in iteration) {
    if (is.null(modpath) == FALSE) {
      iterpath <- file.path(modpath, iter)
      dir.create(iterpath, showWarnings = FALSE)
      ## seed values
      if(is.null(seed) || is.null(seed)){
        seed <- floor(runif(1,1,1000))
      }
      seed <- seed + iter
      set.seed(seed)
    }
    
    ## Fishing mortality - effort - catchability
    ## Different fishing scenarios
    fishing.scenario <- function(Fishscen){
      if(Fishscen == "Constant"){
        harvest_rate <- c(rep(0, burnin/tincr), rep(F_scen, length(fished_t) - burnin/tincr))
      }
      
      if(Fishscen == "Fish_down"){
        harvest_rate <- c(rep(0, burnin/tincr), seq(0, F_high, length.out = length(fished_t) - burnin/tincr))
      }
      
      if(Fishscen == "Two_way"){
        up <- seq(0, F_high, length.out = (length(fished_t)/2))
        down <- seq(F_high, F_low, length.out = (length(fished_t)- burnin/tincr - length(up)))
        harvest_rate <- c(rep(0, burnin/tincr), up, down)
      }
      
      if(Fishscen == "Multi_peak"){
        up <- seq(0, F_high, length.out = (length(fished_t)/6))
        down <- seq(F_high, F_low, length.out = (length(fished_t)/6))
        up2 <- seq(F_low, F_high, length.out = (length(fished_t)/6))
        down2 <- seq(F_high, F_low, length.out = length(fished_t)- burnin/tincr - length(up) - length(down) - length(up2))
        harvest_rate <- c(rep(0, burnin/tincr), up, down, up2, down2)
        # if(Fishscen == "Multi-peak"){
        #   up <- seq(0, 0.8, length.out = (length(fished_t)/8))
        #   down <- seq(0.8, 0.218, length.out = (length(fished_t)/8))
        #   up2 <- seq(0.218, 0.627, length.out = (length(fished_t)/8))
        #   down2 <- seq(0.627, 0.263, length.out = (length(fished_t)/8))
        #   up3 <- seq(0.263, 0.554, length.out = (length(fished_t)/8))
        #   down3 <- seq(0.554, 0.181, length.out = length(fished_t) - burnin/tincr - length(c(up, down, up2, down2, up3)))
        #   harvest_rate <- c(rep(0, burnin/tincr), up, down, up2, down2, up3, down3)
      }
      
      return(harvest_rate)
    }
    
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
        qmat <- matrix(rep(qtf[1,], dim(Emat)[1]),ncol=dim(qtf)[2],byrow=TRUE) # qft -> qtf?
      }else{
        qmat <- qtf
      }
    }else if(length(qtf)>1){
      qmat <- as.matrix(qtf)
    }
    if(Fishscen == "None"){
      ## If no harvest_rate provided assuming that effort * catchability = fishing mortality
      if(!is.na(harvest_rate[1]) & !is.nan(harvest_rate[1])){
        if(length(as.numeric(harvest_rate))==1){
          harvest_rate <- rep(harvest_rate, length(fished_t))
        }else{
          harvest_rate <- matrix(rep(harvest_rate, each = length(fished_t)), 
                                 ncol = length(harvest_rate), nrow = length(fished = fished_t))
        }
      }else{
        harvest_rate <- Emat * qmat
      }
    }
    
    harvest_rate <- fishing.scenario(Fishscen)
    
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
               pSel <- logisticSelect(Lt=Lt, L50=L50X, wqs=wqsX)},
             gillnet={
               if(!is.na(fleetNo)){
                 sel_listX <- sel_list[[fleet_No]]
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
    
    # times
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
      inds$Linf <- exp(rnorm(nrow(inds), log(Linf.mu), Linf.cv))
      inds$Winf <- LWa*inds$Linf^LWb
      # inds$K <- 10^(phiprime.mu - 2*log10(inds$Linf)) * rlnorm(nrow(inds), 0, K.cv)
      inds$K <- exp(rnorm(nrow(inds), log(K.mu), K.cv))
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
    
    
    reproduce.inds <- function(inds, save = FALSE){
      ## reproduction can only occur of population contains >1 mature individual
      if(rec_dyn == "BH"){  
        if(repro > 0 & sum(inds$mat) > 0){
          ## calc. SSB
          SSB <- sum(inds$W*inds$mat, na.rm = TRUE)
          n.recruits <- ceiling(srrBH(rmaxBH, betaBH, SSB) * repro)
          ## add noise to recruitment process
          n.recruits <- exp(rnorm(1, log(n.recruits), SigmaR))
          ## save SSB + n.recruits for stock recruitment plot
          if(save) stockRec <<- rbind(stockRec, data.frame(SSB = SSB, recruits = n.recruits, time = tj))
          ## make recruits
          offspring <- make.inds(
            id = seq(lastID+1, length.out= n.recruits)
          )
          # express genes in recruits
          offspring <- express.inds(offspring)
          ##combine all individuals
          inds <- rbind(inds, offspring)
        }
      }
      if(rec_dyn == "Constant"){
        if(repro > 0 & sum(inds$mat) >0){
          betaBH = 1
          ## calculate biomass
          SSB <- sum(inds$W * inds$mat, na.rm = TRUE)
          n.recruits <- ceiling(srrBH(rmaxBH, betaBH, SSB) * repro)
          ## 
          ## add noise to recruitment process
          n.recruits <- exp(rnorm(1, log(n.recruits), SigmaR))
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
      }
      if(rec_dyn == "AR"){
        if(repro > 0 & sum(inds$mat) >0){
          ## calculate biomass
          SSB <- sum(inds$W * inds$mat)
          n.recruits <- ceiling(srrBH(rmaxBH, betaBH, SSB) * repro)
          ## add noise to recruitment process
          RecDev <- rnorm(1, -(SigmaR ^ 2) / 2, SigmaR)
          RecDev_AR <- rho + sqrt(1 - rho ^ 2) * RecDev
          R_t <- exp(RecDev_AR)
          n.recruits <- n.recruits * R_t
          ## save SSB + n.recruits for stock recruitment plot
          if(save) stockRec <<- rbind(stockRec, data.frame(SSB = SSB, recruits = n.recruits))
          ## make recruits
          offspring <- make.inds(
            id = seq(lastID+1, length.out = n.recruits)
          )
          ## express genes in recruits
          offspring <- express.inds(offspring)
          ##combine all individuals
          inds <- rbind(inds, offspring)
        }
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
      inds$FdSurvey <- as.numeric(pSelSurvey * 1)
      
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
    
    spmOpt <- function(x, B0){
      K = x[1]
      r = x[2]
      n = x[3]
      Bthat <- rep(NA, length(resf0$pop$B))
      Bthat[1] <- B0
      for(i in 2:length(Bthat)){
        Bthat[i] <- Bthat[i-1] + ((r / (n - 1)) * Bthat[i-1] * (1 - (Bthat[i-1] / K)^(n-1)))*tincr
      }
      sum((resf0$pop$B - Bthat)^2)
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
    
    # Initial population
    lastID <- 0
    inds <- make.inds(
      id=seq(N0)
    )
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
    colnames(stockRec) <- c("SSB", "recruits", "time")
    
    ## For simulation of unfished population
    ## Initial population
    lastID <- 0
    indsf0 <- make.inds(
      id=seq(N0)
    )
    indsf0 <- express.inds(inds = indsf0)
    ## results object 
    resf0 <- list()
    resf0$pop <- list(
      dates = yeardec2date(date2yeardec(timemin.date) + (timeseq - timemin)),
      N = NaN*timeseq,
      B = NaN*timeseq,
      SSB = NaN*timeseq
    )
    
    ## for total catches per time step (e.g. for correction factor in VPA)
    catches <- vector("numeric",length(timeseq))
    
    ## simulation
    if(progressBar) pb <- txtProgressBar(min=1, max=length(timeseq), style=3)
    
    for(j in seq(timeseq)){
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
      
      # population processes
      inds <- grow.inds(inds)
      inds <- mature.inds(inds)
      inds <- reproduce.inds(inds = inds, save = TRUE)
      inds <- death.inds(inds)
      
      
      ## sample lfq data
      if(lfqSamp){
        if(!is.na(numSamp) & !is.null(numSamp)){
          monthlySamples <- min(sum(inds$Fd),numSamp) ## cannot sample fish than which died
          ## writeLines(noquote("The number of length measurements exceeds number of fish being caught."))
          catches[j] <- sum(inds$Fd)
        }else{
          monthlySamples <- ceiling(sum(inds$Fd)*lfqFrac)
          catches[j] <- sum(inds$Fd)
        }
        tmp <- inds[inds$L > 0,] ## have to have positive length
        if(nrow(tmp) > 0){
          samp <- try(sample(seq(tmp$L), monthlySamples, prob = tmp$F/max(tmp$F)), silent = TRUE)
          if(class(samp) != "try-error"){
            lfq[[j]] <- tmp$L[samp]
            indsSamp[[j]] <- tmp[samp,]
          }
          rm(samp)
        }
      }
      
      ## survey samples indivudals without removing them
      if(lfqSampSurvey){
        tmp <- inds[inds$L > 0,] ## have to have positive length
        if(nrow(tmp) > 0){
          sampSurvey <- try(sample(seq(tmp$L), numSampSurvey, prob = tmp$FSurvey), silent = TRUE)
          if(class(sampSurvey) != "try-error"){
            lfqSurvey[[j]] <- tmp$L[sampSurvey]
            indsSampSurvey[[j]] <- tmp[sampSurvey,]
          }
          rm(sampSurvey)
        }
      }
      
      inds <- remove.inds(inds)
      
      # update results
      res$pop$N[j] <- nrow(inds)
      res$pop$B[j] <- sum(inds$W, na.rm = TRUE)
      res$pop$SSB[j] <- sum(inds$W*inds$mat, na.rm = TRUE)
      
      ## simulate unfished population for K, r and SSB_F=0
      # population processes
      indsf0 <- grow.inds(indsf0)
      indsf0 <- mature.inds(indsf0)
      indsf0 <- reproduce.inds(inds = indsf0, save = FALSE)
      indsf0 <- death.inds(indsf0, f0 = TRUE)
      indsf0 <- remove.inds(indsf0)
      
      # update results
      resf0$pop$N[j] <- nrow(indsf0)
      resf0$pop$B[j] <- sum(indsf0$W, na.rm = TRUE)
      resf0$pop$SSB[j] <- sum(indsf0$W * indsf0$mat, na.rm = TRUE)
      
      ## update progressbar
      if(progressBar) setTxtProgressBar(pb, j)
    }
    if(progressBar) close(pb)
    
    
    ## Estimate carrying capacity
    ## Alternative way build into loop above
    startyear  <- as.POSIXlt(timemin.date)
    startyear$year <- startyear$year + spmYears
    year10 <- as.Date(startyear)
    cutoff <- which.min(abs(yeardec2date( date2yeardec(timemin.date) + (timeseq - timemin)) - year10))
    cc_years <- seq(timeseq)[-(1:cutoff)]
    if(length(cc_years) > 3){
      mod <- lm(resf0$pop$B[cc_years] ~ 1)
      resf0$pop$K <- as.numeric(coefficients(mod))
    }
    if(length(cc_years) > 3){
      mod <- lm(resf0$pop$SSB[cc_years] ~ 1)
      resf0$pop$SSBf0 <- as.numeric(coefficients(mod))
    }
    if(length(cc_years) > 3){
      mod <- lm(res$pop$SSB[cc_years] ~ 1)
      res$pop$SSBf <- as.numeric(coefficients(mod))
    }
    
    ## estimate K, r, n 
    resSPM <- optim(par = c(resf0$pop$K, 1, 2), fn = spmOpt, B = resf0$pop$B, hessian = F, method = "BFGS")
    resSPM <- optim(par = resSPM$par, fn = spmOpt, B = resf0$pop$B, hessian = F, method = "BFGS")
    resSPM <- optim(par = resSPM$par, fn = spmOpt, B = resf0$pop$B, hessian = F, method = "CG")
    ## alternatively SANN but takes longer
    
    
    Kest <- resSPM$par[1]
    rest <- resSPM$par[2]
    nest <- resSPM$par[3]
    gammal <- nest^(nest/(nest-1))/(nest-1)
    m <- rest * Kest / (nest^(nest/(nest-1)))
    ## Deterministic reference levels
    if(nest == 1){
      ## Fox reference levels
      Bdmsy <- Kest/exp(1)
      msyd <- rest * Kest / exp(1)
      Fdmsy <- rest
    }else{
      ## Pella and Tomlinson reference levels
      Bdmsy <- nest^(1/(1-nest)) * Kest
      msyd <- m
      Fdmsy <- m/Bdmsy
    }
    
    ## save parameters
    resf0$pop$K2 <- Kest
    resf0$pop$r <- rest
    resf0$pop$n <- nest
    resf0$pop$m <- m
    resf0$pop$gamma <- gammal
    
    ## Export data -------------------------------------------------------------
    
    
    ## for simulation of population without exploitation, necessary to make the lfq export optional:
    if(any(!is.na(fished_t[1]) & !is.nan(fished_t[1])) & (lfqFrac != 0 & !is.na(lfqFrac) & !is.nan(lfqFrac))){
      ## Trim and Export 'lfq'
      lfq2 <- lfq[which(sapply(lfq, length) > 0)]
      ## binned version of lfq
      dates <- yeardec2date( date2yeardec(timemin.date) + (as.numeric(names(lfq2)) - timemin) )
      Lran <- range(unlist(lfq2))
      Lran[1] <- floor(Lran[1])
      Lran[2] <- (ceiling(Lran[2])%/%bin.size + ceiling(Lran[2])%%bin.size + 1) * bin.size
      bin.breaks <- seq(Lran[1], Lran[2], by=bin.size)
      bin.mids <- bin.breaks[-length(bin.breaks)] + bin.size/2
      res$lfqbin <- list(
        sample.no = seq(bin.mids),
        midLengths = bin.mids,
        dates = dates,
        catch = sapply(lfq2, FUN = function(x){
          hist(x, breaks=bin.breaks, plot = FALSE, include.lowest = TRUE)$counts
        })
      )
    }
    
    
    ## export of survey data
    if(addSurvey & any(!is.na(survey_t[1]) & !is.nan(survey_t[1])) &
       (!is.na(numSampSurvey) & !is.na(numSampSurvey))){
      ## Trim and Export 'lfq'
      lfqS2 <- lfqSurvey[which(sapply(lfqSurvey, length) > 0)]
      ## binned version of lfq
      datesS <- yeardec2date( date2yeardec(timemin.date) + (as.numeric(names(lfqS2)) - timemin) )
      Lran <- range(unlist(lfqS2))
      Lran[1] <- floor(Lran[1])
      Lran[2] <- (ceiling(Lran[2])%/%bin.size + ceiling(Lran[2])%%bin.size + 1) * bin.size
      bin.breaks <- seq(Lran[1], Lran[2], by=bin.size)
      bin.mids <- bin.breaks[-length(bin.breaks)] + bin.size/2
      res$lfqbinSurvey <- list(
        sample.no = seq(bin.mids),
        midLengths = bin.mids,
        dates = datesS,
        catch = sapply(lfqS2, FUN = function(x){
          hist(x, breaks=bin.breaks, plot = FALSE, include.lowest = TRUE)$counts
        })
      )
    }
    
    ## Saving reference levels
    ##--------------------------------------------------------------------------------------------------    
    # fished_dates <- res$lfqbin$dates ## yeardec2date(date2yeardec(timemin.date) + (fished_t - timemin))
    # ## fished_index <- seq(fished_t[1],fished_t[length(fished_t)],1)
    # fished_years <- unique(format(fished_dates, "%Y")) ## format(yeardec2date(date2yeardec(timemin.date) + (fished_index - timemin)),"%Y")
    # res$refLev <- list()
    # temp1 <- aggregate(res$pop$SSB, 
    #                    by = list(years = format(res$pop$dates, "%Y")), 
    #                    FUN = mean, na.rm=TRUE) 
    # temp2 <- aggregate(resf0$pop$SSB, 
    #                    by = list(years = format(res$pop$dates, "%Y")), 
    #                    FUN = mean, na.rm=TRUE)
    # 
    # 
    # 
    # 
    # ## SPR
    # ##--------------------------------------------------------------------------------------------------    
    # res$refLev$years <- fished_years
    # res$refLev$SPR <- temp1$x[temp1$years %in% as.numeric(fished_years)] / 
    #   temp2$x[temp1$years %in% as.numeric(fished_years)] 
    # ## SPM refs
    # res$refLev$Bdmsy <- Bdmsy
    # res$refLev$msyd <- msyd
    # res$refLev$Fdmsy <- Fdmsy
    # ## yearly harvest rate and biomass
    # temp <- aggregate(list(bio = res$pop$B), by=list(years = format(res$pop$dates,"%Y")), mean, na.rm=TRUE)
    # bioYear <- temp$bio[temp$years %in% as.numeric(fished_years)]
    # if(any(!is.na(as.numeric(harvest_rate)))){
    #   temp <- aggregate(list(x = as.numeric(harvest_rate)), 
    #                     by=list(year = floor(fished_t)), mean, na.rm=TRUE)
    #   harvest_rateYear <- temp$x
    # }
    # temp <- unique(format(yeardec2date(date2yeardec(timemin.date) + (fished_t - timemin)),"%Y"))
    # harvest_rateYear <- harvest_rateYear[temp %in% fished_years]
    # 
    # ## Production curve
    # bioPlot <- spmPlot(c(res$pop$B[1], Kest,  rest, nest))
    # Prod <- (rest / (nest - 1)) * bioPlot * (1 - (bioPlot / Kest)^(nest-1))
    # 
    # ## only works if fished_t is not empty meaning when lfq data is collected! fix!!
    # 
    # ## for ypr and LBIs
    # ## aggregate lfq data per year
    # c_sum <- with(res$lfqbin,by(t(catch), format(dates, "%Y"), FUN = colSums))
    # c_list <- lapply(as.list(c_sum), c)
    # midLengths <- res$lfqbin$midLengths
    # datesUni <- unique(as.Date(paste0(format(res$lfqbin$dates, "%Y"), "-01-01")))
    # ## cumulative and percentage cumulative catches
    # cumSum <- lapply(c_list, cumsum)
    # cumSum_perc <- vector("list",length(datesUni))
    # for(i in 1:length(c_list)){
    #   cumSum_perc[[i]] <- cumSum[[i]] / sum(c_list[[i]])
    # }
    # 
    # 
    # ## YPR
    # ##--------------------------------------------------------------------------------------------------
    # yprList <- vector("list", length(c_list))
    # for(i in 1:length(c_list)){  ## loop through years
    #   
    #   lfqi <- structure(list(dates = datesUni[i],
    #                          midLengths = midLengths,
    #                          catch = as.matrix(c_list[[i]])),
    #                     class = "lfq")
    #   
    #   lfqi <- lfqModify(lfqi, vectorise_catch = TRUE)  ## only to remove 0s catches of small fish (for Lr)
    #   
    #   lfqi <- c(lfqi,
    #             Linf = Linf.mu,
    #             K = K.mu,
    #             t0 = 0,
    #             C = C,
    #             ts = ts,
    #             t_anchor = weighted.mean(date2yeardec(as.Date(
    #               paste("2015",which(repro_wt != 0),"15",sep="-"))) %% 1,
    #               ## 2015 only as example year to get the decimial of spawning months
    #               w = repro_wt[which(repro_wt != 0)]),  ## weighted mean of t_anchor
    #             M = M,
    #             a = LWa,
    #             b = LWb)
    #   class(lfqi) <- "lfq"
    #   
    #   if(length(lfqi$midLengths) > 1){
    #     lfqi <- lfqModify(lfqi, bin_size = bin.size)
    #     lfqi <- lfqModify(lfqi, plus_group = "Linf")
    #     plus_group <- TRUE
    #     lfqi$Lr <- min(lfqi$midLengths)
    #     fm <- harvest_rateYear[i]
    #     lfqi$FM <- fm * selfunc(lfqi$midLengths, NA)
    #     resi <- predict_mod(param = lfqi,
    #                         type = "ThompBell",
    #                         FM_change = seq(0,3,0.05),
    #                         plot = FALSE, hide.progressbar = TRUE)
    #     
    #     indi <- names(resi$df_Es) %in% c("F01","Fmax","F05")
    #     tmp <- resi$df_Es[indi]
    #     
    #     if(!"F05" %in% names(resi$df_Es)){
    #       tmp <- unlist(c(tmp, F05 = NA))
    #     }
    #   }else{
    #     tmp = data.frame(F01 = NA, Fmax = NA, F05 = NA)
    #   }
    #   
    #   yprList[[i]] <- tmp
    # }
    # yprRes <- do.call(rbind, yprList)
    # 
    # 
    # 
    # ## LBIs
    # ##--------------------------------------------------------------------------------------------------
    # ## mean length of largest 5% (/Linf >.8)
    # numb <- lapply(c_list, function(x) x[rev(order(midLengths))])    # from largest starting
    # midLengthsRev <- midLengths[rev(order(midLengths))]
    # Lmax5 <- vector('numeric',length(datesUni))
    # for(i in 1:length(datesUni)){
    #   numbcum <- cumsum(numb[[i]]) 
    #   numbcumperc <- round(numbcum / sum(numb[[i]]),5)
    #   numbnum5 <- rep(0, length(numbcumperc))
    #   numbnum5[numbcumperc <= 0.05] <- numb[[i]][numbcumperc <= 0.05]
    #   numbnum5[max(which(numbcumperc <= 0.05),na.rm = TRUE) + 1] <- (0.05 - numbcumperc[max(which(numbcumperc <= 0.05),na.rm = TRUE)]) * sum(numb[[i]])
    #   Lmax5[i] <- sum(numbnum5 * midLengthsRev, na.rm = TRUE) / sum(numbnum5, na.rm = TRUE)
    # }
    # res$refLev$Lmax5 <- round(Lmax5,2)
    # ## 95th percentile (/Linf >.8)
    # L95 <- unlist(lapply(cumSum_perc, function(x) min(midLengths[which(x >= 0.95)],na.rm=TRUE)))
    # res$refLev$L95 <- L95
    # ## Pmega (Lopt + 10%) (>.3)
    # Lopt <- (2 / 3) * Linf.mu
    # Pmega <- vector('numeric', length(datesUni))
    # for(i in 1:length(datesUni)){
    #   Pmega[i] <- sum(c_list[[i]][which(midLengths >= (Lopt + 0.1 * Lopt))], na.rm = TRUE) / 
    #     sum(c_list[[i]], na.rm = TRUE)
    # }
    # res$refLev$Pmega <- Pmega
    # res$refLev$Lopt <- Lopt
    # ## Popt
    # Popt <- vector('numeric', length(datesUni))
    # for(i in 1:length(datesUni)){
    #   Popt[i] <- sum(c_list[[i]][which(midLengths > Lopt*0.9 & midLengths < Lopt*1.1)], na.rm = TRUE) /
    #     sum(c_list[[i]], na.rm = TRUE)
    # }
    # ## Pmat
    # Pmat <- vector('numeric', length(datesUni))
    # for(i in 1:length(datesUni)){
    #   Pmat[i] <- sum(c_list[[i]][which(midLengths > Lmat.f)], na.rm = TRUE) /
    #     sum(c_list[[i]],na.rm = TRUE)
    # }
    # ## Pobj
    # Pobj <- vector('numeric', length(datesUni))
    # for(i in 1:length(datesUni)){
    #   Pobj[i] <- Pmega[i] + Popt[i] + Pmat[i]
    # }
    # res$refLev$Popt <- Popt
    # res$refLev$Pmat <- Pmat
    # res$refLev$Pobj <- Pobj
    # ## 25th percentile of length distribution (/Lmat >1)
    # L25 <- unlist(lapply(cumSum_perc, function(x) min(midLengths[which(x >= 0.25)],na.rm=TRUE)))
    # res$refLev$L25 <- L25
    # Lmat <- mean(c(Lmat.f,Lmat.m))
    # ## Lc (/Lmat >1)  ICES: Lc (Length at first catch = 50% of mode)
    # Lc <- L50   ## check again with gillnet selectivity, then ICES formula
    # res$refLev$Lc <- Lc
    # ## alternatively: 50% of mode
    # modes <- unlist(lapply(c_list, function(x) which.max(x)))   ## deleted midLengths[which...
    # Lcalt <- vector('numeric',length(datesUni))
    # for(i in 1:length(datesUni)){
    #   temp <- as.numeric(cumSum_perc[[i]][modes[i]])
    #   Lcalt[i] <- midLengths[which.min(abs(cumSum_perc[[i]] - 0.5 * temp))]
    # }
    # res$refLev$Lc.alt <- Lcalt
    # ## mean length of individuals > Lc (/Lopt ~1 ; /LF=M >= 1)
    # c_listLC <- lapply(c_list, function(x) x[midLengths >= Lc])
    # midLengthsLC <- midLengths[midLengths >= Lc]
    # Lmean <- vector('numeric', length(datesUni))
    # for(i in 1:length(c_listLC)){
    #   Lmean[i] <- sum(midLengthsLC * c_listLC[[i]], na.rm = TRUE) / sum(c_listLC[[i]], na.rm = TRUE)
    # }
    # LFeM <- 0.75 * Lc + 0.25 * Linf.mu
    # res$refLev$Lmean <- round(Lmean,2)
    # res$refLev$LFeM <- LFeM
    # ## length class with maximum in biomass (/Lopt ~1)
    # 
    # midWeights <- LWa * midLengths ^ LWb
    # bio_list <- lapply(c_list, function(x) x * midWeights)
    # Lmaxy <- unlist(lapply(bio_list, function(x) midLengths[which.max(x)]))
    # res$refLev$Lmaxy <- as.numeric(Lmaxy)
    # 
    # 
    # ## ypr absolute refs
    # res$refLev$YPR <- yprRes
    # 
    ##--------------------------------------------------------------------------------------------------
    # res$refLev$states <- data.frame("Lmax5/Linf" = round(Lmax5/Linf.mu,2),
    #                                 "L95/Linf" = round(L95/Linf.mu,2),
    #                                 "Pmega" = round(Pmega,2),
    #                                 "L25/Lmat" = round(L25/Lmat,2),
    #                                 "Lc/Lmat" = round(Lcalt/Lmat,2),
    #                                 "Lmean/Lopt" = round(Lmean/Lopt,2),
    #                                 "Lmaxy/Lopt" = round(Lmaxy/Lopt,2),
    #                                 "Lmean/LFeM" = round(Lmean/LFeM,2),
    #                                 "SPR" = round(res$refLev$SPR,2),
    #                                 "F/Fmsy" = round(harvest_rateYear/Fdmsy,2),
    #                                 "B/Bmsy" = round(bioYear/Bdmsy,2),
    #                                 "F/F01" = round(harvest_rateYear/yprRes[,1],2),
    #                                 "F/Fmax" = round(harvest_rateYear/yprRes[,2],2),
    #                                 "F/F05" = round(harvest_rateYear/yprRes[,3],2))
    # 
    # res$refLev$statesRefPoint <- data.frame("Lmax5/Linf" = ">0.8",
    #                                         "L95/Linf" = ">0.8",
    #                                         "Pmega" = ">0.3",
    #                                         "L25/Lmat" = ">1",
    #                                         "Lc/Lmat" = ">1",
    #                                         "Lmean/Lopt" = "=1",
    #                                         "Lmaxy/Lopt" = "=1",
    #                                         "Lmean/LFeM" = ">=1",
    #                                         "SPR" = ">0.3",
    #                                         "F/Fmsy" = "<=1",
    #                                         "B/Bmsy" = ">=1",
    #                                         "F/F01" = "<=1",
    #                                         "F/Fmax" = "<=1",
    #                                         "F/F05" = "<=1")
    # 
    # ## estimate LBIs again specific for males and females!
    
    
    
    # individuals
    indsSamp <- indsSamp[which(sapply(indsSamp, length) > 0)]
    res$inds <- indsSamp
    
    
    ## individuals in survey
    indsSampSurvey <- indsSampSurvey[which(sapply(indsSampSurvey, length) > 0)]
    res$indsSurvey <- indsSampSurvey
    
    
    # record mean parameters
    res$growthpars <- list(
      K = K.mu,
      Linf = Linf.mu,
      t0 = t0,
      C = C,
      ts = ts,
      t_anchor = weighted.mean(date2yeardec(as.Date(paste("2015",which(repro_wt != 0),"15",sep="-"))) %% 1,
                               w = repro_wt[which(repro_wt != 0)]),  ## weighted mean of t_anchor
      phiprime = phiprime.mu,
      tmaxrecr = tmaxrecr
    )
    
    ## fisheries dependent information
    ## if fisheries are simulated
    if(any(!is.na(fished_t) & !is.nan(fished_t))){
      res$fisheries <- list(
        fished_t = yeardec2date(date2yeardec(timemin.date) +
                                  (timeseq[timeseq %in% fished_t] - timemin)),
        E = Emat,
        q = qmat,
        F = harvest_rate,
        C = catches
      )
    }
    
    
    
    if(plot){
      layout(matrix(c(1,2,3,4,5,6,7,7), ncol=2, byrow=TRUE), heights=c(4,4,4,1))
      par(mar=c(3,4,3,2))
      ## Numbers
      with(res$pop, plot(dates, N, type='l', lwd=2,
                         xlab="",ylab="Numbers",
                         main = "Population trajectory",
                         ylim=c(0,max(resf0$pop$N,na.rm=TRUE))))
      with(resf0$pop, lines(dates, N, col='dodgerblue2', lwd=2))
      points(res$pop$dates[1], N0,pch=4, lwd=2, col='darkred')
      ## Biomass  + SSB
      with(res$pop, plot(dates, B, type='l', lwd=2,
                         xlab="",ylab="Biomass",
                         main = "Biomass trajectory",
                         ylim=c(0,max(resf0$pop$B,na.rm=TRUE))))
      with(res$pop, lines(dates, SSB, lwd=2, lty=3))
      with(resf0$pop, lines(dates, B, col='dodgerblue2', lwd=2))
      with(resf0$pop, lines(dates, SSB, col='dodgerblue2', lwd=2, lty=3))
      abline(h = resf0$pop$K, col='darkred', lwd=2)
      abline(h = resf0$pop$SSBf0, col='darkred', lwd=2,lty=3)
      ## SPM
      with(resf0$pop, plot(dates, B, type='b', lwd=2,
                           xlab="",ylab="Biomass",col='dodgerblue2',
                           pch=16,
                           main = "Surplus production model",
                           ylim=c(0,max(resf0$pop$B,na.rm=TRUE))))
      abline(h = resf0$pop$K2,lwd=2, lty=3, col = 'darkred')
      with(resf0$pop, lines(dates, bioPlot, lwd=2, lty=1, col = 'darkred'))
      ## Production curve
      plot(bioPlot, Prod, type='l', lwd=2, col = 'dodgerblue2',
           main = "Production curve",
           xlab="Biomass",ylab="Surplus production", ylim = c(0,max(Prod,na.rm=TRUE)*1.1))
      segments(x0 = 0, y0 = msyd, x1 = Bdmsy, y1 = msyd, col = "darkred", lty=3, lwd=2)
      segments(x0 = Bdmsy, y0 = 0, x1 = Bdmsy, y1 = msyd, col = "darkred", lty=3, lwd=2)
      
      ## growth curve
      Lplot <- seq(0,Linf.mu+10,0.1)
      suppressWarnings(agePlot <- TropFishR::VBGF(res$growthpars[1:5], L = Lplot))
      ages <- unlist(lapply(indsSamp, function(x) x[["A"]]))
      lengths <- unlist(lapply(indsSamp, function(x) x[["L"]]))
      plot(ages, lengths, pch=16,
           ylim=c(0,Linf.mu+30),
           main = "VBGF",col='dodgerblue2',
           ylab="Length", xlab = "Age")
      lines(agePlot, Lplot, col='darkred',lwd=2)
      
      ## Stock recruitment relationship
      SSBplot <- seq(0,100,0.1)
      n_recruits <- srrBH(rmaxBH,betaBH,SSBplot)
      n_recruits <- n_recruits * rlnorm(length(n_recruits),0, sdlog = srr.cv)
      plot(SSBplot, n_recruits, type = "l",
           xlab="SSB", ylab ="Recruits",
           main = "Stock recruitment relationship",
           lwd=2, col='dodgerblue2',ylim = c(0,max(n_recruits)*1.2))
      
      ## Legend
      par(mar=c(0,0,0,0))
      plot.new()
      legend("center", legend=c("fishing", "no fishing", "Reference levels"),
             ncol = 3,
             col=c("black",'dodgerblue2', "darkred"), lwd=2, lty=1, bty='n',
             x.intersp = 0.4, seg.len = 0.5,cex=1.1)
      
      on.exit(par(mfrow=c(1,1), mar=c(5,4,4,2)))
    }
    
    
    {
      if (is.null(modpath) == FALSE)
        saveRDS(res, file.path(iterpath, "True.rds"))
      if (is.null(modpath))
        return(res)
      rm(res)
      rm(iterpath)
    }
  }
  if (is.null(modpath) == FALSE)
    return(modpath)
  
  return(res)
  
} # end of function

