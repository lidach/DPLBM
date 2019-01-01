virtualPop_par <- function(
  burnin = 10,
  tincr = 1/12,
  K.mu = 0.5, K.cv = 0.1,
  Linf.mu = 80, Linf.cv = 0.1,
  t0 = -0.03,
  ts = 0, C = 0.85,
  LWa = 0.01, LWb = 3,
  Lmat.f = 0.5*Linf.mu, wmat.f = Lmat.f*0.2,
  Lmat.m = 0.45*Linf.mu, wmat.m = Lmat.m*0.15,
  rec_dyn = "BH",
  rmaxBH = 1000,
  betaBH = 1, srr.cv = 0.1,
  SigmaR = 0.737,
  rho = 0.43,
  repro_wt = c(0,0,0,1,0,0,0,0,0,0,0,0),
  M = 0.7,
  Etf = 500,
  qtf = 0.001,
  Fishscen = "None",
  SigmaF = 0.2,
  F_scen = 0.5,
  F_high = 1,
  F_low = 0.01,
  harvest_rate = NaN,
  gear_types = "trawl",   
  L50 = 0.25*Linf.mu,
  wqs = L50*0.2,
  sel_list = list(mesh_size=100, mesh_size1=60,select_dist="lognormal",select_p1=3, select_p2=0.5),  
  bin.size = 1,
  timemin = 0, timemax = 25, timemin.date = as.Date("1980-01-01"),
  N0 = 1000,
  fished_t = seq(17,25,tincr),
  lfqFrac = 1,
  numSamp = NA,
  addSurvey = FALSE,
  survey_t = seq(17,25,1/6),
  survey_L50 = 0.1 * Linf.mu,
  survey_wqs = survey_L50 * 0.2,
  numSampSurvey = 500,
  progressBar = TRUE,
  plot = FALSE,
  seed = NULL,
  iteration = 1,
  modpath = NULL,
  parallel = TRUE,
  no_cores = detectCores() - 1,
  clusterType = "PSOCK"){
  
  
  
  
  if(parallel){ # Parallel version
    
    
    # for (iter in iteration) {
    if (is.null(modpath) == FALSE) {
      
      ## seed values
      if(is.null(seed) || is.null(seed)){
        seed <- floor(runif(1,1,1000))
      }
      seed <- seed + 1
      set.seed(seed)
      # }
      
      ARGS <- list(
        "burnin",
        "tincr",
        "K.mu", "K.cv",
        "Linf.mu", "Linf.cv",
        "t0",
        "ts", "C",
        "LWa", "LWb",
        "Lmat.f", "wmat.f",
        "Lmat.m", "wmat.m",
        "rec_dyn",
        "rmaxBH",
        "betaBH",
        "SigmaR",
        "rho",
        "repro_wt",
        "M",
        "Etf",
        "qtf",
        "Fishscen",
        "SigmaF",
        "F_scen",
        "F_high",
        "F_low",
        "harvest_rate",
        "gear_types",
        "L50",
        "wqs",
        "sel_list",
        "bin.size",
        "timemin", "timemax", "timemin.date",
        "N0",
        "fished_t",
        "lfqFrac",
        "numSamp",
        "addSurvey",
        "survey_t",
        "survey_L50",
        "survey_wqs",
        "numSampSurvey",
        "progressBar",
        "plot",
        "seed")
      
      parFun <- function(x){
        
        
        # call virtualPop (local directory)
        require(fishdynr)
        require(TropFishR)
        source("~/DPLBM/virtualPop2.R")
        fit <- virtualPop(
          burnin = burnin,
          tincr = tincr,
          K.mu = K.mu, K.cv = K.cv,
          Linf.mu = Linf.mu, Linf.cv = Linf.cv,
          t0 = t0,
          ts = ts, C = C,
          LWa = LWa, LWb = LWb,
          Lmat.f = Lmat.f, wmat.f = wmat.f,
          Lmat.m = Lmat.m, wmat.m = wmat.m,
          rec_dyn = rec_dyn,
          rmaxBH = rmaxBH,
          betaBH = betaBH,
          SigmaR = SigmaR,
          rho = rho,
          repro_wt = repro_wt,
          M = M,
          Etf = Etf,
          qtf = qtf,
          Fishscen = Fishscen,
          SigmaF = SigmaF,
          F_scen = F_scen,
          F_high = F_high,
          F_low = F_low,
          harvest_rate = harvest_rate,
          gear_types = gear_types, 
          L50 = L50,
          wqs = wqs,
          sel_list = sel_list, 
          bin.size = bin.size,
          timemin = timemin, timemax = timemax, timemin.date = timemin.date,
          N0 = N0,
          fished_t = fished_t,
          lfqFrac = lfqFrac,
          numSamp = numSamp,
          addSurvey = addSurvey,
          survey_t = survey_t,
          survey_L50 = survey_L50,
          survey_wqs = survey_wqs,
          numSampSurvey = numSampSurvey,
          progressBar = progressBar,
          plot = plot,
          seed = NULL)
        
        # return result
        
        return(fit)
        
      }
      
      
      cl <- parallel::makeCluster(no_cores, type=clusterType)
      parallel::clusterExport(cl, varlist = ARGS, envir=environment())
      res <- parLapply(cl, iteration, parFun)
      stopCluster(cl)
      
      for(iter in iteration){
        iterpath <- file.path(modpath, iter)
        dir.create(iterpath, showWarnings = FALSE)
        
        if (is.null(modpath) == FALSE)
          saveRDS(res[[iter]], file.path(iterpath, "True.rds"))
        if (is.null(modpath))
          return(res)
      }
    }
    if (is.null(modpath) == FALSE)
      return(modpath)
    
    return(res)
  }
  
  if(!parallel){ # Non-parallel version
    
    if(is.null(modpath) & length(iteration) > 1)
      stop("must specify path (modpath) to save simulation iterations")
    if(is.null(modpath)){
      iteration <- 1
      seed <- NULL}
    
    # call ELEFAN_GA
    require(TropFishR)
    require(fishdynr)
    source("~/DPLBM/virtualPop2.R")
    fit <- virtualPop(
      
      burnin = burnin,
      tincr = tincr,
      K.mu = K.mu, K.cv = K.cv,
      Linf.mu = Linf.mu, Linf.cv = Linf.cv,
      t0 = t0,
      ts = ts, C = C,
      LWa = LWa, LWb = LWb,
      Lmat.f = Lmat.f, wmat.f = wmat.f,
      Lmat.m = Lmat.m, wmat.m = wmat.m,
      rec_dyn = rec_dyn,
      rmaxBH = rmaxBH,
      betaBH = betaBH,
      SigmaR = SigmaR,
      rho = rho,
      repro_wt = repro_wt,
      M = M,
      Etf = Etf,
      qtf = qtf,
      Fishscen = Fishscen,
      SigmaF = SigmaF,
      F_scen = F_scen,
      F_high = F_high,
      F_low = F_low,
      harvest_rate = harvest_rate,
      gear_types = gear_types, 
      L50 = L50,
      wqs = wqs,
      sel_list = sel_list, 
      bin.size = bin.size,
      timemin = timemin, timemax = timemax, timemin.date = timemin.date,
      N0 = N0,
      fished_t = fished_t,
      lfqFrac = lfqFrac,
      numSamp = numSamp,
      addSurvey = addSurvey,
      survey_t = survey_t,
      survey_L50 = survey_L50,
      survey_wqs = survey_wqs,
      numSampSurvey = numSampSurvey,
      progressBar = progressBar,
      plot = plot,
      seed = NULL)
    
    # return result
    return(fit)
  }
  
  
  
} # end of function