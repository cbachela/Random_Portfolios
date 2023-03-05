  
  
  ############################################################################
  ### RANDOM PORTFOLIOS - MSCI COUNTRY INDICES
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     13.01.2023
  # First version:    13.01.2023
  # --------------------------------------------------------------------------
  

  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(volesti)
  require(stringr)
  require(garcholz)
  require(slolz)
  require(simolz)
  require(covolz)
  require(BSS)
  require(BBSolz)
  require(Backtest)
  require(RP)
  wd_data <- "R:/Asset_Management/Research_Projects/External_Research_Projects/GeomScale/Random_Portfolios/"
  wd <- "H:/R/github/Random_Portfolios/"
  source( paste0(wd, "Source/class_BacktestCustom.R") )
  source( paste0(wd, "Source/custom_functions.R") )
  
  
  
  
  # --------------------------------------------------------------------------
  # FUNCTIONS
  # --------------------------------------------------------------------------
  
  simFUN <- function(wmat) 
  { 
    simPortfolio( X = BT0$data$X_sim, 
                  wghts = wmat,
                  fc = 0, vc = 0,
                  language = "C" )
  }
  
  
  # --------------------------------------------------------------------------
  # LOAD MSCI COUNTRY INDEX DATA
  # --------------------------------------------------------------------------
  
  X_msci <- getMSCIData( universe = "dm", frqncy = "d" )
  X_msci <- X_msci[ ,-which(colnames(X_msci) %in% c("IL", "GR"))]  
  X_msci_est <- X_msci[isWeekday(time(X_msci)), ]
  
  X_bm <- getMSCIData( universe = "bm", frqncy = "d" )
  
  
  # --------------------------------------------------------------------------
  # RANDOM PORTFOLIO - DIRICHLET
  # --------------------------------------------------------------------------
  
  n <- ncol(X_msci)
  n_sim <- 10^3
  Names <- paste0("dirichlet", 1:n_sim)
  lambda <- 0.1
  alpha <- rep(1, n) * lambda
  rebdates <- rownames(X_msci_est)[ seq(from = 100, to = nrow(X_msci_est), by = 63) ]
  lWeights <- list()
  lSim <- list()
  for ( j in 1:n_sim ) {
    samples <- rdirichlet( n = length(rebdates), alpha = alpha )
    colnames(samples) <- colnames(X_msci)
    samples <- as.timeSeries( samples, rebdates )
    lWeights[[j]] <- samples
    lSim[[j]] <- simPortfolio( X = X_msci, 
                               wghts = samples, 
                               fc = 0, vc = 0 )
  }
  names(lWeights) <- names(lSim) <- Names
  
  w_array <- list2array( lWeights )
  sim <- do.call( cbind, lSim )
    
  
  lStats <- descStats( sim )
  colors <- fBasics::divPalette(n = ncol(sim), "RdYlGn")
  statsolz:::plot.stats( lStats, sortby = NULL, col = colors )

  
  
  # --------------------------------------------------------------------------
  # SORTED BACKTESTS
  # --------------------------------------------------------------------------
  
  
  # Random sorts
  lStats_rnd <- descStats( sim )
  
  
  # Sort by momentum
  
  # Compute momentum scores
  customPortfolio <- function(GPS)
  {
    X <- getData(GPS)
    mu <- meanGeo(X)
    GPS@solver$portfolio <- "minvariance"
    GPO <- GPO::gpo(GPS)
    GPO@weights <- mu
    return( GPO )
  }
  BT_mom <- BacktestBase$new()
  BT_mom$setCtrl( GPS = GPO::gps( Data = X_msci_est ),
                  width = 252,
                  rebdates = rebdates )
  BT_mom$spec$GPS@solver$portfolio <- "custom"
  BT_mom$data <- list( X_est = X_msci_est,
                       X_sim = X_msci )
  BT_mom$initRBE()
  # debugonce( customPortfolio )
  BT_mom$run()

      
  
  Score_mom <- BT_mom$output$weights
  # debugonce( simSortedBy )
  sim_mom_is <- simSortedBy( sim = sim,
                             BT0 = BT_mom,
                             Score = Score_mom,
                             w_array = w_array,
                             insample = TRUE )
  sim_mom_is_score <- attr(sim_mom_is, "score")
  sim_mom_is_ordering <- attr(sim_mom_is, "ordering")
  sim_mom_oos <- simSortedBy( sim = sim,
                              BT0 = BT_mom,
                              Score = Score_mom,
                              w_array = w_array,
                              insample = FALSE )
  sim_mom_oos_score <- attr(sim_mom_oos, "score")
  sim_mom_oos_ordering <- attr(sim_mom_oos, "ordering")
  lStats_mom_is <- descStats( sim_mom_is )
  stats_mom_is <- lStats_mom_is$stats
  lStats_mom_oos <- descStats( sim_mom_oos )
  stats_mom_oos <- lStats_mom_oos$stats
  
  
  # Sort by variance
  
  # Compute variance scores
  customPortfolio <- function(GPS)
  {
    X <- getData(GPS)
    sds <- apply( X, 2, sd )
    GPS@solver$portfolio <- "minvariance"
    GPO <- GPO::gpo(GPS)
    GPO@weights <- sds
    return( GPO )
  }
  BT_vol <- BacktestBase$new()
  BT_vol$setCtrl( GPS = GPO::gps( Data = X_msci_est ),
                  width = 252,
                  rebdates = rebdates )
  BT_vol$spec$GPS@solver$portfolio <- "custom"
  BT_vol$data <- list( X_est = X_msci_est,
                       X_sim = X_msci )
  BT_vol$initRBE()
  BT_vol$run()
  
  
  Score_vol <- BT_vol$output$weights
  sim_vol_is <- simSortedBy( sim = sim,
                             BT0 = BT_vol,
                             Score = Score_vol,
                             w_array = w_array,
                             insample = TRUE )
  sim_vol_is_score <- attr(sim_vol_is, "score")
  sim_vol_is_ordering <- attr(sim_vol_is, "ordering")
  sim_vol_oos <- simSortedBy( sim = sim,
                              BT0 = BT_vol,
                              Score = Score_vol,
                              w_array = w_array,
                              insample = FALSE )
  sim_vol_oos_score <- attr(sim_vol_oos, "score")
  sim_vol_oos_ordering <- attr(sim_vol_oos, "ordering")
  lStats_vol_is <- descStats( sim_vol_is )
  stats_vol_is <- lStats_vol_is$stats
  lStats_vol_oos <- descStats( sim_vol_oos )
  stats_vol_oos <- lStats_vol_oos$stats
  
  
  idx <- 2
  plot( x = as.numeric(sim_vol_oos_ordering[idx, ]), y = as.numeric(sim_mom_oos_ordering[idx, ]) )
  plot( x = as.numeric(sim_vol_is_ordering[idx, ]), y = as.numeric(sim_mom_is_ordering[idx, ]) )
  
   
  unlist( lapply( 1:nrow(sim_vol_oos_ordering), FUN = function(idx) { cor( as.numeric(sim_vol_oos_ordering[idx, ]), as.numeric(sim_mom_oos_ordering[idx, ]) ) } ) )
  unlist( lapply( 1:(nrow(sim_vol_oos_score)-1), FUN = function(idx) { cor( as.numeric(sim_vol_oos_ordering[idx, ]), as.numeric(sim_mom_oos_ordering[idx+1, ]) ) } ) )
  
  
  
  
  
  
  
  lSim_is <- list( momentum = sim_mom_is,
                   lowvola = sim_vol_is,
                   rnd_sort = sim )
  lSim_oos <- list( momentum = sim_mom_oos,
                    lowvola = sim_vol_oos,
                    rnd_sort = sim )
  lStats_is <- list( momentum = lStats_mom_is,
                     lowvola = lStats_vol_is,
                     rnd_sort = lStats_rnd )
  lStats_oos <- list( momentum = lStats_mom_oos,
                      lowvola = lStats_vol_oos,
                      rnd_sort = lStats_rnd )
  
  
  factor_names <- names(lStats_oos)
  
  par( mfrow = c(3, 2) )
  for ( factor_name in factor_names ) {
    plot( x = 1:ncol(sim),
          y = lStats_is[[ factor_name ]]$stats["means", ], 
          pch = 19, col = colors, main = factor_name )
  }
  
  
  par( mfrow = c(3, 2) )
  for ( factor_name in factor_names ) {
    plot( x = 1:ncol(sim),
          y = lStats_oos[[ factor_name ]]$stats["means", ], 
          pch = 19, col = colors, main = factor_name )
  }
  
  par( mfrow = c(3, 2) )
  for ( factor_name in  factor_names ) {
    plot( x = lStats_is[[ factor_name ]]$stats["sds", ], 
          y = lStats_is[[ factor_name ]]$stats["means", ], 
          pch = 19, col = colors, main = factor_name )
  }
  
  par( mfrow = c(3, 2) )
  for ( factor_name in factor_names ) {
    plot( x = lStats_oos[[ factor_name ]]$stats["sds", ], 
          y = lStats_oos[[ factor_name ]]$stats["means", ], 
          pch = 19, col = colors, main = factor_name )
  }
  
  
  # Correlations
  tmp <- lapply( factor_names, FUN = function(factor_name) {
    cor( x = lStats_is[[ factor_name ]]$stats["sds", ], 
         y = lStats_is[[ factor_name ]]$stats["means", ] )
  } )
  correlations_is <- setNames( unlist( tmp ), factor_names )
  tmp <- lapply( factor_names, FUN = function(factor_name) {
    cor( x = lStats_oos[[ factor_name ]]$stats["sds", ], 
         y = lStats_oos[[ factor_name ]]$stats["means", ] )
  } )
  correlations_oos <- setNames( unlist( tmp ), factor_names )
  
  correlations_is
  correlations_oos
  
  
  
  # --------------------------------------------------------------------------
  # FACTOR CONTROL
  # --------------------------------------------------------------------------

  customPortfolio <- function(GPS)
  {
    X <- getData(GPS)
    sds <- apply( X, 2, sd )
    # Add variance constraint
    addConstraint(GPS) <- linearConstraint( name = "vola", 
                                            sense = "=", 
                                            rhs = mean(sds), 
                                            Amat = matrix(sds, nrow = 1, 
                                                          dimnames = list("vola", 
                                                                          names(sds))) )
    GPS@solver$portfolio <- "minvariance"
    GPO <- GPO::gpo(GPS)
    RBE <- rebalenv( rebdate = as.character(GPS@date) )
    RBE$GPS <- GPS
    RBE$GPO <- GPO
    return( GPO )
  }
  
  BT_Fctrl <- BacktestBase$new()
  BT_Fctrl$setCtrl( GPS = GPO::gps( Data = X_msci_est ),
                    width = 252,
                    rebdates = rebdates )
  BT_Fctrl$data <- list( X_est = X_msci_est,
                         X_sim = X_msci )
  BT_Fctrl$spec$GPS@solver$portfolio <- "custom"
  BT_Fctrl$initRBE()
  # debugonce( customPortfolio )
  BT_Fctrl$run()
  sim_fctrl <- BT_Fctrl$simulate( fc = 0, vc = 0 )
  
  
  
  
  tmp <- na.omit( cbind( X_bm, 
                         sim_fctrl, 
                         apply(sim_mom_oos, 1, mean),
                         apply(sim_vol_oos, 1, mean) ) )
  descStats( tmp ) 
  
  plotSimTS( as.simTS(tmp) )
  
  
  
  
  GPS <- rebalenv( rebdate = "2001-01-26" )$GPS
  Cons <- getConstraints(GPS)
  Cons@linear$Amat
  P <- lincon2Polytope( object = Cons )
  
  idx_eq <- which(P$sense == "=")
  A <- P$A[-idx_eq, ]
  b <- P$b[-idx_eq]
  Aeq <- matrix( P$A[idx_eq, ], nrow = length(idx_eq), ncol = ncol(P$A), byrow = TRUE )
  beq <- P$b[idx_eq]
  
  # Sample with volesti
  # pre_proc_list = preprocess_with_quadprog(A, b, Aeq, beq)
  # debugonce( preprocess_with_quadprog )
  result_list <- samples_uniform_portfolios( A = A, b = b, Aeq = Aeq, beq = beq, ess = 10^3 )
  # samples <- t(result_list$samples)
  samples <- t(result_list$random_portfolios)
  colnames(samples) <- Cons@selection
  samples <- samples[ sample(1:nrow(samples))[1:n_sim], ]
  
  apply( samples, 1, sum )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BEAR BULL ANALYSIS
  # --------------------------------------------------------------------------
  
  
  BBSObj <- BBSRC$new()
  BBSObj$setCtrl()
  BBSObj$data$X_level <- cumulated(X_bm, "discrete")
  factor_names <- names(lSim_oos)
  lPhase_stats <- list()
  for ( factor_name in factor_names ) {
    BBSObj$data$X <- na.omit( cbind( X_bm, 
                                     lSim_oos[[ factor_name ]] ) )
    BBSObj$runRobust()
    lPhase_stats[[ factor_name ]] <- BBSObj$phaseStats()
  }
  
  
  par( mfrow = c(2, 1) )
  colors <- fBasics::divPalette(n = ncol(sim), "RdYlGn")
  factor_name <- "momentum"
  bm_names <- colnames(X_bm)
  plot( x = lPhase_stats[[ factor_name ]]$`states==1`["sds", ], 
        y = lPhase_stats[[ factor_name ]]$`states==1`["means", ], pch = 19, col = colors )
  for ( i in seq(along = bm_names) ) {
    points( x = lPhase_stats[[ factor_name ]]$`states==1`["sds", bm_names[i]], 
            y = lPhase_stats[[ factor_name ]]$`states==1`["means", bm_names[i]], 
            pch = 19, cex = 2, col = i )
  }
  plot( x = lPhase_stats[[ factor_name ]]$`states==-1`["sds", ], 
        y = lPhase_stats[[ factor_name ]]$`states==-1`["means", ], pch = 19, col = colors )
  for ( i in seq(along = bm_names) ) {
    points( x = lPhase_stats[[ factor_name ]]$`states==-1`["sds", bm_names[i]], 
            y = lPhase_stats[[ factor_name ]]$`states==-1`["means", bm_names[i]], 
            pch = 19, cex = 2, col = i )
  }

  
  
  
  
  stats_field <- "sharpe"
  ldens_bull <- lapply( lPhase_stats, FUN = function(x) { density(x$`states==1`[stats_field, ]) })
  ldens_bear <- lapply( lPhase_stats, FUN = function(x) { density(x$`states==-1`[stats_field, ]) })
  par( mfrow = c(2, 1) )
  slolz:::plot.ldensity( ldens_bull, fillin = TRUE, main = "Bull" )
  slolz:::plot.ldensity( ldens_bear, fillin = TRUE, main = "Bear" )
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # NO REBALANCING, I.E., PURELY IN-SAMPLE ANALYSIS
  # --------------------------------------------------------------------------
  
  mu <- meanGeo(X_msci)
  covmat <- cov(X_msci)
  
  mu_rp <- as.numeric( samples %*% mu )
  var_rp <- as.numeric( apply( samples, 1, function(w) { t(w) %*% covmat %*% w } ) )
  
  colors <- fBasics::divPalette(n = length(mu_rp), "RdYlGn")
  plot( x = var_rp, y = mu_rp, pch = 19, col = colors )
  
  
  # Sort by momentum
  ordering <- order( mu_rp )
  plot( x = var_rp[ordering], y = mu_rp[ordering], pch = 19, col = colors )

  
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
