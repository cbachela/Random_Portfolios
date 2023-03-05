
  
  ############################################################################
  ### RANDOM PORTFOLIOS - RANDOM PORTOLIOS WITHIN 5% UPPER BOUNDS - USA (MOMETNUM)
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     13.11.2022
  # First version:    13.11.2022
  # --------------------------------------------------------------------------
  
  
  # We sample random portfolios within a 5% upper bound which is the only constraint
  # used in the MSCI Momentum strategy.
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------

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
  # PARAMETERS
  # --------------------------------------------------------------------------
  
  universe <- "usa"
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  
  
  
  # --------------------------------------------------------------------------
  # Initialize Backtest object
  # --------------------------------------------------------------------------
  
  
  BT0 <- loadBacktest( wd = paste0(wd_data, "Data/"),
                       name = "usa_maxw5" )
  BT0$initRBE()
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BENCHMARK SERIES
  # --------------------------------------------------------------------------
  
  tickers <- setNames( c("NDDUUS", "M1USMMT", "M1USEV", "M1US000V", "M1USQU", "M00IMVSO", "M1USSZT", "M1WODV"),
                       c("Capw", "Momentum", "Enhanced_Value", "Value", "Quality", "Min_Var", "Size", "Multifactor") )
  X_fact <- rodbcGetOLZDBReturns( assetName = tickers,
                                  refCcy = "USD",
                                  frqncy = "daily",
                                  na.rm = "r" )
  X_fact <- X_fact[isWeekday(time(X_fact)), tickers]
  colnames(X_fact) <- names(tickers)
  
 
  
  # --------------------------------------------------------------------------
  # BACKTEST QUASI CAPW, I.E., BM WITHIN OLZ CONSTRAINTS
  # --------------------------------------------------------------------------
  
  BacktestCustom.customPortfolio <- function( rebdate = NULL )
  {
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    .self$benchmarkWeights( rebdate = rebdate )
    w_bm <- getBenchmarkWeights( RBE$GPS )
    # w_bm <- w_bm[selection] / sum(w_bm[selection])  # Not necessary, same result without the rescaling
    setInitialWeights(RBE$GPS) <- w_bm
    GPO <- minturnoverPortfolio( GPS = RBE$GPS )
    w <- getWeights(GPO) 
    wghts <- timeSeries( matrix( w, nrow = 1, dimnames = list(NULL, names(w)) ),
                         rebdate )
    .self$appendOutput( list( weights = wghts ) )
    return( TRUE )
  }
  BacktestCustom$methods( customPortfolio = BacktestCustom.customPortfolio )
  
  
  BT_qcapw <- BacktestCustom$new()
  BT_qcapw$data <- BT0$data
  BT_qcapw$spec <- BT0$spec
  BT_qcapw$spec$portfolio <- "custom"
  BT_qcapw$spec$GPS@covariance$method <- "duv"
  BT_qcapw$spec$selection_filter <- c("db_flag")
  BT_qcapw$spec$cons_fun <- "box"
  BT_qcapw$runLoop()
  
  sim_qcapw <- BT_qcapw$simulate( fc = 0, vc = 0 )
  sim_qcapw <- sim_qcapw[isWeekday(time(sim_qcapw)), ]
  
  
  
  
  # --------------------------------------------------------------------------
  # BACKTEST QUASI EQW PORTFOLIO
  # --------------------------------------------------------------------------
  
  BT_qeqw <- BacktestCustom$new()
  BT_qeqw$data <- BT0$data
  BT_qeqw$spec <- BT0$spec
  BT_qeqw$spec$portfolio <- NULL
  BT_qeqw$spec$GPS@covariance$method <- "duv"
  BT_qeqw$spec$selection_filter <- "db_flag"
  BT_qeqw$spec$cons_fun <- "box"
  # BT_qeqw$spec$selection_filter <- "bm"
  # BT_qeqw$spec$cons_fun <- NULL
  BT_qeqw$runLoop()
  sim_qeqw <- BT_qeqw$simulate( fc = 0, vc = 0 )
  sim_qeqw <- sim_qeqw[isWeekday(time(sim_qeqw)), ]
  
  
  # Add qeqw and qcapw to X_bm
  X_bm <- na.omit( cbind( X_fact, QEQW = sim_qeqw, QCapw = sim_qcapw ) )
  
  
  
  # --------------------------------------------------------------------------
  # RANDOM PORTFOLIOS - UNIFORM
  # --------------------------------------------------------------------------
  
  BT2 <- BacktestCustom$new()
  BT2$data <- BT0$data
  BT2$spec <- BT0$spec
  BT2$spec$portfolio <- "momentum_rp"
  BT2$spec$n_sim <- 10^2 * 5
  BT2$spec$sampling_dist = "uniform"
  # BT2$spec$m <- NULL # means m = length(selection)
  BT2$spec$m <- 30
  BT2$spec$keep_all_rebdates <- TRUE
  BT2$spec$shadow_dirichlet_transform <- FALSE
  BT2$spec$scl_by_capw <- FALSE
  BT2$runLoop()  
  
  # Simulate
  tmp <- lapply( BT2$output, FUN = simFUN )
  sim_2 <- do.call( cbind, tmp )
  sim_2 <- sim_2[isWeekday(time(sim_2)), ]
  colnames(sim_2) <- paste0("rp_uniform_", 1:ncol(sim_2))
  
  
  # Project weights to boundary of feasible set
  lWeights_eqw <- BT2$output
  for ( i in 1:length(lWeights_eqw) ) {
    for ( today in rownames(lWeights_eqw[[i]]) ) {
      RBE <- rebalenv( rebdate = today )
      # Cons <- getConstraints(RBE$GPS)
      lincon <- getConstraints( RBE$GPS, "bounds" )
      upper_bounds <- lincon$upper
      w_unc <- lWeights_eqw[[i]][today, ]
      idx_notNA <- which(!is.na(w_unc))
      w_unc <- w_unc[ ,idx_notNA]
      w_unc <- setNames( as.numeric(w_unc), names(w_unc))
      # debugonce( map2boxcon )
      lWeights_eqw[[i]][today, idx_notNA] <- map2boxcon( w_unc = w_unc, upper = upper_bounds, itermax = 10^3 )
    }
  }
  
  # Simulate
  tmp <- lapply( lWeights_eqw, FUN = simFUN )
  sim_2_cons <- do.call( cbind, tmp )
  sim_2_cons <- sim_2_cons[isWeekday(time(sim_2_cons)), ]
  colnames(sim_2_cons) <- paste0("rp_uniform_cons", 1:ncol(sim_2_cons))
  
  
  
  # --------------------------------------------------------------------------
  # RANDOM PORTFOLIOS - CAPW
  # --------------------------------------------------------------------------
  
  
  BT3 <- BacktestCustom$new()
  BT3$data <- BT0$data
  BT3$spec <- BT0$spec
  BT3$spec$portfolio <- "momentum_rp"
  BT3$spec$n_sim <- 10^2 * 5
  # BT3$spec$sampling_dist = "uniform"
  BT3$spec$sampling_dist = "capw"
  # BT3$spec$m <- NULL # means m = length(selection)
  BT3$spec$m <- 30
  # BT3$spec$th <- 0.1
  BT3$spec$keep_all_rebdates <- TRUE
  BT3$spec$shadow_dirichlet_transform <- FALSE
  BT3$spec$scl_by_capw <- FALSE
  # BT3$spec$n_sim <- 10^3
  # debugonce( BT3$momentum_rpPortfolio )
  # debugonce( BT3$gps )
  # debugonce( BT3$selection.db_flag )
  BT3$runLoop()  
  
  # Simulate
  tmp <- lapply( BT3$output, FUN = simFUN )
  sim_3 <- do.call( cbind, tmp )
  sim_3 <- sim_3[isWeekday(time(sim_3)), ]
  colnames(sim_3) <- paste0("rp_capw", 1:ncol(sim_3))
  
  
  # Project weights to boundary of feasible set
  lWeights_capw <- BT3$output
  for ( i in 1:length(lWeights_capw) ) {
    for ( today in rownames(lWeights_capw[[i]]) ) {
      RBE <- rebalenv( rebdate = today )
      # Cons <- getConstraints(RBE$GPS)
      lincon <- getConstraints( RBE$GPS, "bounds" )
      upper_bounds <- lincon$upper
      w_unc <- lWeights_capw[[i]][today, ]
      idx_notNA <- which(!is.na(w_unc))
      w_unc <- w_unc[ ,idx_notNA]
      w_unc <- setNames( as.numeric(w_unc), names(w_unc))
      # debugonce( map2boxcon )
      lWeights_capw[[i]][today, idx_notNA] <- map2boxcon( w_unc = w_unc, upper = upper_bounds, itermax = 10^3 )
    }
  }
  
  # Simulate
  tmp <- lapply( lWeights_capw, FUN = simFUN )
  sim_3_cons <- do.call( cbind, tmp )
  sim_3_cons <- sim_3_cons[isWeekday(time(sim_3_cons)), ]
  colnames(sim_3_cons) <- paste0("rp_capw_cons", 1:ncol(sim_3_cons))
  
  
  
  
  
  
  
  do.call( cbind, lapply( BT3$output, FUN = function(X) { apply( X, 1, max, na.rm = TRUE ) } ) )
  do.call( cbind, lapply( BT3$output, FUN = function(X) { apply( X, 1, FUN = function(x) { length(na.omit(x) > 0)} ) } ) )
  tmp <- do.call( cbind, lapply( BT3$output, FUN = function(X) { apply( X, 1, FUN = function(x) { sum(na.omit(x))} ) } ) )
  tmp[is.na(tmp)] <- 0
  plot(as.numeric(tmp))
  
  tmp <- do.call( cbind, lapply( lWeights_capw, FUN = function(X) { apply( X, 1, max, na.rm = TRUE ) } ) )
  tmp[is.na(tmp)] <- 0
  tmp
  which(tmp > 0.11, arr.ind = TRUE )
  do.call( cbind, lapply( lWeights_capw, FUN = function(X) { apply( X, 1, FUN = function(x) { length(na.omit(x) > 0)} ) } ) )
  tmp <- do.call( cbind, lapply( lWeights_capw, FUN = function(X) { apply( X, 1, FUN = function(x) { sum(na.omit(x))} ) } ) )
  tmp[is.na(tmp)] <- 0
  plot(as.numeric(tmp))
  
  
  weightsBarPlot( lWeights_capw[[1]] )
  
  
  
  
  
  
  
  # Combine
  end_date <- "2022-03-06"
  strategy_names <- colnames(X_bm)
  sim_all <- na.omit( cbind( X_bm[ ,strategy_names], sim_3_cons) )
  sim_all <- window( sim_all, start(sim_all), end_date )
  lStats_all <-  descStats( sim_all, descStatsSpec(annualize = TRUE) )
  
  
  
  
  # sim_all_2 <- na.omit( cbind( X_bm, sim_2_cons) )
  # lStats_2 <- descStats( sim_all_2, descStatsSpec(annualize = TRUE) )
  
  # First decade
  sim_all_decade1 <- window( sim_all, start(sim_all), "2011-01-01" )
  lStats_decade1 <-  descStats( sim_all_decade1, descStatsSpec(annualize = TRUE) )
  
  # Second decade
  sim_all_decade2 <- window( sim_all, "2011-01-01", end(sim_all) )
  lStats_decade2 <-  descStats( sim_all_decade2, descStatsSpec(annualize = TRUE) )
  
  # Bear Bull 
  BBSObj <- BBSRC$new()
  BBSObj$setCtrl()
  BBSObj$data$X_level <- cumulated( BT0$data$X_bm, "discrete" )
  BBSObj$runRobust()
  BBSObj$data$X <- sim_all
  phase_stats <- BBSObj$phaseStats()
  
  
  
  
  
  
  
  ###
  # Risk return plot
  
  lStats <- lStats_all
  df <- data.frame( x = lStats$stats["sds", ],
                    y = lStats$stats["means", ] )
  # df <- data.frame( x = lStats_decade1$stats["sds", ],
  #                   y = lStats_decade1$stats["means", ] )
  # df <- data.frame( x = lStats_decade2$stats["sds", ],
  #                   y = lStats_decade2$stats["means", ] )
  # df <- data.frame( x = phase_stats$'states==1'["sds", ],
  #                   y = phase_stats$'states==1'["means", ] )
  # df <- data.frame( x = phase_stats$'states==-1'["sds", ],
  #                   y = phase_stats$'states==-1'["means", ] )
  h1 <- hist(df$x, breaks=150, plot=F)
  h2 <- hist(df$y, breaks=150, plot=F)
  top <- max(h1$counts, h2$counts)
  k <- kde2d( df$x, df$y, n = 100 )
  # margins
  oldpar <- par()
  par(mar=c(3,3,1,1))
  layout(matrix(c(2,0,1,3),2,2,byrow=TRUE),c(3,1), c(1,3))
  image(k)
  # abline( h = df$y[1], col = "white", lwd = 3 )
  # abline( v = df$x[1], col = "white", lwd = 3 )
  # abline( h = df$y[2], col = "white", lwd = 3 )
  # abline( v = df$x[2], col = "white", lwd = 3 )
  points( x = df$x, y = df$y, pch = 19, col = "orange", cex = 0.5 )
  tmp <- MASS::kde2d( x = df$x, y = df$y, n = 25, h = c(0.005, 0.02), lims = c(range(df$x), range(df$y)) ) 
  # contour( tmp, xlab = "Portfolio Standard Deviation", ylab = "Portfolio Return", 
  #          add = TRUE, drawlabels = FALSE, col = "red" )
  colors <- fBasics::rainbowPalette(n = ncol(X_bm))
  colors[1:2] <- c("black", "red")
  for ( i in 1:ncol(X_bm) ) {
    points( x = df$x[i], y = df$y[i], pch = 19, cex = 2, col = colors[i] )
  }
  legend( "topright", colnames(X_bm), lwd = 2, col = colors, text.col = colors, bty = "n")
  par(mar=c(0,2,1,0))
  barplot(h1$counts, axes=FALSE, ylim=c(0, top), space=0, col='grey')
  par(mar=c(2,0,0.5,1))
  barplot(h2$counts, axes=FALSE, xlim=c(0, top), space=0, col='grey', horiz=TRUE)
  
  ###
  
  
  
  
  
  
  
  
  
  
  
  # Rolling outperformance
  
  stats_roll <- applyRoll( Data = sim_all, 
                           Width = 252 * 3,
                           By = 1,
                           FUN = meanGeo )
  
  plot( stats_roll, plot.type = "single" )
  abline( h = 0, col = "grey" )
  
  idx_mom <- which(colnames(sim_all) == "Momentum" )
  plot( stats_roll[ ,-idx_mom] - stats_roll[ ,rep(idx_mom, ncol(stats_roll)-1)], plot.type = "single" )
  abline( h = 0, col = "grey" )
  
  
  
  
  
  # Relative performance
  FUN <- function(i) { simOutperformance( x = sim_all[ ,strategy_name], 
                                          y = sim_all[ ,i] ) }
  lSim_delta <- list()
  for ( strategy_name in colnames(X_bm) ) {
    lSim_delta_tmp <- lapply( (ncol(X_bm)+1):ncol(sim_all), FUN = FUN )
    lSim_delta[[strategy_name]] <- do.call( cbind, lSim_delta_tmp )
  }
  names(lSim_delta)
  
  # plot( log(cumulated(sim_delta, "discrete")), plot.type = "single", col = "orange" )
  # abline( h = 0, col = "grey" )
  
 
  
  
  # Rolling relative performance
  sim_delta <- lSim_delta[["Capw"]]
  stats_roll <- applyRoll( Data = cumulated(sim_delta, "discrete"),
                           Width = 252 * 1,
                           By = 1,
                           FUN = function(x) { as.numeric(x[nrow(x), ]) / as.numeric(x[1, ]) - 1 } )
  
  plot( stats_roll, plot.type = "single", col = "orange" )                                                    
  abline( h = 0, col = "grey" )
  
                                                                                                                                                                                                                           
  perc_outperf <- apply( stats_roll, 1, function(x) { sum(x > 0) / length(x) } )
  plot( perc_outperf )
  abline(h = 0.5)
  
  
  
  # Market regimes
   
  BBSObj <- BBSRC$new()
  BBSObj$setCtrl()
  BBSObj$data$X_level <- cumulated(X_bm[ ,"Capw"], "discrete")
  BBSObj$runRobust()
  states <- BBSObj$output$states
  lPhase_stats <- list()
  for ( strategy_name in names(lSim_delta) ) {
    BBSObj$data$X <- lSim_delta[[strategy_name]]
    lPhase_stats[[strategy_name]] <- BBSObj$phaseStats()
  }
  
  strategy_name <- "Momentum"
  phase_stats <- lPhase_stats[[strategy_name]]
  boxplot( as.data.frame( cbind( bull = phase_stats$`states==1`["means", ],
                                 bear = phase_stats$`states==-1`["means", ] ) ), 
           beside = TRUE, col = c("green", "red") )
  abline( h = 0, col = "grey" )
  
  
  tmp <- lapply( lPhase_stats, function(x) { cbind( bull = x$`states==1`["means", ],
                                                    bear = x$`states==-1`["means", ] ) } )
  tmp <- lapply( names(tmp), function(name) { x <- tmp[[name]]; colnames(x) <- paste0(name, "_", colnames(x)); return(x) } )
  tmp <- do.call( cbind, tmp )
  
  boxplot( as.data.frame( tmp[ ,1:16] ), beside = TRUE, col = c("green", "red") )
  abline( h = 0 )
  
  
  
  statistic <- "means"
  strategy_name <- "Capw"
  phase_stats <- lPhase_stats[[strategy_name]]
  ldens <- lapply( phase_stats, FUN = function(x) { density(x[statistic, ]) } )
  slolz:::plot.ldensity( ldens )
  abline( v = 0, col = "grey" )
  colors <- fBasics::divPalette( n = length(ldens), "RdYlGn" )
  abline( v = mean(phase_stats[[1]][statistic, ]), col = colors[1] )
  abline( v = mean(phase_stats[[2]][statistic, ]), col = colors[2] )
  
  
  
  
  
  
  # Market regimes
  
  BBSObj <- BBSRC$new()
  BBSObj$setCtrl()
  BBSObj$data$X_level <- cumulated(X_bm[ ,"Capw"], "discrete")
  BBSObj$data$X <- X_fact
  BBSObj$runRobust()
  phase_stats_bm <- BBSObj$phaseStats()
  
  
  strategy_name <- c("Capw", "Momentum")
  tmp <- do.call( cbind, lapply( phase_stats_bm, FUN = function(x) { x["means", strategy_name] } ) )
  barplot( tmp, beside = TRUE, col = c("green", "red") )
  
  
  sim_delta_mom_vs_capw <- simOutperformance( x = X_bm[ ,"Momentum"],
                                              y = X_bm[ ,"Capw"] )
  BBSObj <- BBSRC$new()
  BBSObj$setCtrl()
  BBSObj$data$X_level <- cumulated(X_bm[ ,"Capw"], "discrete")
  BBSObj$data$X <- sim_delta_mom_vs_capw
  BBSObj$runRobust()
  phase_stats_bm <- BBSObj$phaseStats()
  tmp <- do.call( cbind, lapply( phase_stats_bm, FUN = function(x) { x["means", ] } ) )
  barplot( tmp, beside = TRUE, col = c("green", "red") )
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # ALPHA
  # --------------------------------------------------------------------------
  
  path_ff <- "R:/Asset_Management/R/myRData/Factor/"
  env_ff <- readRDS( file = paste0(path_ff, "data.Rds") )
  FF3 <- env_ff$lFactor$`3F`$USA_daily
  FF5 <- env_ff$lFactor$`5F`$USA_daily
  
  ls(env_ff)  
  ls(env_ff$lFactor)

 

  
  # Alpha's of random portfolios
  
  
  dates <- intersect( rownames(X_bm), rownames(FF5) )
  # dates <- intersect( dates, rownames(BBSObj$output$states)[ BBSObj$output$states == -1 ] )
  # dates <- intersect( dates, rownames(BBSObj$output$states)[ BBSObj$output$states == 1 ] )
  Y <- sim_3_cons[dates, ]
  # Y <- sim_2_cons[dates, ]
  X_train <- FF5[dates, 1:5]
  lReg <- list()
  for ( j in 1:ncol(Y) ) {
    reg <- regression( Y_train = Y[ ,j], 
                       X_train = X_train, 
                       type = "ols" )
    lReg[[j]] <- reg$coeffmat
  }  
  alpha_score <- unlist( lapply( lReg, FUN = function(x) { x[1, 1] } ) )  
  t_score <- unlist( lapply( lReg, FUN = function(x) { x[1, 3] } ) )  
  p_score <- unlist( lapply( lReg, FUN = function(x) { x[1, 4] } ) )
  
  
  # FF Regression
  strategy_name <- "Min_Var"
  dates <- intersect( rownames(X_bm), rownames(FF5) )
  reg <- regression( Y_train = X_bm[dates, strategy_name],
                     # X_train = FF3[dates, 1:3],
                     X_train = FF5[dates, 1:5],
                     type = "ols" )
  summary(reg$reg)
  mean( alpha_score); sd( alpha_score )
  
  
  
  
  plot(density(alpha_score))
  abline( v = reg$coeffmat[1, 1] )
  
  plot(density(t_score))  
  lines(density(alpha_score / sd(alpha_score)))
  abline( v = 2, col = "grey" )
  abline( v = reg$coeffmat[1, 3] )

  plot(density(p_score))  
  abline( v = 0.05, col = "grey" )
  abline( v = reg$coeffmat[1, 4] )
  
  
  sum( alpha_score > reg$coeffmat[1, 1] ) / length(alpha_score)
  sum( alpha_score / sd(alpha_score) > reg$coeffmat[1, 3] ) / length(alpha_score)
  
  
  descStats( X_bm[ ,1:6] )
  dd <- drawDownStats( X_bm[ ,1:2] )
  dd
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # SORT RP'S BY FACTOR SCORE, CHECK RANK VS. PERFORMANCE
  # --------------------------------------------------------------------------
  
  # debugonce( list2array )
  w_array <- list2array( L = BT3$output )
  sim <- sim_3

  w_array <- list2array( L = lWeights_capw )
  sim <- sim_3_cons
  
  w_array <- list2array( L = BT2$output )
  sim <- sim_2
  
  w_array <- list2array( L = lWeights_eqw )
  sim <- sim_2_cons
  
  
  # Random concatenation
  lStats_rnd <- descStats( sim )
  
  
  # Letter score
  env <- readRDS( file = paste0(wd_data, "waRehouse/alphabet_strategy.rds") )
  lScore <- env$lScore

  Score_Z <- lScore$Z
  # debugonce( simSortedBy )
  sim_z_is <- simSortedBy( sim = sim,
                           BT0 = BT0,
                           Score = Score_Z,
                           w_array = w_array,
                           insample = TRUE )
  sim_z_is_score <- attr(sim_z_is, "score")
  sim_z_oos <- simSortedBy( sim = sim,
                            BT0 = BT0,
                            Score = Score_Z,
                            w_array = w_array,
                            insample = FALSE )
  sim_z_oos_score <- attr(sim_z_oos, "score")
  lStats_z_is <- descStats( sim_z_is )
  stats_z_is <- lStats_z_is$stats
  lStats_z_oos <- descStats( sim_z_oos )
  stats_z_oos <- lStats_z_oos$stats
  
  
  # Compare scores of RP's to score of max exposure portfolio
  z_strat_letterscore <- portfolioScore( wmat = env$lWeights$Z, 
                                         score_mat = Score_Z, 
                                         b_rescale = TRUE, 
                                         agg_method = "arithmetic" )
  z_strat_letterscore # way above the range of the RP's
  plot( sort(as.numeric(sim_z_is_score)) )
  
  
  # Compute Momentum scores
  # BT_mom <- BT3$copy()
  # BT_mom$output <- list()
  BT_mom <- BacktestCustom$new()
  BT_mom$data <- BT0$data
  BT_mom$spec <- BT0$spec
  # BT_mom$spec$width <- 13
  BT_mom$spec$width <- 52
  BT_mom$spec$portfolio <- "momScore"
  BT_mom$runLoop()
  
  Score_mom <- BT_mom$output$scores
  sim_mom_is <- simSortedBy( sim = sim,
                             BT0 = BT0,
                             Score = Score_mom,
                             w_array = w_array,
                             insample = TRUE )
  sim_mom_is_score <- attr(sim_mom_is, "score")
  sim_mom_oos <- simSortedBy( sim = sim,
                              BT0 = BT0,
                              Score = Score_mom,
                              w_array = w_array,
                              insample = FALSE )
  sim_mom_oos_score <- attr(sim_mom_oos, "score")
  lStats_mom_is <- descStats( sim_mom_is )
  stats_mom_is <- lStats_mom_is$stats
  lStats_mom_oos <- descStats( sim_mom_oos )
  stats_mom_oos <- lStats_mom_oos$stats
  
  # Input of Gianluca: are the maximal mom scores of the RP's far from the max mom portfolio scores?
  
  # Maximum momentum portfolio
  BT_maxmom <- BacktestCustom$new()
  BT_maxmom$spec <- BT_mom$spec
  BT_maxmom$spec$selection_filter <- c("db_flag", "NA")
  BT_maxmom$data <- BT_mom$data
  BT_maxmom$spec$portfolio <- "maxMom"
  # debugonce( BT_maxmom$maxMomPortfolio )
  BT_maxmom$runLoop()
  
 
  # Compare scores of RP's to score of max exposure portfolio
  max_mom_score <- portfolioScore( wmat = BT_maxmom$output$weights, 
                                   score_mat = Score_mom, 
                                   b_rescale = TRUE, 
                                   agg_method = "arithmetic" )
  max_mom_score # way above the range of the RP's
  
  plot( sort(as.numeric(sim_mom_oos_score)), ylim = range(na.omit(c(sim_mom_oos_score, max_mom_score))) )
  points( sort(as.numeric(max_mom_score)), col = 2 )
  range(na.omit(sim_mom_oos_score))
  range(max_mom_score)
  
  
  
  
  # Size score
  Score_size <- BT0$data$wmat_bm
  sim_size_is <- simSortedBy( sim = sim,
                              BT0 = BT0,
                              Score = Score_size,
                              w_array = w_array,
                              insample = TRUE )
  sim_size_is_score <- attr(sim_size_is, "score")
  sim_size_oos <- simSortedBy( sim = sim,
                               BT0 = BT0,
                               Score = Score_size,
                               w_array = w_array,
                               insample = FALSE )
  sim_size_oos_score <- attr(sim_size_oos, "score")
  sim_mom_is_score <- attr(sim_mom_is, "score")
  lStats_size_is <- descStats( sim_size_is )
  stats_size_is <- lStats_size_is$stats
  lStats_size_oos <- descStats( sim_size_oos )
  stats_size_oos <- lStats_size_oos$stats
  
  
  
  # Compute Vola scores
  BT_vol <- BacktestCustom$new()
  BT_vol$data <- BT0$data
  BT_vol$spec <- BT0$spec
  BT_vol$spec$width <- 13
  BT_vol$spec$portfolio <- "volScore"
  BT_vol$runLoop()
  
  Score_vol <- BT_vol$output$scores
  sim_vol_is <- simSortedBy( sim = sim,
                                BT0 = BT0,
                                Score = Score_vol,
                                w_array = w_array,
                                insample = TRUE )
  sim_vol_is_score <- attr(sim_vol_is, "score")
  sim_vol_oos <- simSortedBy( sim = sim,
                                 BT0 = BT0,
                                 Score = Score_vol,
                                 w_array = w_array,
                                 insample = FALSE )
  sim_vol_oos_score <- attr(sim_vol_oos, "score")
  lStats_vol_is <- descStats( sim_vol_is )
  stats_vol_is <- lStats_vol_is$stats
  lStats_vol_oos <- descStats( sim_vol_oos )
  stats_vol_oos <- lStats_vol_oos$stats
  
  
  # Quality
  Score_qual <- BT0$data$quality_mat
  sim_qual_is <- simSortedBy( sim = sim,
                             BT0 = BT0,
                             Score = Score_qual,
                             w_array = w_array,
                             insample = TRUE )
  sim_qual_is_score <- attr(sim_qual_is, "score")
  sim_qual_oos <- simSortedBy( sim = sim,
                              BT0 = BT0,
                              Score = Score_qual,
                              w_array = w_array,
                              insample = FALSE )
  sim_qual_oos_score <- attr(sim_qual_oos, "score")
  lStats_qual_is <- descStats( sim_qual_is )
  stats_qual_is <- lStats_qual_is$stats
  lStats_qual_oos <- descStats( sim_qual_oos )
  stats_qual_oos <- lStats_qual_oos$stats
  
  
  # Value 
  Score_val <- BT0$data$value_mat 
  sim_val_is <- simSortedBy( sim = sim,
                              BT0 = BT0,
                              Score = Score_val,
                              w_array = w_array,
                              insample = TRUE )
  sim_val_is_score <- attr(sim_val_is, "score")
  sim_val_oos <- simSortedBy( sim = sim,
                               BT0 = BT0,
                               Score = Score_val,
                               w_array = w_array,
                               insample = FALSE )
  sim_val_oos_score <- attr(sim_val_oos, "score")
  lStats_val_is <- descStats( sim_val_is )
  stats_val_is <- lStats_val_is$stats
  lStats_val_oos <- descStats( sim_val_oos )
  stats_val_oos <- lStats_val_oos$stats
  
  
  
  plot(density( na.omit(as.numeric(Score_val)))) # !!
  plot(density( na.omit(as.numeric(Score_qual)))) # !!
  
  
  
  # Plot
  
  colors <- fBasics::divPalette(n = ncol(sim), "RdYlGn")
  
  plot( x = 1:ncol(sim), y = stats_z_is["means", ], pch = 19, col = colors )
  plot( x = 1:ncol(sim), y = stats_z_oos["means", ], pch = 19, col = colors )
  plot( x = stats_z_is["sds", ], y = stats_z_is["means", ], pch = 19, col = colors )
  plot( x = stats_z_oos["sds", ], y = stats_z_oos["means", ], pch = 19, col = colors )
  
  plot( x = 1:ncol(sim), y = stats_mom_is["means", ], pch = 19, col = colors )
  plot( x = 1:ncol(sim), y = stats_mom_oos["means", ], pch = 19, col = colors )
  plot( x = stats_mom_is["sds", ], y = stats_mom_is["means", ], pch = 19, col = colors )
  plot( x = stats_mom_oos["sds", ], y = stats_mom_oos["means", ], pch = 19, col = colors )
  
  plot( x = 1:ncol(sim), y = stats_size_is["means", ], pch = 19, col = colors )
  plot( x = 1:ncol(sim), y = stats_size_oos["means", ], pch = 19, col = colors )
  plot( x = stats_size_is["sds", ], y = stats_size_is["means", ], pch = 19, col = colors )
  plot( x = stats_size_oos["sds", ], y = stats_size_oos["means", ], pch = 19, col = colors )
  
  plot( x = 1:ncol(sim), y = stats_vol_is["means", ], pch = 19, col = colors )
  plot( x = 1:ncol(sim), y = stats_vol_oos["means", ], pch = 19, col = colors )
  plot( x = stats_vol_is["sds", ], y = stats_vol_is["means", ], pch = 19, col = colors )
  plot( x = stats_vol_oos["sds", ], y = stats_vol_oos["means", ], pch = 19, col = colors )
  
  plot( x = 1:ncol(sim), y = stats_qual_is["means", ], pch = 19, col = colors )
  plot( x = 1:ncol(sim), y = stats_qual_oos["means", ], pch = 19, col = colors )
  plot( x = stats_qual_is["sds", ], y = stats_qual_is["means", ], pch = 19, col = colors )
  plot( x = stats_qual_oos["sds", ], y = stats_qual_oos["means", ], pch = 19, col = colors )
  
  plot( x = 1:ncol(sim), y = stats_val_is["means", ], pch = 19, col = colors )
  plot( x = 1:ncol(sim), y = stats_val_oos["means", ], pch = 19, col = colors )
  plot( x = stats_val_is["sds", ], y = stats_val_is["means", ], pch = 19, col = colors )
  plot( x = stats_val_oos["sds", ], y = stats_val_oos["means", ], pch = 19, col = colors )
  
  
  lSim_is <- list(  momentum = sim_mom_is,
                    value = sim_val_is,
                    quality = sim_qual_is,
                    size = sim_size_is,
                    lowvola = sim_vol_is,
                    rnd_sort = sim )
  lSim_oos <- list(  momentum = sim_mom_oos,
                     value = sim_val_oos,
                     quality = sim_qual_oos,
                     size = sim_size_oos,
                     lowvola = sim_vol_oos,
                     rnd_sort = sim )
  lStats_is <- list( momentum = lStats_mom_is,
                     value = lStats_val_is,
                     quality = lStats_qual_is,
                     size = lStats_size_is,
                     lowvola = lStats_vol_is,
                     rnd_sort = lStats_rnd )
  lStats_oos <- list( momentum = lStats_mom_oos,
                      value = lStats_val_oos,
                      quality = lStats_qual_oos,
                      size = lStats_size_oos,
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
  
  
  
  
  # Bear-Bull analysis
  
  BBSObj <- BBSRC$new()
  BBSObj$setCtrl()
  BBSObj$data$X_level <- cumulated(X_bm[ ,"Capw"], "discrete")
  
  lPhase_stats <- list()
  for ( factor_name in factor_names ) {
    BBSObj$data$X <- na.omit( cbind( X_bm[ ,c("QCapw", "Momentum")], lSim_oos[[ factor_name ]] ) )
    BBSObj$runRobust()
    lPhase_stats[[ factor_name ]] <- BBSObj$phaseStats()
  }
  

  
  par( mfrow = c(6, 2) )
  colors <- fBasics::divPalette(n = ncol(sim), "RdYlGn")
  for ( factor_name in factor_names ) {
    plot( x = lPhase_stats[[ factor_name ]]$`states==1`["sds", ], 
          y = lPhase_stats[[ factor_name ]]$`states==1`["means", ], pch = 19, col = colors )
    plot( x = lPhase_stats[[ factor_name ]]$`states==-1`["sds", ], 
          y = lPhase_stats[[ factor_name ]]$`states==-1`["means", ], pch = 19, col = colors )
  }
  
  
  par( mfrow = c(2, 1) )
  colors <- fBasics::divPalette(n = ncol(sim), "RdYlGn")
  factor_name <- "momentum"
  plot( x = lPhase_stats[[ factor_name ]]$`states==1`["sds", ], 
        y = lPhase_stats[[ factor_name ]]$`states==1`["means", ], pch = 19, col = colors )
  points( x = lPhase_stats[[ factor_name ]]$`states==1`["sds", "QCapw"], 
          y = lPhase_stats[[ factor_name ]]$`states==1`["means", "QCapw"], pch = 19, cex = 2, col = 1 )
  points( x = lPhase_stats[[ factor_name ]]$`states==1`["sds", "Momentum"], 
          y = lPhase_stats[[ factor_name ]]$`states==1`["means", "Momentum"], pch = 19, cex = 2, col = 4 )
  plot( x = lPhase_stats[[ factor_name ]]$`states==-1`["sds", ], 
          y = lPhase_stats[[ factor_name ]]$`states==-1`["means", ], pch = 19, col = colors )
  points( x = lPhase_stats[[ factor_name ]]$`states==-1`["sds", "QCapw"], 
          y = lPhase_stats[[ factor_name ]]$`states==-1`["means", "QCapw"], pch = 19, cex = 2, col = 1 )
  points( x = lPhase_stats[[ factor_name ]]$`states==-1`["sds", "Momentum"], 
          y = lPhase_stats[[ factor_name ]]$`states==-1`["means", "Momentum"], pch = 19, cex = 2, col = 4 )
  
  

  
  stats_field <- "sharpe"
  ldens_bull <- lapply( lPhase_stats, FUN = function(x) { density(x$`states==1`[stats_field, ]) })
  ldens_bear <- lapply( lPhase_stats, FUN = function(x) { density(x$`states==-1`[stats_field, ]) })
  par( mfrow = c(2, 1) )
  slolz:::plot.ldensity( ldens_bull, fillin = TRUE, main = "Bull" )
  slolz:::plot.ldensity( ldens_bear, fillin = TRUE, main = "Bear" )
  
  
  
  
  # # Use pre-computed factor matrices from Shiny tool
  #
  # #// The factor scores seem sometimes upside down ...? 
  # #// lowvola is reversed. 
  # #// Not so value, why?
  # 
  # 
  # env_fact <- readRDS( file = paste0("R:/Asset_Management/R/Shiny/Factor_Portfolios/Data/", universe, ".rds"))
  # 
  # factor_names <- c("momentum", "value", "quality", "lowvola", "lowsize", "highsize")
  # lFactor <- lapply( paste0(factor_names, "_mat"), FUN = get, envir = env_fact$data_cons )
  # names(lFactor) <- factor_names
  # 
  # # factor_names <- c("momentum", "value", "quality", "lowvola", "lowsize")
  # # lFactor <- lapply( paste0(factor_names, "_stdz_mat"), FUN = get, envir = env_fact$data_cons )
  # # names(lFactor) <- factor_names
  # 
  # lSim_is <- list()
  # lSim_oos <- list()
  # lScore_is <- list()
  # for ( factor_name in factor_names ) {
  #   Score_mat <- lFactor[[ factor_name ]]
  #   sim_is <- simSortedBy( sim = sim,
  #                          BT0 = BT0,
  #                          Score = Score_mat,
  #                          w_array = w_array,
  #                          insample = TRUE )
  #   sim_is_score <- attr(sim_is, "score")
  #   sim_oos <- simSortedBy( sim = sim,
  #                           BT0 = BT0,
  #                           Score = Score_mat,
  #                           w_array = w_array,
  #                           insample = FALSE )
  #   
  #   lSim_is[[ factor_name ]] <- sim_is
  #   lSim_oos[[ factor_name ]] <- sim_oos
  #   lScore_is[[ factor_name ]] <- sim_is_score
  # }
  # 
  # lStats_is <- lapply( lSim_is, FUN = descStats )
  # lStats_oos <- lapply( lSim_oos, FUN = descStats )
  # 
  # 
  # factor_name <- "momentum"
  # factor_name <- "value"
  # factor_name <- "quality"
  # factor_name <- "lowvola"
  # factor_name <- "lowsize"
  # factor_name <- "highsize"
  # 
  # plot( x = 1:ncol(sim), y = lStats_is[[ factor_name ]]$stats["means", ], pch = 19, col = colors )
  # plot( x = 1:ncol(sim), y = lStats_oos[[ factor_name ]]$stats["means", ], pch = 19, col = colors )
  # plot( x = lStats_is[[ factor_name ]]$stats["sds", ], y = lStats_is[[ factor_name ]]$stats["means", ], pch = 19, col = colors )
  # plot( x = lStats_oos[[ factor_name ]]$stats["sds", ], y = lStats_oos[[ factor_name ]]$stats["means", ], pch = 19, col = colors )
  # 
  # 
  # 
  # par( mfrow = c(3, 2) )
  # for ( factor_name in factor_names ) {
  #   plot( x = 1:ncol(sim),
  #         y = lStats_is[[ factor_name ]]$stats["means", ], 
  #         pch = 19, col = colors, main = factor_name )
  # }
  # 
  # par( mfrow = c(3, 2) )
  # for ( factor_name in factor_names ) {
  #   plot( x = lStats_is[[ factor_name ]]$stats["sds", ], 
  #         y = lStats_is[[ factor_name ]]$stats["means", ], 
  #         pch = 19, col = colors, main = factor_name )
  # }
  # 
  # par( mfrow = c(3, 2) )
  # for ( factor_name in factor_names ) {
  #   plot( x = lStats_oos[[ factor_name ]]$stats["sds", ], 
  #         y = lStats_oos[[ factor_name ]]$stats["means", ], 
  #         pch = 19, col = colors, main = factor_name )
  # }
  # 
  # lapply( factor_names, FUN = function(factor_name) {
  #   cor( x = lStats_oos[[ factor_name ]]$stats["sds", ], 
  #        y = lStats_oos[[ factor_name ]]$stats["means", ] )
  # } )
  
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # High Momentum predicts low variance
  # --------------------------------------------------------------------------
  
  
  
  ###
  # Risk return plot
  
  # Combine
  end_date <- "2022-03-06"
  strategy_names <- colnames(X_bm)
  sim_mom_all <- na.omit( cbind( X_bm[ ,strategy_names], sim_mom_oos) )
  sim_mom_all <- window( sim_mom_all, start(sim_all), end_date )
  lStats_mom_all <- descStats( sim_mom_all, descStatsSpec(annualize = TRUE) )
  
  
  lStats <- lStats_mom_all
  df <- data.frame( x = lStats$stats["sds", ],
                    y = lStats$stats["means", ] )
  # df <- data.frame( x = lStats_decade1$stats["sds", ],
  #                   y = lStats_decade1$stats["means", ] )
  # df <- data.frame( x = lStats_decade2$stats["sds", ],
  #                   y = lStats_decade2$stats["means", ] )
  # df <- data.frame( x = phase_stats$'states==1'["sds", ],
  #                   y = phase_stats$'states==1'["means", ] )
  # df <- data.frame( x = phase_stats$'states==-1'["sds", ],
  #                   y = phase_stats$'states==-1'["means", ] )
  h1 <- hist(df$x, breaks=150, plot=F)
  h2 <- hist(df$y, breaks=150, plot=F)
  top <- max(h1$counts, h2$counts)
  k <- kde2d( df$x, df$y, n = 100 )
  # margins
  oldpar <- par()
  par(mar=c(3,3,1,1))
  layout(matrix(c(2,0,1,3),2,2,byrow=TRUE),c(3,1), c(1,3))
  image(k)
  # abline( h = df$y[1], col = "white", lwd = 3 )
  # abline( v = df$x[1], col = "white", lwd = 3 )
  # abline( h = df$y[2], col = "white", lwd = 3 )
  # abline( v = df$x[2], col = "white", lwd = 3 )
  points( x = df$x, y = df$y, pch = 19, col = "orange", cex = 0.5 )
  tmp <- MASS::kde2d( x = df$x, y = df$y, n = 25, h = c(0.005, 0.02), lims = c(range(df$x), range(df$y)) ) 
  # contour( tmp, xlab = "Portfolio Standard Deviation", ylab = "Portfolio Return", 
  #          add = TRUE, drawlabels = FALSE, col = "red" )
  colors <- fBasics::rainbowPalette(n = ncol(X_bm))
  colors[1:2] <- c("black", "red")
  for ( i in 1:ncol(X_bm) ) {
    points( x = df$x[i], y = df$y[i], pch = 19, cex = 2, col = colors[i] )
  }
  legend( "topright", colnames(X_bm), lwd = 2, col = colors, text.col = colors, bty = "n")
  par(mar=c(0,2,1,0))
  barplot(h1$counts, axes=FALSE, ylim=c(0, top), space=0, col='grey')
  par(mar=c(2,0,0.5,1))
  barplot(h2$counts, axes=FALSE, xlim=c(0, top), space=0, col='grey', horiz=TRUE)
  
  ###
  
  
  colors <- c(rep(1, ncol(X_bm)), fBasics::divPalette(n = ncol(sim), "RdYlGn"))
  stats_mom_oos_all <- lStats_mom_all$stats
  plot( x = stats_mom_oos_all["sds", ], y = stats_mom_oos_all["means", ], pch = 19, col = colors )
  points( x = stats_mom_oos_all["sds", "Momentum"], y = stats_mom_oos_all["means", "Momentum"], pch = 19, cex = 3 )
  

  
  # Bear-Bull analysis
  
  BBSObj <- BBSRC$new()
  BBSObj$setCtrl()
  BBSObj$data$X_level <- cumulated(X_bm[ ,"Capw"], "discrete")
  BBSObj$data$X <- na.omit( cbind( X_bm[ ,c("QCapw", "Momentum")], sim_mom_oos ) )
  BBSObj$runRobust()
  phase_stats_mom <- BBSObj$phaseStats()
  
  boxplot( cbind( phase_stats_mom$`states==1`["means", ],
                  phase_stats_mom$`states==-1`["means", ] ), 
           beside = TRUE )
  abline( h = 0 )
  
  
  par( mfrow = c(1, 2) )
  colors <- fBasics::divPalette(n = ncol(sim), "RdYlGn")
  plot( x = phase_stats_mom$`states==1`["sds", ], y = phase_stats_mom$`states==1`["means", ], pch = 19, col = colors )
  points( x = phase_stats_mom$`states==1`["sds", "Momentum"], y = phase_stats_mom$`states==1`["means", "Momentum"], pch = 19, cex = 3 )
  points( x = phase_stats_mom$`states==1`["sds", "QCapw"], y = phase_stats_mom$`states==1`["means", "QCapw"], pch = 19, cex = 2 )
  plot( x = phase_stats_mom$`states==-1`["sds", ], y = phase_stats_mom$`states==-1`["means", ], pch = 19, col = colors )
  points( x = phase_stats_mom$`states==-1`["sds", "Momentum"], y = phase_stats_mom$`states==-1`["means", "Momentum"], pch = 19, cex = 3 )
  points( x = phase_stats_mom$`states==-1`["sds", "QCapw"], y = phase_stats_mom$`states==-1`["means", "QCapw"], pch = 19, cex = 2 )
  
  
  
    
  
  # # --------------------------------------------------------------------------
  # # SAVE
  # # --------------------------------------------------------------------------
  # 
  # env <- new.env()
  # env$lSim <- lSim
  # env$lWeights <- lWeights
  # env$llStats <- llStats
  # env$lPortf_scores <- lPortf_scores
  # # saveRDS( env, file = paste0(wd, "waRehouse/alphabet_strategy_rps.rds") )
  # saveRDS( env, file = paste0(wd, "waRehouse/alphabet_strategy_rps_bmcentered.rds") )
  # 
  # 
  # 
  # 
  # 
  # # --------------------------------------------------------------------------
  # # LOAD
  # # --------------------------------------------------------------------------
  # 
  # # env <- readRDS( file = paste0(wd, "waRehouse/alphabet_strategy_rps.rds") )
  # env <- readRDS( file = paste0(wd, "waRehouse/alphabet_strategy_rps_bmcentered.rds") )
  # lSim <- env$lSim
  # lWeights <- env$lWeights
  # llStats <- env$llStats
  # lPortf_scores <- env$lPortf_scores
  
  
  
  
  
  
 
  