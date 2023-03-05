
  
  ############################################################################
  ### RANDOM PORTFOLIOS - MAXRC SAMPLING - FACTOR CONTROL - USA - ANALYZE
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     05.01.2023
  # First version:    05.01.2023
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
  source( paste0(wd, "Source/class_BacktestPanelFactor.R") )
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
  end_date <- "2022-03-06"
  
  
  
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
  # Initialize Backtest object
  # --------------------------------------------------------------------------
  
  BT0 <- loadBacktest( wd = paste0(wd_data, "Data/"),
                       name = "usa_maxw5" )
  rebdates <- BT0$spec$rebdates
  BT0$spec$rebdates <- rebdates[ seq(from = 1, to = length(rebdates), by = 2) ]
  BT0$initRBE()
  
  
  
  
 
  # --------------------------------------------------------------------------
  # LOAD BACKTESTS
  # --------------------------------------------------------------------------
  
  BT_Fctrl <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "multifac_fctrl" )
  
  BT_qcapw <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "qcapw" )
  sim_qcapw <- BT_qcapw$output$simulations
  
  BT_qeqw <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "qeqw" )
  sim_qeqw <- BT_qeqw$output$simulations
  
  BT1 <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "maxRC_multifac_usa" )
  sim_1 <- BT1$output$simulations

  BT2 <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "maxRC_multifac_usa_valctrl" )
  sim_2 <- BT2$output$simulations
  
  BT3 <- loadBacktest( wd =paste0(wd_data, "waRehouse/"), name = "maxRC_multifac_usa_qualctrl" )
  sim_3 <- BT3$output$simulations
  
  BT4 <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "maxRC_multifac_usa_sizectrl" )
  sim_4 <- BT4$output$simulations

  BT5 <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "maxRC_multifac_usa_sizestdzctrl" )
  sim_5 <- BT5$output$simulations
 
  BT6 <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "maxRC_multifac_usa_valqualsizectrl" )
  sim_6 <- BT6$output$simulations

  BT7 <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "maxRC_multifac_usa_momctrl" )
  sim_7 <- BT7$output$simulations

  BT8 <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "maxRC_multifac_usa_valqualsizemomctrl" )
  sim_8 <- BT8$output$simulations
  
  
  # Add qeqw and qcapw to X_bm
  X_bm <- na.omit( cbind( X_fact, QEQW = sim_qeqw, QCapw = sim_qcapw ) )
  
  descStats( X_bm )
  
  
  
  
  # --------------------------------------------------------------------------
  # SORT RP'S BY FACTOR SCORE, CHECK RANK VS. PERFORMANCE
  # --------------------------------------------------------------------------
  

  w_array_1 <- list2array( L = BT1$output[ grepl("maxRC", names(BT1$output)) ] )
  w_array_2 <- list2array( L = BT2$output[ grepl("maxRC", names(BT2$output)) ] )
  w_array_3 <- list2array( L = BT3$output[ grepl("maxRC", names(BT3$output)) ] )
  w_array_4 <- list2array( L = BT4$output[ grepl("maxRC", names(BT4$output)) ] )
  w_array_5 <- list2array( L = BT5$output[ grepl("maxRC", names(BT5$output)) ] )
  w_array_6 <- list2array( L = BT6$output[ grepl("maxRC", names(BT6$output)) ] )
  w_array_7 <- list2array( L = BT7$output[ grepl("maxRC", names(BT7$output)) ] )
  w_array_8 <- list2array( L = BT8$output[ grepl("maxRC", names(BT8$output)) ] )
  
 
  lW_array <- list( mfac = w_array_1,
                    mfac_val_ctrl = w_array_2,
                    mfac_qual_ctrl = w_array_3,
                    mfac_size_ctrl = w_array_4,
                    mfac_size_stdz_ctrl = w_array_5,
                    mfac_val_qual_size_ctrl = w_array_6,
                    mfac_mom_ctrl = w_array_7,
                    mfac_val_qual_size_mom_ctrl = w_array_8 )
  
  lSim_backtests <- list( mfac = BT1$output$simulations,
                          mfac_val_ctrl = BT2$output$simulations,
                          mfac_qual_ctrl = BT3$output$simulations,
                          mfac_size_ctrl = BT4$output$simulations,
                          mfac_size_stdz_ctrl = BT5$output$simulations,
                          mfac_val_qual_size_ctrl = BT6$output$simulations,
                          mfac_mom_ctrl = BT7$output$simulations,
                          mfac_val_qual_size_mom_ctrl = BT8$output$simulations )
  
  
  lSim_is <- list()
  lSim_oos <- list()
  bt_names <- names(lW_array)
  factor_names <- names(BT_Fctrl$data$lFactor)
  # bt_names <- c("mfac", "mfac_size_stdz_ctrl")
  # factor_names <- c("momentum", "momentum_stdz", "value", "quality",  "size_stdz", "lowvola_stdz")
  for ( bt_name in bt_names ) {
    
    sim <- lSim_backtests[[ bt_name ]]
    lSim_is_tmp <- lSim_oos_tmp <- list()
    lSim_is_tmp[["rnd"]] <- sim
    lSim_oos_tmp[["rnd"]] <- sim
    
    for ( factor_name in factor_names ) {
      
      w_array <- lW_array[[ bt_name ]]
      score_mat <- BT_Fctrl$data$lFactor[[ factor_name ]]
      
      sim_is <- simSortedBy( sim = sim,
                             # BT0 = BT_Fctrl,
                             BT0 = BT0,
                             Score = score_mat,
                             w_array = w_array,
                             insample = TRUE )
      sim_oos <- simSortedBy( sim = sim,
                              # BT0 = BT_Fctrl,
                              BT0 = BT0,
                              Score = score_mat,
                              w_array = w_array,
                              insample = FALSE )
      lSim_is_tmp[[ factor_name ]] <- sim_is
      lSim_oos_tmp[[ factor_name ]] <- sim_oos
    }
    lSim_is[[ bt_name ]] <- lSim_is_tmp
    lSim_oos[[ bt_name ]] <- lSim_oos_tmp
  }

  lStats_is <- lapply( lSim_is, FUN = function(X) { lapply( X, FUN = descStats ) } )
  lStats_oos <- lapply( lSim_oos, FUN = function(X) { lapply( X, FUN = descStats ) } )
  
  
  
  
  
  lapply( lStats_is, FUN = function(X) { lapply( X, FUN = function(x) { range(x$stats["sds", ]) } ) } )
  lapply( lStats_is, FUN = function(X) { lapply( X, FUN = function(x) { range(x$stats["means", ]) } ) } )
  
  
  
  # --------------------------------------------------------------------------
  # CHECK FACTOR EXPOSURE - EX POST
  # --------------------------------------------------------------------------
  
  lFactorScores <- list()
  # factor_names <- c("momentum_stdz", "value", "quality", "size_stdz", "lowvola_stdz")
  factor_names <- names(BT_Fctrl$data$lFactor)
  for ( i in 1:length(lW_array) ) {
    lFactorScores_tmp <- list()
    for ( fac_name in factor_names ) {
      lWmat <- lapply( 1:dim(lW_array[[i]])[[3]], FUN = function(j) { lW_array[[i]][ ,,j] } )
      lFacScore <- portfolioScore( wmat = lWmat, 
                                   score_mat = BT_Fctrl$data$lFactor[[fac_name]] )
      fac_scores <- do.call( cbind, lFacScore )
      lFactorScores_tmp[[fac_name]] <- fac_scores
    }
    lFactorScores[[i]] <- lFactorScores_tmp
  }
  names(lFactorScores) <- names(lW_array)
  
  
  
  factor_names <- names(lStats_is[[1]])[-length(lStats_is[[1]])]
  lStats_is <- lapply( lStats_is, function(x) { x[ factor_names ] } )
  lStats_oos <- lapply( lStats_oos, function(x) { x[ factor_names ] } )
  
  
  
  
  # --------------------------------------------------------------------------
  # SAVE
  # --------------------------------------------------------------------------
  
  # env <- new.env()
  # env$lW_array <- lW_array
  # env$lSim_backtests <- lSim_backtests
  # env$lSim_is <- lSim_is
  # env$lSim_oos <- lSim_oos
  # env$lStats_is <- lStats_is
  # env$lStats_oos <- lStats_oos
  # env$lFactorScores <- lFactorScores
  # 
  # saveRDS( object = env, file = paste0(wd_data, "waRehouse/rp_maxRC_multifac_usa_factor_control_analze.rds") )
  
  
  
  # --------------------------------------------------------------------------
  # LOAD
  # --------------------------------------------------------------------------
  
  env <- readRDS( file = paste0(wd_data, "waRehouse/rp_maxRC_multifac_usa_factor_control_analze.rds") )
  lW_array <- env$lW_array
  lSim_backtests <- env$lSim_backtests
  lSim_is <- env$lSim_is
  lSim_oos <- env$lSim_oos
  lStats_is <- env$lStats_is
  lStats_oos <- env$lStats_oos
  lFactor_scores <- env$lFactorScores
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # ANALYZE
  # --------------------------------------------------------------------------
  
  
  
  
  bt_name <- "mfac"
  bt_name <- "mfac_val_ctrl"
  factor_names <- names(lStats_oos[[ bt_name ]])
  sim <- lSim_backtests[[ bt_name ]]
  colors <- fBasics::divPalette(n = ncol(sim), "RdYlGn")
  mfrow <- c(4, 2)
  
  par( mfrow = mfrow )
  for ( factor_name in factor_names ) {
    plot( x = 1:ncol(sim),
          y = lStats_is[[ bt_name ]][[ factor_name ]]$stats["means", ], 
          pch = 19, col = colors, main = factor_name )
  }
  par( mfrow = mfrow )
  for ( factor_name in factor_names ) {
    plot( x = 1:ncol(sim),
          y = lStats_oos[[ bt_name ]][[ factor_name ]]$stats["means", ], 
          pch = 19, col = colors, main = factor_name )
  }
  par( mfrow = mfrow )
  for ( factor_name in  factor_names ) {
    plot( x = lStats_is[[ bt_name ]][[ factor_name ]]$stats["sds", ], 
          y = lStats_is[[ bt_name ]][[ factor_name ]]$stats["means", ], 
          pch = 19, col = colors, main = factor_name )
  }
  par( mfrow = mfrow )
  for ( factor_name in factor_names ) {
    plot( x = lStats_oos[[ bt_name ]][[ factor_name ]]$stats["sds", ], 
          y = lStats_oos[[ bt_name ]][[ factor_name ]]$stats["means", ], 
          pch = 19, col = colors, main = factor_name )
  }
  
  
  
  par( mfrow = mfrow )
  for ( factor_name in factor_names ) {
    plot( x = 1:ncol(sim),
          y = lStats_is[[ bt_name ]][[ factor_name ]]$stats["sds", ], 
          pch = 19, col = colors, main = factor_name )
  }
  par( mfrow = mfrow )
  for ( factor_name in factor_names ) {
    plot( x = 1:ncol(sim),
          y = lStats_oos[[ bt_name ]][[ factor_name ]]$stats["sds", ], 
          pch = 19, col = colors, main = factor_name )
  }
  
  
  
  
  
  
  fac_name <- "momentum_stdz"
  fac_name <- "lowvola"
  fac_name <- "size_stdz"
  fac_name <- "value"
  fac_name <- "quality"
  plot( x = lStats_oos[["mfac"]][[fac_name]]$stats["sds", ], 
        y = lStats_oos[["mfac"]][[fac_name]]$stats["means", ], pch = 19,
        xlim = range( unlist( lapply( lStats_oos, FUN = function(x) { x[[fac_name]]$stats["sds", ]}))),
        ylim = range( unlist( lapply( lStats_oos, FUN = function(x) { x[[fac_name]]$stats["means", ]}))) )
  # points( x = lStats_oos[["mfac_val_ctrl"]][[fac_name]]$stats["sds", ], 
  #         y = lStats_oos[["mfac_val_ctrl"]][[fac_name]]$stats["means", ], pch = 19, col = 2 )
  # points( x = lStats_oos[["mfac_qual_ctrl"]][[fac_name]]$stats["sds", ], 
  #         y = lStats_oos[["mfac_qual_ctrl"]][[fac_name]]$stats["means", ], pch = 19, col = 3 )
  # points( x = lStats_oos[["mfac_size_ctrl"]][[fac_name]]$stats["sds", ], 
  #         y = lStats_oos[["mfac_size_ctrl"]][[fac_name]]$stats["means", ], pch = 19, col = 4 )
  points( x = lStats_oos[["mfac_size_stdz_ctrl"]][[fac_name]]$stats["sds", ], 
          y = lStats_oos[["mfac_size_stdz_ctrl"]][[fac_name]]$stats["means", ], pch = 19, col = 5 )
  points( x = lStats_oos[["mfac_val_qual_size_ctrl"]][[fac_name]]$stats["sds", ], 
          y = lStats_oos[["mfac_val_qual_size_ctrl"]][[fac_name]]$stats["means", ], pch = 19, col = 6 )
  # points( x = lStats_oos[["mfac_mom_ctrl"]][[fac_name]]$stats["sds", ], 
  #         y = lStats_oos[["mfac_mom_ctrl"]][[fac_name]]$stats["means", ], pch = 19, col = 7 )
  # points( x = lStats_oos[["mfac_val_qual_size_mom_ctrl"]][[fac_name]]$stats["sds", ], 
  #         y = lStats_oos[["mfac_val_qual_size_mom_ctrl"]][[fac_name]]$stats["means", ], pch = 19, col = 8 )
  
  
  
  mu <- do.call( cbind, lapply( lStats_oos, FUN = function(x) { x[[fac_name]]$stats["means", ] } ) )
  sds <- do.call( cbind, lapply( lStats_oos, FUN = function(x) { x[[fac_name]]$stats["sds", ] } ) )
  sharpe <- do.call( cbind, lapply( lStats_oos, FUN = function(x) { x[[fac_name]]$stats["sharpe", ] } ) )
  
  boxplot( mu, beside = TRUE )
  boxplot( sds, beside = TRUE )
  boxplot( sharpe, beside = TRUE )
  
  
  
  
  
  fac_name <- "lowvola"
  sim <- lSim_backtests[[1]]
  colors <- fBasics::divPalette(n = ncol(sim), "RdYlGn")
  xlab <- ylab <- ""
  xlim_is <- range( unlist( lapply( lStats_is, FUN = function(x) { x[[fac_name]]$stats["sds", ] } ) ) )
  xlim_oos <- range( unlist( lapply( lStats_oos, FUN = function(x) { x[[fac_name]]$stats["sds", ] } ) ) )
  statistic <- "sds"
  ylim_is <- range( unlist( lapply( lStats_is, FUN = function(x) { x[[fac_name]]$stats[statistic, ] } ) ) )
  ylim_oos <- range( unlist( lapply( lStats_oos, FUN = function(x) { x[[fac_name]]$stats[statistic, ] } ) ) )
  mfrow <- c(4, 2)
  
  par( mfrow = mfrow )
  for ( bt_name in names(lStats_is) ) {
    plot( x = 1:ncol(sim),
          y = lStats_is[[ bt_name ]][[ fac_name ]]$stats[statistic, ], 
          pch = 19, col = colors, main = bt_name,
          xlab = xlab, ylab = ylab, ylim = ylim_is )
  }
  par( mfrow = mfrow )
  for ( bt_name in names(lStats_oos) ) {
    plot( x = 1:ncol(sim),
          y = lStats_oos[[ bt_name ]][[ fac_name ]]$stats[statistic, ], 
          pch = 19, col = colors, main = bt_name,
          xlab = xlab, ylab = ylab, ylim = ylim_oos )
  }
  par( mfrow = mfrow )
  for ( bt_name in names(lStats_is) ) {
    plot( x = lStats_is[[ bt_name ]][[ fac_name ]]$stats["sds", ], 
          y = lStats_is[[ bt_name ]][[ fac_name ]]$stats["means", ], 
          pch = 19, col = colors, main = bt_name,
          xlab = xlab, ylab = ylab, xlim = xlim_is, ylim = ylim_is )
  }
  par( mfrow = mfrow )
  for ( bt_name in names(lStats_oos) ) {
    plot( x = lStats_oos[[ bt_name ]][[ fac_name ]]$stats["sds", ], 
          y = lStats_oos[[ bt_name ]][[ fac_name ]]$stats["means", ], 
          pch = 19, col = colors, main = bt_name,
          xlab = xlab, ylab = ylab, xlim = xlim_oos, ylim = ylim_oos )
  }
  
  
  
  
  
  # --------------------------------------------------------------------------
  # CORRELATIONS 
  # --------------------------------------------------------------------------
  
  # Correlations between means and sds (low-vola anomaly)
  
  bt_names <- names(lStats_oos)
  factor_names <- names(lStats_oos[[1]])
  
  lCorr_is <- list()
  lCorr_oos <- list()
  for ( bt_name in bt_names ) {
    tmp <- lapply( factor_names, FUN = function(factor_name) {
      cor( x = lStats_is[[ bt_name ]][[ factor_name ]]$stats["sds", ], 
           y = lStats_is[[ bt_name ]][[ factor_name ]]$stats["means", ] )
    } )
    lCorr_is[[ bt_name ]] <- setNames( unlist( tmp ), factor_names )
    tmp <- lapply( factor_names, FUN = function(factor_name) {
      cor( x = lStats_oos[[ bt_name ]][[ factor_name ]]$stats["sds", ], 
           y = lStats_oos[[ bt_name ]][[ factor_name ]]$stats["means", ] )
    } )
    lCorr_oos[[ bt_name ]] <- setNames( unlist( tmp ), factor_names )
  }
  correlations_is <- t( do.call( cbind, lCorr_is ) )
  correlations_oos <- t( do.call( cbind, lCorr_oos ) )
  
  correlations_is
  correlations_oos

  
  
  # Correlations between some statistic and rank
  
  statistic <- "means"
  statistic <- "sds"
  statistic <- "sharpe"
  
  bt_names <- names(lStats_oos)
  factor_names <- names(lStats_oos[[1]])
  
  lCorr_is <- list()
  lCorr_oos <- list()
  for ( bt_name in bt_names ) {
    tmp <- lapply( factor_names, FUN = function(factor_name) {
      cor( x = lStats_is[[ bt_name ]][[ factor_name ]]$stats[statistic, ], 
           y = 1:ncol(lStats_is[[1]][[1]]$stats) )
    } )
    lCorr_is[[ bt_name ]] <- setNames( unlist( tmp ), factor_names )
    tmp <- lapply( factor_names, FUN = function(factor_name) {
      cor( x = lStats_oos[[ bt_name ]][[ factor_name ]]$stats[statistic, ], 
           y = 1:ncol(lStats_oos[[1]][[1]]$stats)  )
    } )
    lCorr_oos[[ bt_name ]] <- setNames( unlist( tmp ), factor_names )
  }
  correlations_is <- t( do.call( cbind, lCorr_is ) )
  correlations_oos <- t( do.call( cbind, lCorr_oos ) )
  
  correlations_is
  correlations_oos

  
  
  
  # --------------------------------------------------------------------------
  # CHECK FACTOR EXPOSURE - EX POST
  # --------------------------------------------------------------------------

  
  bt_name <- "mfac_size_stdz_ctrl"
  fac_name <- "size_stdz"
  dates <- intersect( rownames(BT_Fctrl$data$factor_scores_bm), 
                      rownames( lFactorScores[[bt_name]][[fac_name]]) )
  tmp <- cbind( BT_Fctrl$data$factor_scores_bm[ ,fac_name][dates, ], 
                BT_Fctrl$data$factor_scores_bm[ ,fac_name][dates, ] - 0.1,
                BT_Fctrl$data$factor_scores_bm[ ,fac_name][dates, ] + 0.1,
                lFactorScores[[bt_name]][[fac_name]][dates, ] )
  plot( tmp, plot.type = "single", col = "grey" )
  lines( tmp[ ,1], lwd = 2 )
  lines( tmp[ ,2], col = 2, lwd = 2 )
  lines( tmp[ ,3], col = 2, lwd = 2 )
  
  
  
  
  
  # bt_name <- "mfac_size_stdz_ctrl"
  # fac_name <- "size_stdz"
  # bt_name <- "mfac_val_ctrl"
  # fac_name <- "value"
  bt_name <- "mfac_qual_ctrl"
  fac_name <- "quality"
  lWmat <- lapply( 1:dim(lW_array[[bt_name]])[[3]], FUN = function(j) { lW_array[[bt_name]][ ,,j] } )
  # lWmat <- lapply( 1:10, FUN = function(j) { lW_array[[bt_name]][ ,,j] } )
  lFacScore <- portfolioScore( wmat = lWmat, 
                               score_mat = BT_Fctrl$data$lFactor[[fac_name]] )
  fac_scores <- do.call( cbind, lFacScore )
  tmp <- cbind( BT_Fctrl$data$factor_scores_bm[ ,fac_name], 
                BT_Fctrl$data$factor_scores_bm[ ,fac_name] - 0.1,
                BT_Fctrl$data$factor_scores_bm[ ,fac_name] + 0.1 )
  dates <- intersect( rownames(fac_scores), rownames(tmp) )
  plot( fac_scores[dates, ], plot.type = "single", col = "grey" )
  lines( tmp[dates ,1], lwd = 2 )
  lines( tmp[dates ,2], col = 2, lwd = 2 )
  lines( tmp[dates ,3], col = 2, lwd = 2 )
  
  
  plot( na.omit(BT_Fctrl$data$factor_scores_bm) )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # PLOT RP-PERFORMANCE
  # --------------------------------------------------------------------------
  
 
  # Risk return plot
  riskReturnPlot <- function( lStats, lStats_sim1 = NULL, X_bm )
  {
    # lStats <- lStats_all
    df <- data.frame( x = lStats$stats["sds", ],
                      y = lStats$stats["means", ] )
    rownames(df) <- colnames(lStats$stats)
    h1 <- hist(df$x, breaks=250, plot=F)
    h2 <- hist(df$y, breaks=250, plot=F)
    top <- max(h1$counts, h2$counts) * 2
    lims_x <- range(c(df$x * 0.995, df$x * 1.005))
    # lims_x <- range(df$x)
    lims_y <- range(df$y)
    k <- kde2d( df$x, df$y, n = 250, lims = c(lims_x, lims_y) )
    # margins
    oldpar <- par()
    par( mar = c(3, 3, 1, 1) )
    layout( mat = matrix(c(2,0,1,3), 2, 2, byrow = TRUE), c(3, 1), c(1, 3) )
    image(k)
    idx_x <- c( which(df$x > quantile(df$x, 0.9)), which(df$x < quantile(df$x, 0.1)) )
    idx_y <- c( which(df$y > quantile(df$y, 0.9)), which(df$y < quantile(df$y, 0.1)) )
    idx <- c( idx_x, idx_y )
    points( x = df$x[idx], y = df$y[idx], pch = 19, col = "orange", cex = 0.5 )
    # tmp <- MASS::kde2d( x = df$x, y = df$y, n = 25, h = c(0.005, 0.02), lims = c(range(df$x), range(df$y)) ) 
    # contour( tmp, xlab = "Portfolio Standard Deviation", ylab = "Portfolio Return",
    #          add = TRUE, drawlabels = FALSE, col = "red" )
    if ( !is.null(lStats_sim1) ) {
      df1 <- data.frame( x = lStats_sim1$stats["sds", ],
                         y = lStats_sim1$stats["means", ] )
      k1 <- kde2d( df1$x, df1$y, n = 250 )
      tmp <- MASS::kde2d( x = df1$x+0.001, y = df1$y, n = 25, lims = c(range(df1$x+0.001), range(df1$y)) ) 
      contour( tmp, xlab = "Portfolio Standard Deviation", ylab = "Portfolio Return",
               add = TRUE, drawlabels = FALSE, col = "blue" )
    }
    colors <- fBasics::rainbowPalette(n = ncol(X_bm))
    colors[1:2] <- c("black", "red")
    for ( i in 1:ncol(X_bm) ) {
      points( x = df$x[i], y = df$y[i], pch = 19, cex = 2, col = colors[i] )
    }
    legend( "topright", colnames(X_bm), lwd = 2, col = colors, text.col = colors, bty = "n")
    par(mar = c(0,2,1,0))
    barplot(h1$counts, axes=FALSE, ylim=c(0, top), space=0, col='grey')
    par(mar = c(2,0,0.5,1))
    barplot(h2$counts, axes=FALSE, xlim=c(0, top), space=0, col='grey', horiz=TRUE)
  }
  
  
  X_bm <- X_fact
  # sim_1 <- BT1$output$simulations
  sim_1 <- BT5$output$simulations
  dates <- intersect( rownames(X_bm), rownames(sim_1) )
  dates <- dates[ dates <= end_date ]
  X_bm <- X_bm[dates, ]
  sim_1 <- sim_1[dates, ]
  strategy_names <- colnames(X_bm)
  sim_all <- na.omit( cbind( X_bm[ ,strategy_names], lSim_backtests[["mfac_val_qual_size_ctrl"]]) )
  lStats_all <-  descStats( sim_all, descStatsSpec(annualize = TRUE) )
  lStats_sim1 <- descStats( sim_1, descStatsSpec(annualize = TRUE) )
  
 
  
  riskReturnPlot( lStats = lStats_all, 
                  lStats_sim1 = lStats_sim1, 
                  X_bm = X_bm[ ,strategy_names] )
  
  riskReturnPlot( lStats = lStats_all, 
                  X_bm = X_bm[ ,strategy_names] )
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BEARS AND BULLS
  # --------------------------------------------------------------------------
  
  BBSObj <- BBSRC$new()
  BBSObj$setCtrl()
  BBSObj$data$X_level <- cumulated(X_bm[ ,"Capw"], "discrete")
  
  bm_names <- c("Capw", "QCapw", "QEQW", "Momentum", "Value", "Quality", "Size", "Enhanced_Value", "Multifactor" )
  bt_names <- names(lSim_oos)
  factor_names <- names(lStats_oos[[1]])
  lPhase_stats <- list()
  for( bt_name in bt_names ) {
    lPhase_stats_tmp <- list()
    for ( factor_name in factor_names ) {
      BBSObj$data$X <- na.omit( cbind( X_bm[ ,bm_names], 
                                       lSim_oos[[ bt_name ]][[ factor_name ]] ) )
      BBSObj$runRobust()
      lPhase_stats_tmp[[ factor_name ]] <- BBSObj$phaseStats()
    }
    lPhase_stats[[ bt_name ]] <- lPhase_stats_tmp
  }
 
  
  
  
  par( mfrow = c(2, 1) )
  colors <- fBasics::divPalette(n = ncol(lPhase_stats$mfac$rnd$`states==1`), "RdYlGn")
  bt_name <- "mfac_size_stdz_ctrl"
  factor_name <- "lowvola"
  plot( x = lPhase_stats[[ bt_name ]][[ factor_name ]]$`states==1`["sds", ], 
        y = lPhase_stats[[ bt_name ]][[ factor_name ]]$`states==1`["means", ], pch = 19, col = colors )
  for ( i in seq(along = bm_names) ) {
    points( x = lPhase_stats[[ bt_name ]][[ factor_name ]]$`states==1`["sds", bm_names[i]], 
            y = lPhase_stats[[ bt_name ]][[ factor_name ]]$`states==1`["means", bm_names[i]], 
            pch = 19, cex = 2, col = i )
  }
  plot( x = lPhase_stats[[ bt_name ]][[ factor_name ]]$`states==-1`["sds", ], 
        y = lPhase_stats[[ bt_name ]][[ factor_name ]]$`states==-1`["means", ], pch = 19, col = colors )
  for ( i in seq(along = bm_names) ) {
    points( x = lPhase_stats[[ bt_name ]][[ factor_name ]]$`states==-1`["sds", bm_names[i]], 
            y = lPhase_stats[[ bt_name ]][[ factor_name ]]$`states==-1`["means", bm_names[i]], 
            pch = 19, cex = 2, col = i )
  }
  colors <- 1:length(bm_names)
  legend( "topleft", bm_names, lwd = 2, col = colors, text.col = colors, bty = "n")
  
  
  
  
  bt_name <- "mfac"
  stats_field <- "sharpe"
  ldens_bull <- lapply( lPhase_stats[[bt_name]], FUN = function(x) { density(x$`states==1`[stats_field, ]) })
  ldens_bear <- lapply( lPhase_stats[[bt_name]], FUN = function(x) { density(x$`states==-1`[stats_field, ]) })
  par( mfrow = c(2, 1) )
  slolz:::plot.ldensity( ldens_bull, fillin = TRUE, main = "Bull" )
  slolz:::plot.ldensity( ldens_bear, fillin = TRUE, main = "Bear" )
  
  
  
  
  
  
  
  
  
  