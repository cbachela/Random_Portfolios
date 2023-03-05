  
  
  ############################################################################
  ### MAXIMUM FACTOR EXPOSURE PORTFOLIOS - MSCI MULFIACTOR CONSTRAINTS - USA
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     09.02.2023
  # First version:    09.02.2023
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
  # PARAMETERS
  # --------------------------------------------------------------------------
  
  universe <- "usa"
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  end_date <- "2022-03-06"
  fc <- 0
  vc <- 0
  
  
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
  
  
  
  env <- readRDS( file = paste0(wd_data, "waRehouse/rp_maxRC_multifac_usa_factor_control_analze.rds") )
  lW_array <- env$lW_array
  lSim_backtests <- env$lSim_backtests
  lSim_is <- env$lSim_is
  lSim_oos <- env$lSim_oos
  lStats_is <- env$lStats_is
  lStats_oos <- env$lStats_oos
  lFactor_scores <- env$lFactorScores
  
  
  
  # --------------------------------------------------------------------------
  # MAXIMUM FACTOR EXPOSURE PORTFOLIOS
  # --------------------------------------------------------------------------
  
  # Prepare default specifications
  BT_maxExp <- BT_Fctrl$copy()
  BT_maxExp$spec$cons_fun <- c("box", "sector", "sectorLB", "vola") #, "controlbmex10")
  # BT_maxExp$spec$factor_names_constraints <- c("value", "quality", "size_stdz", "momentum_stdz")
  BT_maxExp$spec$portfolio <- "maxExposure"
  BT_maxExp$spec$direction <- -1
  BT_maxExp$spec$fc <- fc
  BT_maxExp$spec$vc <- vc
  
  # Maximum momentum portfolio
  factor_name <- "momentum_stdz"
  BT_maxmom <- BT_maxExp$copy()
  BT_maxmom$spec$factor_name <- factor_name
  BT_maxmom$runLoop()
  sim_maxmom <- BT_maxmom$simulate()
  
  BT_minmom <- BT_maxExp$copy()
  BT_minmom$spec$factor_name <- factor_name
  BT_minmom$spec$direction <- 1
  BT_minmom$runLoop()
  sim_minmom <- BT_minmom$simulate()
  
  # Maximum value portfolio
  factor_name <- "value"
  BT_maxval <- BT_maxExp$copy()
  BT_maxval$spec$factor_name <- factor_name
  BT_maxval$runLoop()
  sim_maxval <- BT_maxval$simulate()
  
  BT_minval <- BT_maxExp$copy()
  BT_minval$spec$factor_name <- factor_name
  BT_minval$spec$direction <- 1
  BT_minval$runLoop()
  sim_minval <- BT_minval$simulate()
  
  # Maximum quality portfolio
  factor_name <- "quality"
  BT_maxqual <- BT_maxExp$copy()
  BT_maxqual$spec$factor_name <- factor_name
  BT_maxqual$runLoop()
  sim_maxqual <- BT_maxqual$simulate()
  
  BT_minqual <- BT_maxExp$copy()
  BT_minqual$spec$factor_name <- factor_name
  BT_minqual$spec$direction <- 1
  BT_minqual$runLoop()
  sim_minqual <- BT_minqual$simulate()
  
  # Maximum size portfolio
  factor_name <- "size_stdz"
  BT_maxsize <- BT_maxExp$copy()
  BT_maxsize$spec$factor_name <- factor_name
  BT_maxsize$runLoop()
  sim_maxsize <- BT_maxsize$simulate()
  
  BT_minsize <- BT_maxExp$copy()
  BT_minsize$spec$factor_name <- factor_name
  BT_minsize$spec$direction <- 1
  BT_minsize$runLoop()
  sim_minsize <- BT_minsize$simulate()
  
  
  sim <- na.omit( cbind( QCapw = sim_qcapw,
                         QEQW = sim_qeqw,
                         maxMomentum = sim_maxmom, 
                         maxValue = sim_maxval, 
                         maxQuality = sim_maxqual, 
                         maxSize = sim_maxsize,
                         minMomentum = sim_minmom, 
                         minValue = sim_minval, 
                         minQuality = sim_minqual, 
                         minSize = sim_minsize ) )
  X_bm <- na.omit( cbind( X_fact, sim ) )
  statsolz:::plot.stats( descStats(X_bm), sortby = NULL )
  
  
  
  # --------------------------------------------------------------------------
  # PLOT RP-PERFORMANCE
  # --------------------------------------------------------------------------
  
  
  # Risk return plot
  riskReturnPlot <- function( lStats, lStats_sim1 = NULL, X_bm, colors = NULL )
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
    if ( is.null(colors) ) {
      colors <- fBasics::rainbowPalette( n = ncol(X_bm) )
      colors[1:4] <- c("black", "red", "darkgrey", "lightgrey")
    }
    for ( i in 1:ncol(X_bm) ) {
      points( x = df$x[i], y = df$y[i], pch = 19, cex = 2, col = colors[i] )
    }
    legend( "topleft", colnames(X_bm), lwd = 2, col = colors, text.col = colors, bty = "n", cex = 1 )
    par(mar = c(0,2,1,0))
    barplot(h1$counts, axes=FALSE, ylim=c(0, top), space=0, col='grey')
    par(mar = c(2,0,0.5,1))
    barplot(h2$counts, axes=FALSE, xlim=c(0, top), space=0, col='grey', horiz=TRUE)
  }
  
  bt_name <- "mfac"
  bm_names <- c("Capw", "Multifactor", colnames(sim) )
  factor_names <- c("momentum_stdz", "value", "quality", "size_stdz", "lowvola_stdz")
  sim_oos <- do.call( cbind, lSim_oos[[ bt_name ]][ factor_names ] )
  sim_all <- na.omit( cbind( X_bm[ ,bm_names], sim_oos ) )
  lStats_all <-  descStats( sim_all, descStatsSpec(annualize = TRUE) )
  colors <- c("black", "brown", "darkgrey", "lightgrey", 
              "red", "steelblue", "lightgreen", "yellow", 
              "darkred", "darkblue", "darkgreen", "orange")
  
  riskReturnPlot( lStats = lStats_all, 
                  X_bm = X_bm[ ,bm_names],
                  colors = colors )
  
  t( lStats_all$stats[ ,bm_names] )
  
  
  
  # lStats_sim1 <- descStats( sim_1, descStatsSpec(annualize = TRUE) )
  
  # riskReturnPlot( lStats = lStats_all, 
  #                 lStats_sim1 = lStats_sim1, 
  #                 X_bm = X_bm[ ,strategy_names] )
  
  
  
  # --------------------------------------------------------------------------
  # BEAR BULL ANALYSIS
  # --------------------------------------------------------------------------
  
  BBSObj <- BBSRC$new()
  BBSObj$setCtrl()
  BBSObj$data$X_level <- cumulated(X_bm[ ,"Capw"], "discrete")
  
  bm_names <- c("Capw", "Multifactor", colnames(sim) )
  bt_names <- c("mfac", "mfac_val_qual_size_mom_ctrl", "mfac_size_stdz_ctrl")
  factor_names <- setNames( c("rnd", "momentum_stdz", "value", "quality", "size_stdz", "lowvola_stdz"),
                            c("random", "momentum", "value", "quality", "size", "lowvola") )
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
  # bm_names <- c("Capw", "Multifactor", colnames(sim) )
  # bt_name <- "mfac"
  # factor_names <- setNames( c("rnd", "momentum_stdz", "value", "quality", "size_stdz", "lowvola_stdz"),
  #                           c("random", "momentum", "value", "quality", "size", "lowvola") )
  # lPhase_stats <- list()
  # for ( factor_name in factor_names ) {
  #   BBSObj$data$X <- na.omit( cbind( X_bm[ ,bm_names],
  #                                    lSim_oos[[ bt_name ]][[ factor_name ]] ) )
  #   BBSObj$runRobust()
  #   lPhase_stats[[ factor_name ]] <- BBSObj$phaseStats()
  # }

  
  
  par( mfrow = c(2, 1) )
  colors <- fBasics::divPalette(n = ncol(lPhase_stats$mfac$rnd$`states==1`), "RdYlGn")
  colors_bm <-  c( c("black", "darkgrey", "lightgrey"), fBasics::rainbowPalette(n = ncol(sim)) )
  bt_name <- "mfac"
  # bt_name <- "mfac_val_qual_size_mom_ctrl"
  # bt_name <- "mfac_size_stdz_ctrl"  
  factor_name <- "momentum_stdz"
  plot( x = lPhase_stats[[ bt_name ]][[ factor_name ]]$`states==1`["sds", ], 
        y = lPhase_stats[[ bt_name ]][[ factor_name ]]$`states==1`["means", ], pch = 19, col = colors, xlab = "", ylab ="" )
  for ( i in seq(along = bm_names) ) {
    points( x = lPhase_stats[[ bt_name ]][[ factor_name ]]$`states==1`["sds", bm_names[i]], 
            y = lPhase_stats[[ bt_name ]][[ factor_name ]]$`states==1`["means", bm_names[i]], 
            pch = 19, cex = 2, col = colors_bm[i] )
  }
  plot( x = lPhase_stats[[ bt_name ]][[ factor_name ]]$`states==-1`["sds", ], 
        y = lPhase_stats[[ bt_name ]][[ factor_name ]]$`states==-1`["means", ], pch = 19, col = colors, xlab = "", ylab ="" )
  for ( i in seq(along = bm_names) ) {
    points( x = lPhase_stats[[ bt_name ]][[ factor_name ]]$`states==-1`["sds", bm_names[i]], 
            y = lPhase_stats[[ bt_name ]][[ factor_name ]]$`states==-1`["means", bm_names[i]], 
            pch = 19, cex = 2, col = colors_bm[i] )
  }
  legend( "topleft", bm_names, lwd = 2, col = colors_bm, text.col = colors_bm, bty = "n", cex = 0.7 )
  
    
  
  
  
  
  
  
  
  
  
  
  