
  
  ############################################################################
  ### RANDOM PORTFOLIOS - GEOMETRIC RANDOM WALK SAMPLING - USA
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
  # Initialize Backtest object
  # --------------------------------------------------------------------------
  
  # # Instantiate default Backtest class and prepare data for given universe
  # BT0 <- BacktestCustom$new()
  # 
  # # Default specifications for given universe
  # BT0$setCtrl( universe = universe,
  #              name = paste0(universe, "_maxw5"),
  #              selection_filter = "db_flag",
  #              # selection_filter = "bm",
  #              # cons_fun = c("box", "country", "sector", "reit", "esgWhenAvailable"),
  #              cons_fun = "box",
  #              width = 52 * 3,
  #              wd_data = paste0(wd, "Data/"),
  #              wd_warehouse = paste0(wd, "waRehouse/"),
  #              clean_RBE = FALSE )
  # 
  # # Default data for given universe
  # BT0$inputData()
  # BT0$data$upper_mat <- BT0$data$upper_mat * 0 + 0.05
  # BT0$spec$rebdates <- BT0$spec$rebdates[ which(BT0$spec$rebdates < tail(rownames(BT0$data$X_sim), 1) ) ]
  # 
  # # Save default backtest object
  # BT0$save( wd = BT0$spec$wd_data,
  #           name = BT0$spec$name,
  #           without_data = FALSE )
  
  
  BT0 <- loadBacktest( wd = paste0(wd_data, "Data/"),
                       name = "usa_maxw5" )
  BT0$initRBE()
  
  
  
  
  # --------------------------------------------------------------------------
  # BACKTEST UNIFORM RP'S WITHIN CONSTRAINTS
  # --------------------------------------------------------------------------
  
  BT1 <- BacktestCustom$new()
  BT1$data <- BT0$data
  BT1$spec <- BT0$spec
  BT1$spec$keep_all_rebdates <- FALSE
  BT1$spec$portfolio <- "grw"
  # debugonce( BT1$grwPortfolio )
  BT1$runLoop()  
  
  # Simulate
  tmp <- lapply( BT1$output, FUN = simFUN )
  sim_1 <- do.call( cbind, tmp )
  sim_1 <- sim_1[isWeekday(time(sim_1)), ]
  colnames(sim_1) <- paste0("rp_grw_uniform", 1:ncol(sim_1))
  
  # Apply FEV-bias
  lWeights <- BT1$output
  lWeights_fev <- lapply( lWeights, FUN = function(wmat) { fevBias( x = wmat, q = 4 ) } )
  
  # Project FEV-biased weights to boundary of feasible set
  lWeights_cons <- lWeights
  for ( i in 1:length(lWeights_cons) ) {
    for ( today in rownames(lWeights_cons[[i]]) ) {
      RBE <- rebalenv( rebdate = today )
      # Cons <- getConstraints(RBE$GPS)
      lincon <- getConstraints( RBE$GPS, "bounds" )
      upper_bounds <- lincon$upper
      w_unc <- lWeights_fev[[i]][today, ]
      # idx_notNA <- which(!is.na(w_unc))
      idx_notNA <- which( !is.na( lWeights[[i]][today, ] ) )
      w_unc <- w_unc[idx_notNA]
      w_unc <- setNames( as.numeric(w_unc), names(w_unc))
      # debugonce( map2boxcon )
      lWeights_cons[[i]][today, idx_notNA] <- map2boxcon( w_unc = w_unc, 
                                                          upper = upper_bounds, 
                                                          itermax = 10^3 )
    }
  }
  
  # Simulate
  tmp <- lapply( lWeights_cons, FUN = simFUN )
  sim_1_cons <- do.call( cbind, tmp )
  sim_1_cons <- sim_1_cons[isWeekday(time(sim_1_cons)), ]
  colnames(sim_1_cons) <- paste0("rp_grw_uniform_cons", 1:ncol(sim_1_cons))
 
  
  
  lStats_1 <- descStats( sim_1, descStatsSpec(annualize = TRUE) )
  stats_1 <- lStats_1$stats
  lStats_1_cons <- descStats( sim_1_cons, descStatsSpec(annualize = TRUE) )
  stats_1_cons <- lStats_1_cons$stats  
  
  plot( x = stats_1_cons["sds", ], y = stats_1_cons["means", ], pch = 19, col = 2 )
  points( x = stats_1["sds", ], y = stats_1["means", ], pch = 19, col = 1 )
  
  
  
  
  
  
  # Combine
  end_date <- "2022-03-06"
  strategy_names <- colnames(X_bm)
  sim_all <- na.omit( cbind( X_bm[ ,strategy_names], sim_1_cons) )
  sim_all <- window( sim_all, start(sim_all), end_date )
  lStats_all <-  descStats( sim_all, descStatsSpec(annualize = TRUE) )
  
  
  lStats_sim1 <- descStats( window( sim_1, start(sim_1), end_date), descStatsSpec(annualize = TRUE) )
  
  
  
  
  
  
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
  rownames(df) <- colnames(sim_all)
  h1 <- hist(df$x, breaks=250, plot=F)
  h2 <- hist(df$y, breaks=250, plot=F)
  top <- max(h1$counts, h2$counts) * 2
  lims_x <- range(c(df$x * 0.995, df$x * 1.005))
  # lims_x <- range(df$x)
  lims_y <- range(df$y)
  k <- kde2d( df$x, df$y, n = 250, lims = c(lims_x, lims_y) )
  # margins
  oldpar <- par()
  par(mar=c(3,3,1,1))
  layout(matrix(c(2,0,1,3), 2, 2, byrow = TRUE), c(3,1), c(1,3))
  image(k)
  idx_x <- c( which(df$x > quantile(df$x, 0.9)), which(df$x < quantile(df$x, 0.1)) )
  idx_y <- c( which(df$y > quantile(df$y, 0.9)), which(df$y < quantile(df$y, 0.1)) )
  idx <- c( idx_x, idx_y )
  points( x = df$x[idx], y = df$y[idx], pch = 19, col = "orange", cex = 0.5 )
  # tmp <- MASS::kde2d( x = df$x, y = df$y, n = 25, h = c(0.005, 0.02), lims = c(range(df$x), range(df$y)) ) 
  # contour( tmp, xlab = "Portfolio Standard Deviation", ylab = "Portfolio Return",
  #          add = TRUE, drawlabels = FALSE, col = "red" )
  df1 <- data.frame( x = lStats_sim1$stats["sds", ],
                     y = lStats_sim1$stats["means", ] )
  k1 <- kde2d( df1$x, df1$y, n = 250 )
  tmp <- MASS::kde2d( x = df1$x+0.001, y = df1$y, n = 25, lims = c(range(df1$x+0.001), range(df1$y)) ) 
  contour( tmp, xlab = "Portfolio Standard Deviation", ylab = "Portfolio Return",
           add = TRUE, drawlabels = FALSE, col = "blue" )
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
  
  ###
  
  
  
  
  
  
  
  
  