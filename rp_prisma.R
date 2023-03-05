
  
  ############################################################################
  ### RANDOM PORTFOLIOS - RANDOM PORTOLIOS WITHIN UPPER BOUNDS - AKTIEN CH
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     13.11.2022
  # First version:    13.11.2022
  # --------------------------------------------------------------------------
  
  
  
  
  
  
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
  wd <- "R:/Asset_Management/Research_Projects/External_Research_Projects/GeomScale/Random_Portfolios/"
  source( paste0(wd, "Source/class_BacktestCustom.R") )
  source( paste0(wd, "Source/custom_functions.R") )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # FUNCTIONS
  # --------------------------------------------------------------------------
  
  simFUN <- function(wmat) 
  { 
    simPortfolio( X = BT0$data$X_sim, 
                  wghts = wmat,
                  language = "C" )
  }
  
  
  
  
  # --------------------------------------------------------------------------
  # PARAMETERS
  # --------------------------------------------------------------------------
  
  # universe <- "aktien_ch"
  universe <- "aktien_ch_fub"
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  
  
  # --------------------------------------------------------------------------
  # LOAD DATA AND INITIALIZE REBALANCING ENVIRONMENTS
  # --------------------------------------------------------------------------
  
  BT0 <- loadBacktest( wd = paste0(wd, "Data/"),
                       name = universe )
  ###
  # Quick fix
  rebdates <- BT0$spec$rebdates
  # rebdates <- rebdates[ rebdates >= "2011-02-23" ]  # if adaptive ubber bounds are used
  BT0$spec$rebdates <- rebdates
  BT0$data$upper_mat <- BT0$data$upper_mat * 1.25
  ###
  BT0$initRBE()
  
  
  
 
  # --------------------------------------------------------------------------
  # LOAD OLZ FUND SERIES BENCHMARKS
  # --------------------------------------------------------------------------
  
  # Benchmark series
  tickers <- setNames( c("OLZSOIR SW", "SPI", "SMI", "SMIMC", "SLIC", "SPIEX"), #, ".SPI_EXB3", "VPBSEBI LE"),
                       c("OLZ", "SPI", "SMI", "SMIM", "SLI", "SPI_Extra") ) # , "SPI_EXB3", "VPB") )
  X_bm <- rodbcGetOLZDBReturns( assetName = tickers,
                                  refCcy = "CHF",
                                  frqncy = "daily",
                                  na.rm = "r" )
  X_bm <- X_bm[isWeekday(time(X_bm)), tickers]
  colnames(X_bm) <- names(tickers)
  
  
  
  
  
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
    RBE$GPO <- GPO
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
  

  
  to <- unlist( lapply( BT_qcapw$spec$rebdates, FUN = function(date) { rebalenv(date)$GPO@turnover } ) )
  to
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BACKTEST QUASI EQW PORTFOLIO
  # --------------------------------------------------------------------------
  
  BT_qeqw <- BacktestCustom$new()
  BT_qeqw$data <- BT0$data
  BT_qeqw$spec <- BT0$spec
  BT_qeqw$spec$portfolio <- NULL
  BT_qeqw$spec$GPS@covariance$method <- "duv"
  BT_qeqw$spec$selection_filter <- c("db_flag", "upperbound")
  BT_qeqw$spec$cons_fun <- "box"
  # BT_qeqw$spec$selection_filter <- "bm"
  # BT_qeqw$spec$cons_fun <- NULL
  # debugonce( BT_qeqw$setSelection )
  BT_qeqw$runLoop()
  sim_qeqw <- BT_qeqw$simulate( fc = 0, vc = 0 )
  
  
  # Add qeqw and qcapw to X_bm
  X_bm <- na.omit( cbind( X_bm, QEQW = sim_qeqw, QCapw = sim_qcapw ) )
  
  
  plotSimTS( as.simTS( X_bm ) )
  lStats <- descStats( window( X_bm, "2021-12-31", "2022-12-31" ) )
  lStats
  # debugonce( statsolz:::plot.stats )
  statsolz:::plot.stats( lStats, sortby = NULL, scl = 0 )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # RANDOM PORTFOLIOS - UNIFORM
  # --------------------------------------------------------------------------
  
  BT2 <- BacktestCustom$new()
  BT2$data <- BT0$data
  BT2$spec <- BT0$spec
  BT2$spec$portfolio <- "momentum_rp"
  BT2$spec$n_sim <- 10^2 * 1
  BT2$spec$sampling_dist = "uniform"
  # BT2$spec$m <- NULL # means m = length(selection)
  BT2$spec$m <- 30
  BT2$spec$shadow_dirichlet_transform <- FALSE
  BT2$spec$scl_by_capw <- FALSE
  BT2$runLoop()  
  
  # Simulate
  tmp <- lapply( BT2$output, FUN = simFUN )
  sim_2 <- do.call( cbind, tmp )
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
  BT3$spec$m <- NULL # means m = length(selection)
  # BT3$spec$m <- 30
  # BT3$spec$th <- 0.1
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
  # sim_all <- na.omit( cbind( X_bm, sim_2, sim_3, sim_2_cons, sim_3_cons) )
  # sim_all <- na.omit( cbind( X_bm, sim_2_cons, sim_3_cons) )
  sim_all <- na.omit( cbind( X_bm, sim_3_cons) )
  lStats <-  descStats( sim_all, descStatsSpec(annualize = TRUE) )
  
  # statsolz:::plot.stats( lStats, sortby = NULL )
  
  colors <- fBasics::rainbowPalette(n = ncol(X_bm))
  plot( x = lStats$stats["sds", ], y = lStats$stats["means", ], pch = 19, col = "grey" )
  points( x = lStats$stats["sds", colnames(X_bm)], y = lStats$stats["means", colnames(X_bm)], pch = 19, col = colors, cex = 2 )
  points( x = lStats$stats["sds", colnames(sim_2_cons)], y = lStats$stats["means", colnames(sim_2_cons)], pch = 8, col = "orange")
  points( x = lStats$stats["sds", colnames(sim_3_cons)], y = lStats$stats["means", colnames(sim_3_cons)], pch = 19, col = "orange")
  legend("bottomleft", colnames(X_bm), lwd = 2, col = colors, text.col = colors, bty = "n" )
  

  lStats$desc
  t( lStats$stats[stats_fields, 1:10] )  
  
  
  
  
  # --------------------------------------------------------------------------
  # PERFORMANCE ANALYSIS
  # --------------------------------------------------------------------------
  
  
  # Year to date
  start_date <- "2022-01-01"
  
  # plot( as.simTS( X_bm[which(rownames(X_bm) > start_date), ] ) )
  
  sim_all_ytd <- sim_all[which(rownames(sim_all) >= start_date), ]
  lStats_ytd <-  descStats( sim_all_ytd, descStatsSpec(annualize = TRUE) )

  colors <- fBasics::rainbowPalette(n = ncol(X_bm))
  plot( x = lStats_ytd$stats["sds", ], y = lStats_ytd$stats["means", ], pch = 19, col = "grey" )
  points( x = lStats_ytd$stats["sds", colnames(X_bm)], y = lStats_ytd$stats["means", colnames(X_bm)], pch = 19, col = colors, cex = 2 )
  points( x = lStats_ytd$stats["sds", colnames(sim_2_cons)], y = lStats_ytd$stats["means", colnames(sim_2_cons)], pch = 8, col = "orange")
  points( x = lStats_ytd$stats["sds", colnames(sim_3_cons)], y = lStats_ytd$stats["means", colnames(sim_3_cons)], pch = 19, col = "orange")
  legend("bottomleft", colnames(X_bm), lwd = 2, col = colors, text.col = colors, bty = "n" )
  
  
  
  t(lStats_ytd$stats[stats_fields, 1:5])
  
  
  
  
  
  ###
  
  df <- data.frame( x = lStats$stats["sds", ],
                    y = lStats$stats["means", ] )
  df <- data.frame( x = lStats_ytd$stats["sds", ],
                    y = lStats_ytd$stats["means", ] )
  h1 <- hist(df$x, breaks = 150, plot = FALSE)
  h2 <- hist(df$y, breaks = 150, plot = FALSE)
  top <- max(h1$counts, h2$counts)
  k <- kde2d( df$x, df$y, n = 300 )
  # margins
  oldpar <- par()
  par(mar=c(3,3,1,1))
  layout( matrix(c(2,0,1,3),2,2, byrow = TRUE), c(3,1), c(1,3))
  image(k)
  # image(k, col = viridis_pal()(200) )
  points( x = df$x, y = df$y, pch = 19, col = "orange", cex = 0.5 )
  tmp <- MASS::kde2d( x = df$x, y = df$y, n = 25, h = c(0.005, 0.02), lims = c(range(df$x), range(df$y)) ) 
  # contour( tmp, xlab = "Portfolio Standard Deviation", ylab = "Portfolio Return", 
  #          add = TRUE, drawlabels = FALSE, col = "red" )
  colors <- fBasics::rainbowPalette(n = ncol(X_bm))
  for ( i in 1:ncol(X_bm) ) {
    points( x = df$x[i], y = df$y[i], pch = 19, cex = 2, col = colors[i] )
  }
  legend( "topright", colnames(X_bm), lwd = 2, col = colors, text.col = colors, bty = "n")
  par(mar=c(0,2,1,0))
  barplot(h1$counts, axes=FALSE, ylim=c(0, top), space=0, col='grey')
  par(mar=c(2,0,0.5,1))
  barplot(h2$counts, axes=FALSE, xlim=c(0, top), space=0, col='grey', horiz=TRUE)
  
  ###
  
  
  
  
  
  
  
  require(ggplot2)
  require(patchwork)
  require(ggpubr)
  
  
  make_xyz <- function(x, y, n=200)
  {
    z <- MASS::kde2d(x, y, n=n)
    d <- matrix(nrow=n*n, ncol=3)
    k = 1
    for (i in 1:n){
      for (j in 1:n){
        d[k,] <- c(z$x[i], z$y[j], z$z[i,j])
        k <- k + 1
      }
    }
    colnames(d) <- c("x", "y", "z")
    d <- as_tibble(d)
    d
  }
  d <- make_xyz(df$x, df$y, 200)
  plt <- ggplot(d) +
    geom_raster(aes(x=x, y=y, fill=z)) +
    coord_fixed() +  
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position='none') +
    scale_fill_viridis()
  
  plt
  
  
  
  
  plt + 
    geom_point(aes(x, y), alpha=0.5, pch=1, data = df, col="white", size=0.5) + 
    geom_point(aes(x, y), alpha=0.5, pch=19, data = df[1:5, ], col = "black", size = 3) + 
    geom_density_2d(aes(x,y), data=xy, alpha=0.2, col="white")
  
  
  
  
  
  
  dens1 <- ggplot(xy, aes(x = x)) + 
    geom_histogram(color = "black", fill = "white") + 
    theme_void()
  
  dens2 <- ggplot(xy, aes(x = y)) + 
    geom_histogram(color = "black", fill = "white") + 
    theme_void() + 
    coord_flip()
  
  dens1 + plot_spacer() + plt + dens2 + 
    plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
  
  
  dens1 + plt + dens2 + 
    plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # ROLLING PERFORMANCE
  # --------------------------------------------------------------------------
  
  
  mu_roll <- applyRoll( Data = sim_all,
                        Width = 252 * 1,
                        By = 1,
                        FUN = meanGeo, 
                        scalefactor = 252 )
  sds_roll <- applyRoll( Data = sim_all,
                         Width = 252 * 1,
                         By = 1,
                         FUN = function(X) { apply(X, 2, sd) * sqrt(252) } )
  
  
  # 3-d plot
  
  
  
  
  # --------------------------------------------------------------------------
  # RELATIVE PERFORMANCE
  # --------------------------------------------------------------------------
  
  
  # Relative performance
  
  lSim_delta <- lapply( (ncol(X_bm)+1):ncol(sim_all), FUN = function(i) { simOutperformance( x = sim_all[ ,"OLZ"], 
                                                                                             y = sim_all[ ,i] ) } )
  sim_delta <- do.call( cbind, lSim_delta )
  
  plot( log(cumulated(sim_delta, "discrete")), plot.type = "single", col = "orange" )
  abline( h = 0, col = "grey", lwd = 2 )
  
  
  # Rolling relative performance
  
  stats_roll <- applyRoll( Data = cumulated(sim_delta, "discrete"),
                           Width = 252 * 1,
                           By = 1,
                           FUN = function(x) { as.numeric(x[nrow(x), ]) / as.numeric(x[1, ]) - 1 } )
  
  plot( stats_roll, plot.type = "single", col = "orange" )
  abline( h = 0, col = "grey" )
  
  
  perc_outperf <- apply( stats_roll, 1, function(x) { sum(x > 0) / length(x) } )
  plot( perc_outperf )
  
  
  
  # Market regime
  
  BBSObj <- BBSRC$new()
  BBSObj$setCtrl()
  BBSObj$data$X_level <- cumulated(X_bm[ ,"SPI"], "discrete")
  BBSObj$runRobust()
  states <- BBSObj$output$states
  
  Y <- BBSObj$data$X_level
  plot(Y, main = "")
  abline(v = time(states)[which(states == 1)], col = "green")
  abline(v = time(states)[which(states == -1)], col = "tomato")
  lines(Y, col = "blue")
  
  
  
  BBSObj$data$X <- sim_delta
  phase_stats <- BBSObj$phaseStats()
  
  boxplot( as.data.frame( cbind( s1 = phase_stats$`states==1`["means", ],
                                 s2 = phase_stats$`states==-1`["means", ] ) ), 
           beside = TRUE, col = c("green", "red") )
  
  
  
  
  
  
  
  # Calmar ratio
  calmarRatio( X = X_bm, scalefactor = 252 )
  cr <- calmarRatio( X = sim_all, scalefactor = 252 )
  
  plot( density( cr ) )
  plot( density( meanGeo(sim_all, scalefactor = 252) ))
  
  
  # Relative Calmar ratio
  cr_rel <- calmarRatio( X = sim_delta, scalefactor = 252 )
  
  plot( density( cr_rel ) ) 
  abline( v = 0 )
  
  
  # Rolling relative Calmar ratio
  cr_rel_roll <- applyRoll( Data = sim_delta,
                            Width = 252 * 1,
                            By = 1,
                            FUN = calmarRatio,
                            scalefactor = 252 )
  
  
  headleft(cr_rel_roll)
  plot( cr_rel_roll, plot.type = "single" )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # SAVE
  # --------------------------------------------------------------------------

  env <- new.env()
  env$sim_all <- sim_all
  env$lStats <- lStats
  env$lWeights_eqw <- lWeights_eqw
  env$lWeights_capw <- lWeights_capw
  saveRDS( env, file = paste0(wd, "waRehouse/aktien_ch_rp.rds") )





  # --------------------------------------------------------------------------
  # LOAD
  # --------------------------------------------------------------------------

  env <- readRDS( file = paste0(wd, "waRehouse/aktien_ch_rp.rds") )
  sim_all <- env$sim_all
  lStats <- env$lStats
  lWeights <- env$lWeights

  
  
  
  
  
 
  
  
  # --------------------------------------------------------------------------
  # SIMULATIONS OF LETTER-SCORE PORTFOLIOS
  # --------------------------------------------------------------------------
  
  # Full period
  sim_tmp_full <- do.call( cbind, lSim )
  sim_tmp <- na.omit( cbind( bm = BT0$data$X_bm, 
                             eqw = sim_eqw, 
                             sim_tmp_full  ) )
  lStats <- descStats( sim_tmp )
  
    
  # First decade
  sim_tmp_decade1 <- window( sim_tmp, start(sim_tmp), "2011-01-01" )
  lStats_decade1 <- descStats( sim_tmp_decade1 )
  
  sim_tmp_decade1_inv <- sim_tmp_inv[rownames(sim_tmp_decade1), ]
  lStats_decade1_inv <- descStats( sim_tmp_decade1_inv )
 
  # Second decade
  sim_tmp_decade2 <- window( sim_tmp, "2011-01-01", end(sim_tmp) )
  lStats_decade2 <- descStats( sim_tmp_decade2 )
  
  sim_tmp_decade2_inv <- sim_tmp_inv[rownames(sim_tmp_decade2), ]
  lStats_decade2_inv <- descStats( sim_tmp_decade2_inv )
 
  
  
  # Plot risk-return
  
  colors <- c(rep("grey", length(letter_vec)), 2, 4)
  
  statsolz:::plot.stats( lStats, sortby = NULL, col = colors, legend.loc = NULL )
  length( which( lStats$stats["means", ] >  lStats$stats["means", "bm"]) ) / 26
  
  statsolz:::plot.stats( lStats_inv, sortby = NULL, col = colors, legend.loc = NULL )
  length( which( lStats_inv$stats["means", ] >  lStats_inv$stats["means", "bm"]) ) / 26
  
  statsolz:::plot.stats( lStats_decade1, sortby = NULL, col = colors, legend.loc = NULL )
  length( which( lStats_decade1$stats["means", ] >  lStats_decade1$stats["means", "bm"]) ) / 26
  
  statsolz:::plot.stats( lStats_decade1_inv, sortby = NULL, col = colors, legend.loc = NULL )
  length( which( lStats_decade1_inv$stats["means", ] >  lStats_decade1_inv$stats["means", "bm"]) ) / 26
  
  statsolz:::plot.stats( lStats_decade2, sortby = NULL, col = colors, legend.loc = NULL )
  length( which( lStats_decade2$stats["means", ] >  lStats_decade2$stats["means", "bm"]) ) / 26
  
  statsolz:::plot.stats( lStats_decade2_inv, sortby = NULL, col = colors, legend.loc = NULL )
  length( which( lStats_decade2_inv$stats["means", ] >  lStats_decade2_inv$stats["means", "bm"]) ) / 26
  
  
  
  # colors <- fBasics::divPalette(n = ncol(sim_sort), "RdYlGn" )
  # plot( as.simTS( sim_tmp ), col = c(colors, 1, 4) )
  # barplot( lStats$stats["sharpe", ], col = c(colors, 1, 4) )
  # barplot( lStats$stats["cumret", ], col = c(colors, 1, 4) )
  # abline( h = min(lStats$stats["cumret", ]) )
  # 
  # 
  # weightsBarPlot( lWeights[["Z"]] )
  # weightsBarPlot( lWeights[["A"]] )
  
  
  # Outperformance
  sim_outperf <- simOutperformance( x = sim_tmp[ ,"Z"], y = sim_tmp[ ,"bm"] )
  sim_outperf_decade1 <- sim_outperf[rownames(sim_tmp_decade1), ]
  sim_outperf_decade2 <- sim_outperf[rownames(sim_tmp_decade2), ]
  
  plot( as.simTS(sim_outperf) )
  descStats( sim_outperf )
  
  
  # Sector allocation
  sector_alloc <- groupWeights( wmat = lWeights[["Z"]], group_mat = BT0$data$sector_mat )
  weightsBarPlot( sector_alloc )
  boxplot( as.data.frame(sector_alloc) )
  
  
  # Sharpe test
  X <- sim_tmp[ ,c("Z", "bm", "eqw")]
  ST <- sharpeTest( X = X, method = "HAC_PW" )
  ST
 
  
  
  
  # --------------------------------------------------------------------------
  # RISK-RETURN DISTRIBTION CONTOUR PLOT
  # --------------------------------------------------------------------------
  
  # checkout: https://www.r-bloggers.com/2014/09/5-ways-to-do-2d-histograms-in-r/
  
  # Or
  
  # R package LogConcDEAD
  
  
  
  df <- data.frame( x = lStats$stats["sds", ],
                    y = lStats$stats["means", ] )
 
  tmp <- MASS::kde2d( x = df$x, y = df$y, n = 25, lims = c(range(df$x), range(df$y)) )
  statsolz:::plot.stats( lStats, sortby = NULL, col = "grey", txt.color = NULL, legend.loc = NULL, cex = 1 )
  contour( tmp, xlab = "Portfolio Standard Deviation", ylab = "Portfolio Return", 
           add = TRUE, drawlabels = FALSE, col = "red" )
  points( x = df$x[1], y = df$y[1], pch = 19, cex = 2 )
  points( x = df$x[2], y = df$y[2], pch = 19, cex = 2, col = "blue" )
  points( x = df$x[which(grepl("A.", rownames(df)))], y = df$y[which(grepl("A.", rownames(df)))], pch = 19, col = "brown" )
  
  
  
  
  tmp <- MASS::kde2d( x = df$x, y = df$y, n = 25, lims = c(range(df$x), range(df$y)) )
  image( tmp )
  points( x = df$x, y = df$y, col = "grey" )
  
  
  
  
  
  h1 <- hist(df$x, breaks=250, plot=F)
  h2 <- hist(df$y, breaks=250, plot=F)
  top <- max(h1$counts, h2$counts)
  # margins
  oldpar <- par()
  par( mar = c(3, 3, 1, 1) )
  layout( matrix( c(2, 0, 1, 3), 2, 2, byrow = TRUE), c(3, 1), c(1, 3) ) 
  tmp <- MASS::kde2d( x = df$x, y = df$y, n = 15, lims = c(range(df$x), range(df$y)) )
  statsolz:::plot.stats( lStats, sortby = NULL, col = "grey", 
                         txt.color = NULL, legend.loc = NULL, cex = 1 )
  contour( tmp, xlab = "Portfolio Standard Deviation", ylab = "Portfolio Return", 
           add = TRUE, drawlabels = FALSE, col = "red" )
  points( x = df$x[1], y = df$y[1], pch = 19, cex = 2 )
  points( x = df$y[2], y = df$y[2], pch = 19, cex = 2, col = "blue" )
  par( mar = c(0, 2, 1, 0) )
  barplot( h1$counts, axes = FALSE, ylim = c(0, top), space = 0, col = 'grey' )
  par( mar = c(2, 0, 0.5, 1) )
  barplot( h2$counts, axes = FALSE, xlim = c(0, top), space = 0, col = 'grey', horiz = TRUE )
  
  
  
  ###
  
  df <- data.frame( x = lStats$stats["sds", ],
                    y = lStats$stats["means", ] )
  h1 <- hist(df$x, breaks=250, plot=F)
  h2 <- hist(df$y, breaks=250, plot=F)
  top <- max(h1$counts, h2$counts)
  k <- kde2d( df$x, df$y, n = 100 )
  # margins
  oldpar <- par()
  par(mar=c(3,3,1,1))
  layout(matrix(c(2,0,1,3),2,2,byrow=TRUE),c(3,1), c(1,3))
  image(k)
  abline( h = df$y[1], col = "white", lwd = 3 )
  abline( v = df$x[1], col = "white", lwd = 3 )
  abline( h = df$y[2], col = "white", lwd = 3 )
  abline( v = df$x[2], col = "white", lwd = 3 )
  points( x = df$x, y = df$y, pch = 19, col = "grey" )
  tmp <- MASS::kde2d( x = df$x, y = df$y, n = 25, h = c(0.005, 0.02), lims = c(range(df$x), range(df$y)) ) 
  contour( tmp, xlab = "Portfolio Standard Deviation", ylab = "Portfolio Return", 
           add = TRUE, drawlabels = FALSE, col = "red" )
  points( x = df$x[1], y = df$y[1], pch = 19, cex = 2 )
  points( x = df$x[2], y = df$y[2], pch = 19, cex = 2, col = "blue" )
  par(mar=c(0,2,1,0))
  barplot(h1$counts, axes=FALSE, ylim=c(0, top), space=0, col='grey')
  par(mar=c(2,0,0.5,1))
  barplot(h2$counts, axes=FALSE, xlim=c(0, top), space=0, col='grey', horiz=TRUE)
  
  ###
  
 
  
  
  
  
  
  
  
  # Using LogConcDEAD
  # install.packages("LogConcDEAD")
  require(LogConcDEAD)
  
  
  tic <- Sys.time()
  A <- as.matrix( t( lStats$stats[c("sds", "means"), ] ) )
  LCD <- mlelcd( x = A )
  (toc <- Sys.time() - tic)
  
  
  plot( LCD, addp = TRUE, uselog = TRUE, type = "ic", drawlabels = FALSE )
  points( x = lStats$stats["sds", (ncol(lStats$stats)-27):(ncol(lStats$stats)-1)],
          y = lStats$stats["means", (ncol(lStats$stats)-27):(ncol(lStats$stats)-1)],
          col = "grey", cex = 2, pch = 19 )
  points( x = lStats$stats["sds", 2],
          y = lStats$stats["means", 2], col = "blue", cex = 2, pch = 19 )
  points( x = lStats$stats["sds", 1],
          y = lStats$stats["means", 1], col = "red", cex = 2, pch = 19 )
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # FACTOR SCORES
  # --------------------------------------------------------------------------
 
  # roll_factor <- function(wmat, factor_mat)
  # {
  #   FUN <- function ( X, factor_mat ) 
  #   {
  #     selection <- intersect(colnames(X), colnames(factor_mat))
  #     zeit <- intersect(as.character(time(X)), as.character(time(factor_mat)))
  #     factor_mat <- factor_mat[zeit, selection]
  #     X <- X[zeit, selection]
  #     X[is.na(factor_mat)] <- NA
  #     X <- X / rowSums(X, na.rm = TRUE)
  #     factor_mat <- factor_mat * (!is.na(X))
  #     ans <- as.timeSeries(rowSums(X * factor_mat, na.rm = TRUE))
  #     return(ans)
  #   }
  #   ans <- do.call(cbind, lapply(X = wmat, FUN = FUN, factor_mat = factor_mat))
  #   
  #   return(ans)
  # }
  
  
  wmat <- lWeights[[1]]
  score_mat <- BT0$data$wmat_bm

  # debugonce( portfolioScore )
  pfs <- portfolioScore( wmat = wmat, score_mat = score_mat, b_rescale = TRUE, agg_method = "arithmetic" )
  pfs_h <- portfolioScore( wmat = wmat, score_mat = score_mat, b_rescale = TRUE, agg_method = "harmonic" )
  
  tmp <- cbind( pfs, pfs_h )
  plot( tmp )
  
  
  
  
  # Load factor matrices from factor Shiny tool
  wd <- "R:/Asset_Management/R/Shiny/Factor_Portfolios/"
  # global_data <- readRDS( file = paste0(wd, "Data/lOutput.rds") )
  env <- readRDS( file = paste0(wd, "Data/", universe, ".rds"))
  ls(env)
  ls(env$data_cons)
  
  
  # Benchmark portfolio
  wmat_bm <- BT0$data$wmat_bm
  
  
  # Size
  # size_score_mat <- BT0$data$wmat_bm
  size_score_mat <- env$data_cons$lowsize_mat
  lSizeScore  <- portfolioScore( wmat = lWeights,
                           score_mat = size_score_mat,
                           b_rescale = TRUE,
                           agg_method = "arithmetic" )
  size_score <- do.call( cbind, lSizeScore )
  size_score_bm <- portfolioScore( wmat = wmat_bm,
                                   score_mat = size_score_mat,
                                   b_rescale = TRUE,
                                   agg_method = "arithmetic" )
  
  colors <- fBasics::divPalette(n = length(lWeights), "RdYlGn")
  plot( na.omit(cbind(size_score, size_score_bm)), plot.type = "single", col = colors )
  lines( size_score_bm[rownames(size_score), ], col = 1 )
  
  
  
  # Value
  # value_score_mat <- BT0$data$value_mat
  value_score_mat <- env$data_cons$value_mat
  lValueScore  <- portfolioScore( wmat = lWeights,
                                  score_mat = value_score_mat,
                                  b_rescale = TRUE,
                                  # agg_method = "harmonic" )
                                  agg_method = "arithmetic" )
  value_score <- do.call( cbind, lValueScore )
  value_score_bm <- portfolioScore( wmat = wmat_bm,
                                    score_mat = value_score_mat,
                                    b_rescale = TRUE,
                                    # agg_method = "harmonic" )
                                    agg_method = "arithmetic" )
  
  colors <- fBasics::divPalette(n = length(lWeights), "RdYlGn")
  plot( na.omit(cbind(value_score, value_score_bm)), plot.type = "single", col = colors )
  lines( value_score_bm[rownames(value_score), ], col = 1 )
  

  
  # Quality
  # quality_score_mat <- BT0$data$quality_mat
  quality_score_mat <- env$data_cons$quality_mat
  lQualityScore  <- portfolioScore( wmat = lWeights,
                                  score_mat = quality_score_mat,
                                  b_rescale = TRUE,
                                  # agg_method = "harmonic" )
                                  agg_method = "arithmetic" )
  quality_score <- do.call( cbind, lQualityScore )
  quality_score_bm <- portfolioScore( wmat = wmat_bm,
                                    score_mat = quality_score_mat,
                                    b_rescale = TRUE,
                                    # agg_method = "harmonic" )
                                    agg_method = "arithmetic" )
  
  
  # ROE
  roe_score_mat <- makeMTSObj( Data = env$data_index$X,
                               Head.Name = "Series_Id",
                               Head.Ret = "z_roe",
                               Head.Date = "Date" )
  lROEScore  <- portfolioScore( wmat = lWeights,
                                    score_mat = quality_score_mat,
                                    b_rescale = TRUE,
                                    # agg_method = "harmonic" )
                                    agg_method = "arithmetic" )
  roe_score <- do.call( cbind, lROEScore )
  roe_score_bm <- portfolioScore( wmat = wmat_bm,
                                  score_mat = roe_score_mat,
                                  b_rescale = TRUE,
                                  # agg_method = "harmonic" )
                                  agg_method = "arithmetic" )
  
  # Momentum
  # momentum_score_mat <- BT0$data$momentum_mat
  momentum_score_mat <- env$data_cons$momentum_mat
  lMomentumScore  <- portfolioScore( wmat = lWeights,
                                     score_mat = momentum_score_mat,
                                     b_rescale = TRUE,
                                     # agg_method = "harmonic" )
                                     agg_method = "arithmetic" )
  momentum_score <- do.call( cbind, lMomentumScore )
  momentum_score_bm <- portfolioScore( wmat = wmat_bm,
                                       score_mat = momentum_score_mat,
                                       b_rescale = TRUE,
                                       # agg_method = "harmonic" )
                                       agg_method = "arithmetic" )
  
  
  
  
  # Volatility
  # vola_score_mat <- BT0$data$vola_mat
  # vola_score_mat <- env$data_cons$lowvola_mat
  vola_score_mat <- env$data_cons$lowvola_stdz_mat
  lVolaScore  <- portfolioScore( wmat = lWeights,
                                 score_mat = vola_score_mat,
                                 b_rescale = TRUE,
                                 # agg_method = "harmonic" )
                                 agg_method = "arithmetic" )
  vola_score <- do.call( cbind, lVolaScore )
  vola_score_bm <- portfolioScore( wmat = wmat_bm,
                                   score_mat = vola_score_mat,
                                   b_rescale = TRUE,
                                   # agg_method = "harmonic" )
                                   agg_method = "arithmetic" )
  
  
  
  xxx_score <- quality_score
  xxx_score_bm <- quality_score_bm
  xxx_score <- momentum_score
  xxx_score_bm <- momentum_score_bm
  xxx_score <- vola_score
  xxx_score_bm <- vola_score_bm
  xxx_score <- roe_score
  xxx_score_bm <- roe_score_bm
  
  
  colors <- fBasics::divPalette(n = length(lWeights), "RdYlGn")
  plot( na.omit(cbind(xxx_score, xxx_score_bm)), plot.type = "single", col = colors )
  lines( xxx_score_bm[rownames(xxx_score), ], col = 1 )
  
  
  
  
  
  
  
  
  
  
  
    
  
  # --------------------------------------------------------------------------
  # FF REGRESSION
  # --------------------------------------------------------------------------
  
  path_ff <- "R:/Asset_Management/R/myRData/Factor/"
  env_ff <- readRDS( file = paste0(path_ff, "data.Rds") )
  FF3 <- env_ff$lFactor$`3F`$USA_daily
  FF5 <- env_ff$lFactor$`5F`$USA_daily
  
  dates <- intersect(rownames(FF3), rownames(sim_tmp))
  X_train <- FF3[dates, 1:3]
  # dates <- intersect(rownames(FF5), rownames(sim_tmp))
  # X_train <- FF5[dates, 1:5]
  lReg <- lReg_inv <- list()
  lCoeff <- lCoeff_inv <- list()
  
  dates_decade1 <- intersect(rownames(FF5), rownames(sim_tmp_decade1))
  X_train_decade1 <- FF5[dates_decade1, 1:5]
  lReg_decade1 <- lReg_decade1_inv <- list()
  lCoeff_decade1 <- lCoeff_decade1_inv <- list()
  
  dates_decade2 <- intersect(rownames(FF5), rownames(sim_tmp_decade2))
  X_train_decade2 <- FF5[dates_decade2, 1:5]
  lReg_decade2 <- lReg_decade2_inv <- list()
  lCoeff_decade2 <- lCoeff_decade2_inv <- list()
  
  
  for ( j in 1:ncol(sim_sort) ) {
    
    reg <- regression( Y_train = sim_sort[dates, j], 
                       X_train = X_train, 
                       type = "ols" )
    lReg[[j]] <- reg
    lCoeff[[j]] <- reg$coeffmat
    
    reg_inv <- regression( Y_train = sim_sort_inv[dates, j], 
                           X_train = X_train, 
                           type = "ols" )
    lReg_inv[[j]] <- reg_inv
    lCoeff_inv[[j]] <- reg_inv$coeffmat
    
    reg_decade1 <- regression( Y_train = sim_sort[dates, j], 
                               X_train = X_train_decade1, 
                               type = "ols" )
    lReg_decade1[[j]] <- reg_decade1
    lCoeff_decade1[[j]] <- reg_decade1$coeffmat
    
    reg_decade1_inv <- regression( Y_train = sim_sort_inv[dates, j], 
                                   X_train = X_train_decade1, 
                                   type = "ols" )
    lReg_decade1_inv[[j]] <- reg_decade1_inv
    lCoeff_decade1_inv[[j]] <- reg_decade1_inv$coeffmat
    
    reg_decade2 <- regression( Y_train = sim_sort[dates, j], 
                             X_train = X_train_decade2, 
                             type = "ols" )
    lReg_decade2[[j]] <- reg_decade2
    lCoeff_decade2[[j]] <- reg_decade2$coeffmat
    
    reg_decade2_inv <- regression( Y_train = sim_sort_inv[dates, j], 
                                 X_train = X_train_decade2, 
                                 type = "ols" )
    lReg_decade2_inv[[j]] <- reg_decade2_inv
    lCoeff_decade2_inv[[j]] <- reg_decade2_inv$coeffmat
    
  }  
  
  beta_mat <- do.call( cbind, lapply( lCoeff, FUN = function(x) { x[-1, 1] } ) )
  beta_pval_mat <- do.call( cbind, lapply( lCoeff, FUN = function(x) { x[-1, 4] } ) )
  alpha_score <- unlist( lapply( lCoeff, FUN = function(x) { x[1, 1] } ) )
  p_values <- unlist( lapply( lCoeff, FUN = function(x) { x[1, 4] } ) )
  
  beta_mat_inv <- do.call( cbind, lapply( lCoeff_inv, FUN = function(x) { x[-1, 1] } ) )
  beta_pval_mat_inv <- do.call( cbind, lapply( lCoeff_inv, FUN = function(x) { x[-1, 4] } ) )
  alpha_score_inv <- unlist( lapply( lCoeff_inv, FUN = function(x) { x[1, 1] } ) )
  p_values_inv <- unlist( lapply( lCoeff_inv, FUN = function(x) { x[1, 4] } ) )
  
  beta_mat_decade1 <- do.call( cbind, lapply( lCoeff_decade1, FUN = function(x) { x[-1, 1] } ) )
  beta_pval_mat_decade1 <- do.call( cbind, lapply( lCoeff_decade1, FUN = function(x) { x[-1, 4] } ) )
  alpha_score_decade1 <- unlist( lapply( lCoeff_decade1, FUN = function(x) { x[1, 1] } ) )
  p_values_decade1 <- unlist( lapply( lCoeff_decade1, FUN = function(x) { x[1, 4] } ) )
  
  beta_mat_decade1_inv <- do.call( cbind, lapply( lCoeff_decade1_inv, FUN = function(x) { x[-1, 1] } ) )
  beta_pval_mat_decade1_inv <- do.call( cbind, lapply( lCoeff_decade1_inv, FUN = function(x) { x[-1, 4] } ) )
  alpha_score_decade1_inv <- unlist( lapply( lCoeff_decade1_inv, FUN = function(x) { x[1, 1] } ) )
  p_values_decade1_inv <- unlist( lapply( lCoeff_decade1_inv, FUN = function(x) { x[1, 4] } ) )
  
  beta_mat_decade2 <- do.call( cbind, lapply( lCoeff_decade2, FUN = function(x) { x[-1, 1] } ) )
  beta_pval_mat_decade2 <- do.call( cbind, lapply( lCoeff_decade2, FUN = function(x) { x[-1, 4] } ) )
  alpha_score_decade2 <- unlist( lapply( lCoeff_decade2, FUN = function(x) { x[1, 1] } ) )
  p_values_decade2 <- unlist( lapply( lCoeff_decade2, FUN = function(x) { x[1, 4] } ) )
  
  beta_mat_decade2_inv <- do.call( cbind, lapply( lCoeff_decade2_inv, FUN = function(x) { x[-1, 1] } ) )
  beta_pval_mat_decade2_inv <- do.call( cbind, lapply( lCoeff_decade2_inv, FUN = function(x) { x[-1, 4] } ) )
  alpha_score_decade2_inv <- unlist( lapply( lCoeff_decade2_inv, FUN = function(x) { x[1, 1] } ) )
  p_values_decade2_inv <- unlist( lapply( lCoeff_decade2_inv, FUN = function(x) { x[1, 4] } ) )
  
  
  
  
  
  
  
  
  barplot( beta_mat, beside = TRUE, col = 1:nrow(beta_mat) )
  abline( h = 1, col = "grey" )
  barplot( beta_pval_mat, beside = TRUE, col = 1:nrow(beta_pval_mat) )
  
  barplot( beta_mat_inv, beside = TRUE, col = 1:nrow(beta_mat_inv) )
  barplot( beta_pval_mat_inv, beside = TRUE, col = 1:nrow(beta_pval_mat_inv) )
  
  barplot( beta_mat_decade1, beside = TRUE, col = 1:nrow(beta_mat_decade1) )
  barplot( beta_pval_mat_decade1, beside = TRUE, col = 1:nrow(beta_pval_mat_decade1) )
  
  barplot( beta_mat_decade1_inv, beside = TRUE, col = 1:nrow(beta_mat_decade1_inv) )
  barplot( beta_pval_mat_decade1_inv, beside = TRUE, col = 1:nrow(beta_pval_mat_decade1_inv) )
  
  barplot( beta_mat_decade2, beside = TRUE, col = 1:nrow(beta_mat_decade2) )
  barplot( beta_pval_mat_decade2, beside = TRUE, col = 1:nrow(beta_pval_mat_decade2) )
  
  barplot( beta_mat_decade2_inv, beside = TRUE, col = 1:nrow(beta_mat_decade2_inv) )
  barplot( beta_pval_mat_decade2_inv, beside = TRUE, col = 1:nrow(beta_pval_mat_decade2_inv) )
  
  
  descStats( X_train )
  descStats( window(X_train, start(X_train), start(X_train_decade2)) )
  descStats( X_train_decade2 )
  
  
  plot( as.simTS(X_train) )
  plot( log(cumulated(X_train, "discrete")), plot.type = "single" )
  abline(h = 0, col = "grey")
  plot( cumulated(X_train_decade2, "discrete"), plot.type = "single" )
  
  
    
  
  
  barplot( p_values )
  abline( h = 0.05 )
  length( which( p_values <= 0.05 ) ) / length(p_values)
  
  
  j <- 26
  summary( lReg[[j]]$reg )
  
  
  alpha_score
  barplot( alpha_score )
  abline( h = 0.05 )

  
  
  t(beta_mat)
  t(beta_pval_mat)  
  
  
 
  
 
  
  
    
  
    
  
  # --------------------------------------------------------------------------
  # ACCEPTANCE-REJECTION PORTFOLIOS
  # --------------------------------------------------------------------------
  
  BT4 <- BT$copy()
  BT4$spec$portfolio <- "accrej"
  BT4$spec$th <- llPortf_scores_2[[1]] / 2
  BT4$data$Score <- lScore[["A"]]  
  # debugonce( BT4$accrejPortfolio )
  BT4$runLoop()
  
  
  headleft( BT4$output$acc_idx )
  
  
  

  plot( llPortf_scores_2[[26]] )
  plot( 1 / llPortf_scores_2[[1]] )
  

  
  apply( lScore[["Z"]], 1, sum )
  
  
  
  
  
  
  
  
  
  
  
  
      
  