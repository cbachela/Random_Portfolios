
  
  ############################################################################
  ### RANDOM PORTFOLIOS - GEOMETRIC RANDOM WALK SAMPLING - MSCI MULFIACTOR CONSTRAINTS - USA
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     15.01.2023
  # First version:    07.01.2023
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
  # CONSTRAINTS - MSCI MULTIFACTOR INDEX METHODLOLOGY 
  # --------------------------------------------------------------------------
  
  BT_MF <- BT0$copy()
  # BT_MF$spec$selection_filter <- c("db_flag", "upperbound", "sector", "country", "NA")
  BT_MF$spec$selection_filter <- c("upperbound", "sector", "country", "NA")
  BT_MF$spec$width <- 52
  BT_MF$spec$cons_fun <- c("box", "sector", "sectorLB", "vola") # Notice that we approximate the quadratic variance constraint with a linear version.
  BT_MF$spec$vola_multiple <- 1.1
  BT_MF$data <- BT0$data
  BT_MF$MSCIMultifactorConstraints()
  
  
  
  
  
  
  headleft( BT_MF$data$country_UB )
  headleft( BT_MF$data$country_LB )
  headleft( BT_MF$data$sector_UB )
  headleft( BT_MF$data$sector_LB )
  headleft( BT_MF$data$upper_mat )
  headleft( BT_MF$data$lower_mat )
  
  apply( BT_MF$data$upper_mat, 1, function(x) { sum( na.omit(x) ) } )
  apply( BT_MF$data$lower_mat, 1, function(x) { sum( na.omit(x) ) } )
  apply( BT_MF$data$sector_UB, 1, function(x) { sum( na.omit(x) ) } )
  apply( BT_MF$data$sector_LB, 1, function(x) { sum( na.omit(x) ) } )
  
  plot( apply( BT_MF$data$upper_mat, 1, function(x) { sum( na.omit(x) ) } ) )
  plot( apply( BT_MF$data$lower_mat, 1, function(x) { sum( na.omit(x) ) } ) )
  plot( apply( BT_MF$data$sector_UB, 1, function(x) { sum( na.omit(x) ) } ) )
  plot( apply( BT_MF$data$sector_LB, 1, function(x) { sum( na.omit(x) ) } ) )
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BACKTEST QUASI CAPW, I.E., BM WITHIN CONSTRAINTS
  # --------------------------------------------------------------------------
  
  BacktestCustom.qcapwPortfolio <- function( rebdate = NULL )
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
  BacktestCustom$methods( qcapwPortfolio = BacktestCustom.qcapwPortfolio )
  
  
  BT_qcapw <- BT_MF$copy()
  BT_qcapw$spec$portfolio <- "qcapw"
  BT_qcapw$spec$GPS@covariance$method <- "duv"
  BT_qcapw$runLoop()
  sim_qcapw <- BT_qcapw$simulate( fc = 0, vc = 0 )
  sim_qcapw <- sim_qcapw[isWeekday(time(sim_qcapw)), ]
  
  
  
  
  # --------------------------------------------------------------------------
  # BACKTEST QUASI EQW PORTFOLIO
  # --------------------------------------------------------------------------
  
  BT_qeqw <- BT_MF$copy()
  BT_qeqw$spec$portfolio <- NULL
  BT_qeqw$spec$GPS@covariance$method <- "duv"
  BT_qeqw$runLoop()
  sim_qeqw <- BT_qeqw$simulate( fc = 0, vc = 0 )
  sim_qeqw <- sim_qeqw[isWeekday(time(sim_qeqw)), ]
  
  
  # Add qeqw and qcapw to X_bm
  X_bm <- na.omit( cbind( X_fact, QEQW = sim_qeqw, QCapw = sim_qcapw ) )
  
  descStats( X_bm )
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BACKTEST UNIFORM RP'S WITHIN CONSTRAINTS
  # --------------------------------------------------------------------------
  
  # BT1 <- BT_MF$copy()
  # BT1$spec$keep_all_rebdates <- FALSE
  # BT1$spec$portfolio <- "grw"
  # BT1$spec$name <- "grw_multifac_usa"
  # BT1$spec$n_sim <- 10^2 * 5
  # BT1$output <- list()
  # # debugonce( BT1$grwPortfolio )
  # BT1$runLoop()  
  # 
  # # Simulate
  # tmp <- lapply( BT1$output, FUN = simFUN )
  # sim_1 <- do.call( cbind, tmp )
  # sim_1 <- sim_1[isWeekday(time(sim_1)), ]
  # colnames(sim_1) <- paste0("rp_grw_uniform", 1:ncol(sim_1))
  # BT1$data$simulations <- sim_1
  # 
  # # Apply FEV-bias
  # lWeights <- BT1$output
  # lWeights_fev <- lapply( lWeights, FUN = function(wmat) { fevBias( x = wmat, q = 6 ) } )
  # 
  # # Project FEV-biased weights to boundary of feasible set
  # lWeights_cons <- lWeights
  # for ( i in 1:length(lWeights_cons) ) {
  #   for ( today in rownames(lWeights_cons[[i]]) ) {
  #     RBE <- rebalenv( rebdate = today )
  #     # Cons <- getConstraints(RBE$GPS)
  #     lincon <- getConstraints( RBE$GPS, "bounds" )
  #     upper_bounds <- lincon$upper
  #     w_unc <- lWeights_fev[[i]][today, ]
  #     # idx_notNA <- which(!is.na(w_unc))
  #     idx_notNA <- which( !is.na( lWeights[[i]][today, ] ) )
  #     w_unc <- w_unc[idx_notNA]
  #     w_unc <- setNames( as.numeric(w_unc), names(w_unc))
  #     # debugonce( map2boxcon )
  #     lWeights_cons[[i]][today, idx_notNA] <- map2boxcon( w_unc = w_unc, 
  #                                                         upper = upper_bounds, 
  #                                                         itermax = 10^3 )
  #   }
  # }
  # BT1_fev2cons <- BT1$copy()
  # BT1_fev2cons$spec$name <- paste0(BT1$spec$name, "_fev2cons")
  # BT1_fev2cons$output <- lWeights_cons
  # 
  # # Simulate
  # tmp <- lapply( lWeights_cons, FUN = simFUN )
  # sim_1_cons <- do.call( cbind, tmp )
  # sim_1_cons <- sim_1_cons[isWeekday(time(sim_1_cons)), ]
  # colnames(sim_1_cons) <- paste0("rp_grw_uniform_fev2cons", 1:ncol(sim_1_cons))
  # BT1_fev2cons$data$simulations <- sim_1_cons
  # 
  # 
  # # # Save
  # # BT1$save( wd = paste0(wd_data, "waRehouse/"), name = BT1$spec$name, without_data = FALSE )
  # # BT1_fev2cons$save( wd = paste0(wd_data, "waRehouse/"), name = BT1_fev2cons$spec$name, without_data = FALSE )
  
  
  # Load
  BT1 <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "grw_multifac_usa" )
  BT1_fev2cons <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "grw_multifac_usa_fev2cons" )
  sim_1 <- BT1$data$simulations
  sim_1_cons <- BT1_fev2cons$data$simulations
  
  
  
  
  
  
  
  
  
  # Combine
  end_date <- "2022-03-06"
  strategy_names <- colnames(X_bm)
  sim_all <- na.omit( cbind( X_bm[ ,strategy_names], sim_1_cons) )
  sim_all <- window( sim_all, start(sim_all), end_date )
  lStats_all <-  descStats( sim_all, descStatsSpec(annualize = TRUE) )
  lStats_sim1 <- descStats( window( sim_1, start(sim_1), end_date), descStatsSpec(annualize = TRUE) )
  
  

  # Risk return plot
  
  lStats <- lStats_all
  df <- data.frame( x = lStats$stats["sds", ],
                    y = lStats$stats["means", ] )
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
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # SORT RP'S BY FACTOR SCORE, CHECK RANK VS. PERFORMANCE
  # --------------------------------------------------------------------------
  
  # debugonce( list2array )
  w_array <- list2array( L = BT1$output )
  sim <- sim_1
  
  w_array <- list2array( L = BT1_fev2cons$output )
  sim <- sim_1_cons
  
 
  
  # Random concatenation
  lStats_rnd <- descStats( sim )
  
  
  
  # Compute Momentum scores
  BT_mom <- BT_MF$copy()
  BT_mom$spec$width <- 52
  BT_mom$spec$portfolio <- "momScore"
  BT_mom$runLoop()
  
  Score_mom <- BT_mom$output$scores
  sim_mom_is <- simSortedBy( sim = sim,
                             BT0 = BT_MF,
                             Score = Score_mom,
                             w_array = w_array,
                             insample = TRUE )
  sim_mom_is_score <- attr(sim_mom_is, "score")
  sim_mom_oos <- simSortedBy( sim = sim,
                              BT0 = BT_MF,
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
  BT_maxmom <- BT_MF$copy()
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
                              BT0 = BT_MF,
                              Score = Score_size,
                              w_array = w_array,
                              insample = TRUE )
  sim_size_is_score <- attr(sim_size_is, "score")
  sim_size_oos <- simSortedBy( sim = sim,
                               BT0 = BT_MF,
                               Score = Score_size,
                               w_array = w_array,
                               insample = FALSE )
  sim_size_oos_score <- attr(sim_size_oos, "score")
  sim_size_is_score <- attr(sim_size_is, "score")
  lStats_size_is <- descStats( sim_size_is )
  stats_size_is <- lStats_size_is$stats
  lStats_size_oos <- descStats( sim_size_oos )
  stats_size_oos <- lStats_size_oos$stats
  
  # log(Size) score
  Score_size_log <- log(BT0$data$wmat_bm)
  sim_size_log_is <- simSortedBy( sim = sim,
                              BT0 = BT_MF,
                              Score = Score_size_log,
                              w_array = w_array,
                              insample = TRUE )
  sim_size_is_log_score <- attr(sim_size_is, "score")
  sim_size_log_oos <- simSortedBy( sim = sim,
                               BT0 = BT_MF,
                               Score = Score_size_log,
                               w_array = w_array,
                               insample = FALSE )
  sim_size_log_oos_score <- attr(sim_size_log_oos, "score")
  sim_size__log_is_score <- attr(sim_size_log_is, "score")
  lStats_size_log_is <- descStats( sim_size_log_is )
  stats_size_log_is <- lStats_size_log_is$stats
  lStats_size_log_oos <- descStats( sim_size_log_oos )
  stats_size_log_oos <- lStats_size_log_oos$stats
  
  # Size stdz
  standardizeScores <- function(scores)
  {
    idx_notNA <- which( !is.na(scores) )
    x <- scores[idx_notNA]
    # Z-score
    x_stdz <- ( x - mean(x)) / sd(x)
    # Winsorize
    x_stdz[ which(x_stdz < -3) ] <- -3
    x_stdz[ which(x_stdz > 3) ] <- 3
    scores_stdz <- scores
    scores_stdz[ idx_notNA ] <- x_stdz
    return( scores_stdz )
  }
  Score_size_stdz <- timeSeries( t( apply( Score_size_log, 1, FUN = standardizeScores ) ),
                                 time(Score_size) )
  sim_size_stdz_is <- simSortedBy( sim = sim,
                                  BT0 = BT_MF,
                                  Score = Score_size_stdz,
                                  w_array = w_array,
                                  insample = TRUE )
  sim_size_is_stdz_score <- attr(sim_size_is, "score")
  sim_size_stdz_oos <- simSortedBy( sim = sim,
                                   BT0 = BT_MF,
                                   Score = Score_size_stdz,
                                   w_array = w_array,
                                   insample = FALSE )
  sim_size_stdz_oos_score <- attr(sim_size_stdz_oos, "score")
  sim_size__stdz_is_score <- attr(sim_size_stdz_is, "score")
  lStats_size_stdz_is <- descStats( sim_size_stdz_is )
  stats_size_stdz_is <- lStats_size_stdz_is$stats
  lStats_size_stdz_oos <- descStats( sim_size_stdz_oos )
  stats_size_stdz_oos <- lStats_size_stdz_oos$stats
 
  
  
  
  
  
  # Compute Vola scores
  BT_vol <- BT_MF$copy()
  BT_vol$spec$width <- 52
  BT_vol$spec$portfolio <- "volScore"
  BT_vol$runLoop()
  
 
  Score_vol <- BT_vol$output$scores
  sim_vol_is <- simSortedBy( sim = sim,
                             BT0 = BT_MF,
                             Score = Score_vol,
                             w_array = w_array,
                             insample = TRUE )
  sim_vol_is_score <- attr(sim_vol_is, "score")
  sim_vol_oos <- simSortedBy( sim = sim,
                              BT0 = BT_MF,
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
                              BT0 = BT_MF,
                              Score = Score_qual,
                              w_array = w_array,
                              insample = TRUE )
  sim_qual_is_score <- attr(sim_qual_is, "score")
  sim_qual_oos <- simSortedBy( sim = sim,
                               BT0 = BT_MF,
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
                             BT0 = BT_MF,
                             Score = Score_val,
                             w_array = w_array,
                             insample = TRUE )
  sim_val_is_score <- attr(sim_val_is, "score")
  sim_val_oos <- simSortedBy( sim = sim,
                              BT0 = BT_MF,
                              Score = Score_val,
                              w_array = w_array,
                              insample = FALSE )
  sim_val_oos_score <- attr(sim_val_oos, "score")
  lStats_val_is <- descStats( sim_val_is )
  stats_val_is <- lStats_val_is$stats
  lStats_val_oos <- descStats( sim_val_oos )
  stats_val_oos <- lStats_val_oos$stats
  
  
  lScore <- list( momentum = Score_mom,
                  value = Score_val,
                  quality = Score_qual,
                  size = Score_size,
                  size_log = Score_size_log,
                  size_stdz = Score_size_stdz,
                  lowvola = Score_vol )
  
  lSim_is <- list(  momentum = sim_mom_is,
                    value = sim_val_is,
                    quality = sim_qual_is,
                    size = sim_size_is,
                    size_log = sim_size_log_is,
                    size_stdz = sim_size_stdz_is,
                    lowvola = sim_vol_is,
                    rnd_sort = sim )
  lSim_oos <- list(  momentum = sim_mom_oos,
                     value = sim_val_oos,
                     quality = sim_qual_oos,
                     size = sim_size_oos,
                     size_log = sim_size_log_oos,
                     size_stdz = sim_size_stdz_oos,
                     lowvola = sim_vol_oos,
                     rnd_sort = sim )
  lStats_is <- list( momentum = lStats_mom_is,
                     value = lStats_val_is,
                     quality = lStats_qual_is,
                     size = lStats_size_is,
                     size_log = lStats_size_log_is,
                     size_stdz = lStats_size_stdz_is,
                     lowvola = lStats_vol_is,
                     rnd_sort = lStats_rnd )
  lStats_oos <- list( momentum = lStats_mom_oos,
                      value = lStats_val_oos,
                      quality = lStats_qual_oos,
                      size = lStats_size_oos,
                      size_log = lStats_size_log_oos,
                      size_stdz = lStats_size_stdz_oos,
                      lowvola = lStats_vol_oos,
                      rnd_sort = lStats_rnd )
  
  
  factor_names <- names(lStats_oos)
  colors <- fBasics::divPalette(n = ncol(sim), "RdYlGn")
  
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
  
  
  # Correlations btw. variance and means
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
  
  
  boxplot( do.call( cbind, lapply( lStats_oos, FUN = function(x) { x$stats["means", ]} ) ) )
  boxplot( do.call( cbind, lapply( lStats_oos, FUN = function(x) { x$stats["sds", ]} ) ) )
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # RISK - RETURN PLOT WITH SORTED SIMULATIONS
  # --------------------------------------------------------------------------
  
  
  # Combine
  end_date <- "2022-03-06"
  strategy_names <- colnames(X_bm)
  sim_all <- na.omit( cbind( X_bm[ ,strategy_names], sim_vol_oos) )   # ~~~~~~~~~~~~~ 
  # sim_all <- na.omit( cbind( X_bm[ ,strategy_names], sim_mom_oos) )   # ~~~~~~~~~~~~~ 
  sim_all <- window( sim_all, start(sim_all), end_date )
  lStats_all <- descStats( sim_all, descStatsSpec(annualize = TRUE) )
  lStats_sim1 <- descStats( window( sim_1, start(sim_1), end_date), descStatsSpec(annualize = TRUE) )
  
  
  ###
  # Risk return plot
  
  lStats <- lStats_all
  df <- data.frame( x = lStats$stats["sds", ],
                    y = lStats$stats["means", ] )
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
  
  
  
  # --------------------------------------------------------------------------
  # BEAR BULL ANALYSIS
  # --------------------------------------------------------------------------
  
  
  BBSObj <- BBSRC$new()
  BBSObj$setCtrl()
  BBSObj$data$X_level <- cumulated(X_bm[ ,"Capw"], "discrete")
  
  bm_names <- c("Capw", "QCapw", "Momentum", "Multifactor", "Size", "Enhanced_Value", "QEQW")
  factor_names <- names(lSim_oos)
  lPhase_stats <- list()
  for ( factor_name in factor_names ) {
    BBSObj$data$X <- na.omit( cbind( X_bm[ ,bm_names], 
                                     lSim_oos[[ factor_name ]] ) )
    BBSObj$runRobust()
    lPhase_stats[[ factor_name ]] <- BBSObj$phaseStats()
  }
  
  
  
  par( mfrow = c(2, 1) )
  colors <- fBasics::divPalette(n = ncol(sim), "RdYlGn")
  factor_name <- "momentum"
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
  colors <- 1:length(bm_names)
  legend( "topleft", bm_names, lwd = 2, col = colors, text.col = colors, bty = "n")
  
  
  
  
  stats_field <- "sharpe"
  ldens_bull <- lapply( lPhase_stats, FUN = function(x) { density(x$`states==1`[stats_field, ]) })
  ldens_bear <- lapply( lPhase_stats, FUN = function(x) { density(x$`states==-1`[stats_field, ]) })
  par( mfrow = c(2, 1) )
  slolz:::plot.ldensity( ldens_bull, fillin = TRUE, main = "Bull" )
  slolz:::plot.ldensity( ldens_bear, fillin = TRUE, main = "Bear" )
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # FACTOR SCORES BM
  # --------------------------------------------------------------------------
  
  lScore_bm <- list()
  for ( factor_name in names(lScore) ) {
    lScore_bm[[ factor_name ]] <- portfolioScore( wmat = BT_MF$data$wmat_bm, 
                                                  score_mat = lScore[[ factor_name ]] )
  }
  factor_scores_bm <- do.call( cbind, lScore_bm )
  
  plot( factor_scores_bm )
  lapply( lScore, FUN = function(x) { range( na.omit( as.numeric( x ) ) ) } )
  
  lFactor <- lapply( lScore, FUN = function(x) { na.omit(x, method = "z" ) } )
  
  
  

  BT_Fctrl <- BacktestPanelFactor$new()
  BT_Fctrl$spec <- BT_MF$spec
  BT_Fctrl$data <- BT_MF$data
  BT_Fctrl$data$factor_scores_bm <- factor_scores_bm
  BT_Fctrl$data$lFactor <- lFactor
  BT_Fctrl$spec$factor_names_constraints <- "value"
  BT_Fctrl$spec$cons_fun <- c(BT_MF$spec$cons_fun, "controlbmex10")
  # BT_Fctrl$spec$cons_fun <- c(BT_MF$spec$cons_fun, "controlbmex00")
  # BT_Fctrl$spec$cons_fun <- c(BT_MF$spec$cons_fun, "controlbmex25")

  

  # QEQW
  BT_Fctrl$spec$portfolio <- NULL
  BT_Fctrl$spec$GPS@covariance$method <- "duv"
  BT_Fctrl$runLoop()
  sim_fctrl_qeqw <- BT_Fctrl$simulate( fc = 0, vc = 0 )
  sim_fctrl_qeqw <- sim_fctrl_qeqw[isWeekday(time(sim_fctrl_qeqw)), ]
  colnames(sim_fctrl_qeqw) <- "qeqw_fctrl_val"
  
  
  X_tmp <- na.omit( cbind( X_bm[ ,c("Value", "QCapw", "QEQW")], sim_fctrl_qeqw ) )
  descStats( X_tmp ) 
  
  plotSimTS( as.simTS( X_tmp ) )   
  
  
  # Check
  tmp <- portfolioScore( BT_Fctrl$output$weights, score_mat = lFactor[[ "value" ]] )
  check <- na.omit( cbind( tmp, factor_scores_bm[ ,"value"] ) )
  cbind( check, check[ ,2] - 0.1, check[ ,2] + 0.1)
  

  
  
    
  
  
  # --------------------------------------------------------------------------
  # BACKTEST UNIFORM RP'S WITHIN CONSTRAINTS WITH VALUE CONTROL
  # --------------------------------------------------------------------------

  BT2 <- BT_Fctrl$copy()
  BT2$spec$factor_names <- "value"
  BT2$spec$factor_names_constraints <- BT2$spec$factor_names
  BT2$spec$keep_all_rebdates <- FALSE
  BT2$spec$portfolio <- "grw"
  BT2$spec$n_sim <- 10^2 * 5
  BT2$spec$name <- "grw_multifac_usa_valctrl"
  BT2$output <- list()
  BT2$runLoop()

  # Simulate
  tmp <- lapply( BT2$output[ grepl("grw", names(BT2$output)) ], FUN = simFUN )
  sim_2 <- do.call( cbind, tmp )
  sim_2 <- sim_2[isWeekday(time(sim_2)), ]
  colnames(sim_2) <- paste0("rp_grw_uniform_valctrl", 1:ncol(sim_2))
  BT2$data$simulations <- sim_2

  # Apply FEV-bias
  lWeights_ctrl <- BT2$output[ grepl("grw", names(BT2$output)) ]
  lWeights_fev_ctrl <- lapply( lWeights_ctrl, FUN = function(wmat) { fevBias( x = wmat, q = 6 ) } )

  # Project FEV-biased weights to boundary of feasible set
  lWeights_cons_ctrl <- lWeights_ctrl
  for ( i in 1:length(lWeights_cons_ctrl) ) {
    for ( today in rownames(lWeights_cons_ctrl[[i]]) ) {
      RBE <- rebalenv( rebdate = today )
      # Cons <- getConstraints(RBE$GPS)
      lincon <- getConstraints( RBE$GPS, "bounds" )
      upper_bounds <- lincon$upper
      w_unc <- lWeights_fev_ctrl[[i]][today, ]
      # idx_notNA <- which(!is.na(w_unc))
      idx_notNA <- which( !is.na( lWeights_ctrl[[i]][today, ] ) )
      w_unc <- w_unc[idx_notNA]
      w_unc <- setNames( as.numeric(w_unc), names(w_unc))
      # debugonce( map2boxcon )
      lWeights_cons_ctrl[[i]][today, idx_notNA] <- map2boxcon( w_unc = w_unc,
                                                               upper = upper_bounds,
                                                               itermax = 10^3 )
    }
  }
  BT2_fev2cons <- BT2$copy()
  BT2_fev2cons$spec$name <- paste0(BT2$spec$name, "_fev2cons")
  BT2_fev2cons$output <- lWeights_cons_ctrl

  # Simulate
  tmp <- lapply( BT2_fev2cons$output, FUN = simFUN )
  sim_2_cons <- do.call( cbind, tmp )
  sim_2_cons <- sim_2_cons[isWeekday(time(sim_2_cons)), ]
  colnames(sim_2_cons) <- paste0("rp_grw_uniform_cons_valctrl", 1:ncol(sim_2_cons))
  BT2_fev2cons$data$simulations <- sim_2_cons

  # # Save
  # BT2$save( wd = paste0(wd_data, "waRehouse/"), name = BT2$spec$name, without_data = FALSE )
  # BT2_fev2cons$save( wd = paste0(wd_data, "waRehouse/"), name = BT2_fev2cons$spec$name, without_data = FALSE )
  
  
  # Load
  BT2 <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "grw_multifac_usa_valctrl" )
  BT2_fev2cons <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "grw_multifac_usa_valctrl_fev2cons" )
  sim_2 <- BT2$data$simulations
  sim_2_cons <- BT2_fev2cons$data$simulations
  
  
  
  
  
  
  

  
  
  
  
  
  
  
  
  # Check
  tmp <- portfolioScore(BT2$output[[2]], score_mat = lFactor[[ "value" ]] )
  check <- na.omit( cbind( tmp, factor_scores_bm[ ,"value"] ) )
  cbind( check, check[ ,2] - 0.1, check[ ,2] + 0.1)
  
  
  
  
  
  # Combine
  end_date <- "2022-03-06"
  strategy_names <- colnames(X_bm)
  # sim_all <- na.omit( cbind( X_bm[ ,strategy_names], sim_1_cons) )
  sim_all <- na.omit( cbind( X_bm[ ,strategy_names], sim_2_cons) )
  sim_all <- window( sim_all, start(sim_all), end_date )
  lStats_all <-  descStats( sim_all, descStatsSpec(annualize = TRUE) )
  
  lStats_sim1 <- descStats( window( sim_1, start(sim_1), end_date), descStatsSpec(annualize = TRUE) )
  lStats_sim1_cons <- descStats( window( sim_1_cons, start(sim_1_cons), end_date), descStatsSpec(annualize = TRUE) )
  
  lStats_sim2 <- descStats( window( sim_2, start(sim_2), end_date), descStatsSpec(annualize = TRUE) )
  lStats_sim2_cons <- descStats( window( sim_2_cons, start(sim_2_cons), end_date), descStatsSpec(annualize = TRUE) )
  
  plot( x = lStats_sim1$stats["sds", ], y = lStats_sim1$stats["means", ], pch = 19 )
  points( x = lStats_sim2$stats["sds", ], y = lStats_sim2$stats["means", ], pch = 19, col = 2 )
  
  plot( x = lStats_sim1_cons$stats["sds", ], y = lStats_sim1_cons$stats["means", ], pch = 19 )
  points( x = lStats_sim2_cons$stats["sds", ], y = lStats_sim2_cons$stats["means", ], pch = 19, col = 2 )
  
  
  
  
  
  
  ###
  # Risk return plot
  
  lStats <- lStats_all
  df <- data.frame( x = lStats$stats["sds", ],
                    y = lStats$stats["means", ] )
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

  
  
  

  
  # Check
  tmp <- portfolioScore(BT2$output[[2]], score_mat = lFactor[[ "value" ]] )
  check <- na.omit( cbind( tmp, factor_scores_bm[ ,"value"] ) )
  cbind( check, check[ ,2] - 0.1, check[ ,2] + 0.1)
  
  headleft( BT2$output$value )
  headleft( BT2$output$grw1 )
  apply( BT2$output$grw1, 1, function(x) { sum(na.omit(x)) } )  
  apply( BT1$output$grw1, 1, function(x) { sum(na.omit(x)) } )  
  apply( lWeights_cons[[1]], 1, function(x) { sum(na.omit(x)) } )  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BACKTEST UNIFORM RP'S WITHIN CONSTRAINTS WITH QUALITY CONTROL
  # --------------------------------------------------------------------------
   
  BT3 <- BT_Fctrl$copy()
  BT3$spec$factor_names <- "quality"
  BT3$spec$factor_names_constraints <- BT4$spec$factor_names
  BT3$spec$keep_all_rebdates <- FALSE
  BT3$spec$portfolio <- "grw"
  BT3$spec$n_sim <- 10^2 * 5
  BT3$spec$name <- "grw_multifac_usa_qualctrl"
  BT3$output <- list()
  BT3$runLoop()

  # Simulate
  tmp <- lapply( BT3$output[ grepl("grw", names(BT3$output)) ], FUN = simFUN )
  sim_3 <- do.call( cbind, tmp )
  sim_3 <- sim_3[isWeekday(time(sim_3)), ]
  colnames(sim_3) <- paste0("rp_grw_uniform_qualctrl", 1:ncol(sim_3))
  BT3$data$simulations <- sim_3

  # Apply FEV-bias
  lWeights_ctrl <- BT3$output[ grepl("grw", names(BT3$output)) ]
  lWeights_fev_ctrl <- lapply( lWeights_ctrl, FUN = function(wmat) { fevBias( x = wmat, q = 6 ) } )

  # Project FEV-biased weights to boundary of feasible set
  lWeights_cons_ctrl <- lWeights_ctrl
  for ( i in 1:length(lWeights_cons_ctrl) ) {
    for ( today in rownames(lWeights_cons_ctrl[[i]]) ) {
      RBE <- rebalenv( rebdate = today )
      # Cons <- getConstraints(RBE$GPS)
      lincon <- getConstraints( RBE$GPS, "bounds" )
      upper_bounds <- lincon$upper
      w_unc <- lWeights_fev_ctrl[[i]][today, ]
      # idx_notNA <- which(!is.na(w_unc))
      idx_notNA <- which( !is.na( lWeights_ctrl[[i]][today, ] ) )
      w_unc <- w_unc[idx_notNA]
      w_unc <- setNames( as.numeric(w_unc), names(w_unc))
      # debugonce( map2boxcon )
      lWeights_cons_ctrl[[i]][today, idx_notNA] <- map2boxcon( w_unc = w_unc,
                                                               upper = upper_bounds,
                                                               itermax = 10^3 )
    }
  }
  BT3_fev2cons <- BT3$copy()
  BT3_fev2cons$spec$name <- paste0(BT3$spec$name, "_fev2cons")
  BT3_fev2cons$output <- lWeights_cons_ctrl

  # Simulate
  tmp <- lapply( BT3_fev2cons$output, FUN = simFUN )
  sim_3_cons <- do.call( cbind, tmp )
  sim_3_cons <- sim_3_cons[isWeekday(time(sim_3_cons)), ]
  colnames(sim_3_cons) <- paste0("rp_grw_uniform_cons_qualctrl", 1:ncol(sim_3_cons))
  BT3_fev2cons$data$simulations <- sim_3_cons

  # # Save
  # BT3$save( wd = paste0(wd_data, "waRehouse/"), name = BT3$spec$name, without_data = FALSE )
  # BT3_fev2cons$save( wd = paste0(wd_data, "waRehouse/"), name = BT3_fev2cons$spec$name, without_data = FALSE )
  
  
  # Load
  BT3 <- loadBacktest( wd =paste0(wd_data, "waRehouse/"), name = "grw_multifac_usa_qualctrl" )
  BT3_fev2cons <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "grw_multifac_usa_qualctrl_fev2cons" )
  sim_3 <- BT3$data$simulations
  sim_3_cons <- BT3_fev2cons$data$simulations
  
  
  lStats_sim3 <- descStats( window( sim_3, start(sim_3), end_date), descStatsSpec(annualize = TRUE) )
  lStats_sim3_cons <- descStats( window( sim_3_cons, start(sim_3_cons), end_date), descStatsSpec(annualize = TRUE) )
  
  
  plot( x = lStats_sim1$stats["sds", ], y = lStats_sim1$stats["means", ], pch = 19 )
  points( x = lStats_sim2$stats["sds", ], y = lStats_sim2$stats["means", ], pch = 19, col = 2 )
  points( x = lStats_sim3$stats["sds", ], y = lStats_sim3$stats["means", ], pch = 19, col = 3 )
  
  plot( x = lStats_sim1_cons$stats["sds", ], y = lStats_sim1_cons$stats["means", ], pch = 19 )
  points( x = lStats_sim2_cons$stats["sds", ], y = lStats_sim2_cons$stats["means", ], pch = 19, col = 2 )
  points( x = lStats_sim3_cons$stats["sds", ], y = lStats_sim3_cons$stats["means", ], pch = 19, col = 3 )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BACKTEST UNIFORM RP'S WITHIN CONSTRAINTS WITH SIZE CONTROL
  # --------------------------------------------------------------------------
  
  BT4 <- BT_Fctrl$copy()
  BT4$spec$factor_names <- "size"
  BT4$spec$factor_names_constraints <- BT4$spec$factor_names
  BT4$spec$keep_all_rebdates <- FALSE
  BT4$spec$portfolio <- "grw"
  BT4$spec$n_sim <- 10^2 * 5
  BT4$spec$name <- "grw_multifac_usa_sizectrl"
  BT4$output <- list()
  BT4$runLoop()
  
  # Simulate
  tmp <- lapply( BT4$output[ grepl("grw", names(BT4$output)) ], FUN = simFUN )
  sim_4 <- do.call( cbind, tmp )
  sim_4 <- sim_4[isWeekday(time(sim_4)), ]
  colnames(sim_4) <- paste0("rp_grw_uniform_sizectrl", 1:ncol(sim_4))
  BT4$data$simulations <- sim_4
  
  # Apply FEV-bias
  lWeights_ctrl <- BT4$output[ grepl("grw", names(BT4$output)) ]
  lWeights_fev_ctrl <- lapply( lWeights_ctrl, FUN = function(wmat) { fevBias( x = wmat, q = 6 ) } )
  
  # Project FEV-biased weights to boundary of feasible set
  lWeights_cons_ctrl <- lWeights_ctrl
  for ( i in 1:length(lWeights_cons_ctrl) ) {
    for ( today in rownames(lWeights_cons_ctrl[[i]]) ) {
      RBE <- rebalenv( rebdate = today )
      # Cons <- getConstraints(RBE$GPS)
      lincon <- getConstraints( RBE$GPS, "bounds" )
      upper_bounds <- lincon$upper
      w_unc <- lWeights_fev_ctrl[[i]][today, ]
      # idx_notNA <- which(!is.na(w_unc))
      idx_notNA <- which( !is.na( lWeights_ctrl[[i]][today, ] ) )
      w_unc <- w_unc[idx_notNA]
      w_unc <- setNames( as.numeric(w_unc), names(w_unc))
      # debugonce( map2boxcon )
      lWeights_cons_ctrl[[i]][today, idx_notNA] <- map2boxcon( w_unc = w_unc,
                                                               upper = upper_bounds,
                                                               itermax = 10^3 )
    }
  }
  BT4_fev2cons <- BT4$copy()
  BT4_fev2cons$spec$name <- paste0(BT4$spec$name, "_fev2cons")
  BT4_fev2cons$output <- lWeights_cons_ctrl
  
  # Simulate
  tmp <- lapply( BT4_fev2cons$output, FUN = simFUN )
  sim_4_cons <- do.call( cbind, tmp )
  sim_4_cons <- sim_4_cons[isWeekday(time(sim_4_cons)), ]
  colnames(sim_4_cons) <- paste0("rp_grw_uniform_cons_sizectrl", 1:ncol(sim_4_cons))
  BT4_fev2cons$data$simulations <- sim_4_cons
  
  # # Save
  # BT4$save( wd = paste0(wd_data, "waRehouse/"), name = BT4$spec$name, without_data = FALSE )
  # BT4_fev2cons$save( wd = paste0(wd_data, "waRehouse/"), name = BT4_fev2cons$spec$name, without_data = FALSE )
  
  
  # Load
  BT4 <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "grw_multifac_usa_sizectrl" )
  BT4_fev2cons <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "grw_multifac_usa_sizectrl_fev2cons" )
  sim_4 <- BT4$data$simulations
  sim_4_cons <- BT4_fev2cons$data$simulations
  
  
  lStats_sim4 <- descStats( window( sim_4, start(sim_4), end_date), descStatsSpec(annualize = TRUE) )
  lStats_sim4_cons <- descStats( window( sim_4_cons, start(sim_4_cons), end_date), descStatsSpec(annualize = TRUE) )
  
  
  plot( x = lStats_sim1$stats["sds", ], y = lStats_sim1$stats["means", ], pch = 19 )
  points( x = lStats_sim2$stats["sds", ], y = lStats_sim2$stats["means", ], pch = 19, col = 2 )
  points( x = lStats_sim3$stats["sds", ], y = lStats_sim3$stats["means", ], pch = 19, col = 3 )
  points( x = lStats_sim4$stats["sds", ], y = lStats_sim4$stats["means", ], pch = 19, col = 4 )
 
  plot( x = lStats_sim1_cons$stats["sds", ], y = lStats_sim1_cons$stats["means", ], pch = 19 )
  points( x = lStats_sim2_cons$stats["sds", ], y = lStats_sim2_cons$stats["means", ], pch = 19, col = 2 )
  points( x = lStats_sim3_cons$stats["sds", ], y = lStats_sim3_cons$stats["means", ], pch = 19, col = 3 )
  points( x = lStats_sim4_cons$stats["sds", ], y = lStats_sim4_cons$stats["means", ], pch = 19, col = 4 )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BACKTEST UNIFORM RP'S WITHIN CONSTRAINTS WITH SIZE STANDARDIZED CONTROL
  # --------------------------------------------------------------------------
  
  BT5 <- BT_Fctrl$copy()
  BT5$spec$factor_names_constraints <- "size_stdz"
  BT5$spec$keep_all_rebdates <- FALSE
  BT5$spec$portfolio <- "grw"
  BT5$spec$n_sim <- 10^2 * 5
  BT5$spec$name <- "grw_multifac_usa_sizestdzctrl"
  BT5$output <- list()
  BT5$runLoop()
  
  # Simulate
  tmp <- lapply( BT5$output[ grepl("grw", names(BT5$output)) ], FUN = simFUN )
  sim_5 <- do.call( cbind, tmp )
  sim_5 <- sim_5[isWeekday(time(sim_5)), ]
  colnames(sim_5) <- paste0("rp_grw_uniform_sizestdzctrl", 1:ncol(sim_5))
  BT5$data$simulations <- sim_5
  
  # Apply FEV-bias
  lWeights_ctrl <- BT5$output[ grepl("grw", names(BT5$output)) ]
  lWeights_fev_ctrl <- lapply( lWeights_ctrl, FUN = function(wmat) { fevBias( x = wmat, q = 6 ) } )
  
  # Project FEV-biased weights to boundary of feasible set
  lWeights_cons_ctrl <- lWeights_ctrl
  for ( i in 1:length(lWeights_cons_ctrl) ) {
    for ( today in rownames(lWeights_cons_ctrl[[i]]) ) {
      RBE <- rebalenv( rebdate = today )
      # Cons <- getConstraints(RBE$GPS)
      lincon <- getConstraints( RBE$GPS, "bounds" )
      upper_bounds <- lincon$upper
      w_unc <- lWeights_fev_ctrl[[i]][today, ]
      # idx_notNA <- which(!is.na(w_unc))
      idx_notNA <- which( !is.na( lWeights_ctrl[[i]][today, ] ) )
      w_unc <- w_unc[idx_notNA]
      w_unc <- setNames( as.numeric(w_unc), names(w_unc))
      # debugonce( map2boxcon )
      lWeights_cons_ctrl[[i]][today, idx_notNA] <- map2boxcon( w_unc = w_unc,
                                                               upper = upper_bounds,
                                                               itermax = 10^3 )
    }
  }
  BT5_fev2cons <- BT5$copy()
  BT5_fev2cons$spec$name <- paste0(BT5$spec$name, "_fev2cons")
  BT5_fev2cons$output <- lWeights_cons_ctrl
  
  # Simulate
  tmp <- lapply( BT5_fev2cons$output, FUN = simFUN )
  sim_5_cons <- do.call( cbind, tmp )
  sim_5_cons <- sim_5_cons[isWeekday(time(sim_5_cons)), ]
  colnames(sim_5_cons) <- paste0("rp_grw_uniform_cons_sizestdzctrl", 1:ncol(sim_5_cons))
  BT5_fev2cons$data$simulations <- sim_5_cons
  
  # # Save
  # BT5$save( wd = paste0(wd_data, "waRehouse/"), name = BT5$spec$name, without_data = FALSE )
  # BT5_fev2cons$save( wd = paste0(wd_data, "waRehouse/"), name = BT5_fev2cons$spec$name, without_data = FALSE )
  
  
  # Load
  BT5 <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "grw_multifac_usa_sizestdzctrl" )
  BT5_fev2cons <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "grw_multifac_usa_sizestdzctrl_fev2cons" )
  sim_5 <- BT5$data$simulations
  sim_5_cons <- BT5_fev2cons$data$simulations
  
  
  
  lStats_sim5 <- descStats( window( sim_5, start(sim_5), end_date), descStatsSpec(annualize = TRUE) )
  lStats_sim5_cons <- descStats( window( sim_5_cons, start(sim_5_cons), end_date), descStatsSpec(annualize = TRUE) )
  
  
  plot( x = lStats_sim1$stats["sds", ], y = lStats_sim1$stats["means", ], pch = 19 )
  points( x = lStats_sim2$stats["sds", ], y = lStats_sim2$stats["means", ], pch = 19, col = 2 )
  points( x = lStats_sim3$stats["sds", ], y = lStats_sim3$stats["means", ], pch = 19, col = 3 )
  points( x = lStats_sim4$stats["sds", ], y = lStats_sim4$stats["means", ], pch = 19, col = 4 )
  points( x = lStats_sim5$stats["sds", ], y = lStats_sim5$stats["means", ], pch = 19, col = 5 )

  plot( x = lStats_sim1_cons$stats["sds", ], y = lStats_sim1_cons$stats["means", ], pch = 19 )
  points( x = lStats_sim2_cons$stats["sds", ], y = lStats_sim2_cons$stats["means", ], pch = 19, col = 2 )
  points( x = lStats_sim3_cons$stats["sds", ], y = lStats_sim3_cons$stats["means", ], pch = 19, col = 3 )
  points( x = lStats_sim4_cons$stats["sds", ], y = lStats_sim4_cons$stats["means", ], pch = 19, col = 4 )
  points( x = lStats_sim5_cons$stats["sds", ], y = lStats_sim5_cons$stats["means", ], pch = 19, col = 5 )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BACKTEST UNIFORM RP'S WITHIN CONSTRAINTS WITH VALUE, QUALITY, SIZE STANDARDIZED CONTROL
  # --------------------------------------------------------------------------
  
  BT6 <- BT_Fctrl$copy()
  BT6$spec$factor_names_constraints <- c("value", "quality", "size_stdz")
  BT6$spec$keep_all_rebdates <- FALSE
  BT6$spec$portfolio <- "grw"
  BT6$spec$n_sim <- 10^2 * 5
  BT6$spec$name <- "grw_multifac_usa_valqualsizectrl"
  BT6$output <- list()
  BT6$runLoop()
  
  # Simulate
  tmp <- lapply( BT6$output[ grepl("grw", names(BT6$output)) ], FUN = simFUN )
  sim_6 <- do.call( cbind, tmp )
  sim_6 <- sim_6[isWeekday(time(sim_6)), ]
  colnames(sim_6) <- paste0("rp_grw_uniform_valqualsizectrl", 1:ncol(sim_6))
  BT6$data$simulations <- sim_6
  
  # Apply FEV-bias
  lWeights_ctrl <- BT6$output[ grepl("grw", names(BT6$output)) ]
  lWeights_fev_ctrl <- lapply( lWeights_ctrl, FUN = function(wmat) { fevBias( x = wmat, q = 6 ) } )
  
  # Project FEV-biased weights to boundary of feasible set
  lWeights_cons_ctrl <- lWeights_ctrl
  for ( i in 1:length(lWeights_cons_ctrl) ) {
    for ( today in rownames(lWeights_cons_ctrl[[i]]) ) {
      RBE <- rebalenv( rebdate = today )
      # Cons <- getConstraints(RBE$GPS)
      lincon <- getConstraints( RBE$GPS, "bounds" )
      upper_bounds <- lincon$upper
      w_unc <- lWeights_fev_ctrl[[i]][today, ]
      # idx_notNA <- which(!is.na(w_unc))
      idx_notNA <- which( !is.na( lWeights_ctrl[[i]][today, ] ) )
      w_unc <- w_unc[idx_notNA]
      w_unc <- setNames( as.numeric(w_unc), names(w_unc))
      # debugonce( map2boxcon )
      lWeights_cons_ctrl[[i]][today, idx_notNA] <- map2boxcon( w_unc = w_unc,
                                                               upper = upper_bounds,
                                                               itermax = 10^3 )
    }
  }
  BT6_fev2cons <- BT6$copy()
  BT6_fev2cons$spec$name <- paste0(BT6$spec$name, "_fev2cons")
  BT6_fev2cons$output <- lWeights_cons_ctrl
  
  # Simulate
  tmp <- lapply( BT6_fev2cons$output, FUN = simFUN )
  sim_6_cons <- do.call( cbind, tmp )
  sim_6_cons <- sim_6_cons[isWeekday(time(sim_6_cons)), ]
  colnames(sim_6_cons) <- paste0("rp_grw_uniform_cons_valqualsizectrl", 1:ncol(sim_6_cons))
  BT6_fev2cons$data$simulations <- sim_6_cons
  
  # # Save
  # BT6$save( wd = paste0(wd_data, "waRehouse/"), name = BT6$spec$name, without_data = FALSE )
  # BT6_fev2cons$save( wd = paste0(wd_data, "waRehouse/"), name = BT6_fev2cons$spec$name, without_data = FALSE )
  
  
  # Load
  BT6 <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "grw_multifac_usa_valqualsizectrl" )
  BT6_fev2cons <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "grw_multifac_usa_valqualsizectrl_fev2cons" )
  sim_6 <- BT6$data$simulations
  sim_6_cons <- BT6_fev2cons$data$simulations
  
  
  
  lStats_sim6 <- descStats( window( sim_6, start(sim_6), end_date), descStatsSpec(annualize = TRUE) )
  lStats_sim6_cons <- descStats( window( sim_6_cons, start(sim_6_cons), end_date), descStatsSpec(annualize = TRUE) )
  
  
  plot( x = lStats_sim1$stats["sds", ], y = lStats_sim1$stats["means", ], pch = 19 )
  points( x = lStats_sim2$stats["sds", ], y = lStats_sim2$stats["means", ], pch = 19, col = 2 )
  points( x = lStats_sim3$stats["sds", ], y = lStats_sim3$stats["means", ], pch = 19, col = 3 )
  points( x = lStats_sim4$stats["sds", ], y = lStats_sim4$stats["means", ], pch = 19, col = 4 )
  points( x = lStats_sim5$stats["sds", ], y = lStats_sim5$stats["means", ], pch = 19, col = 5 )
  points( x = lStats_sim6$stats["sds", ], y = lStats_sim6$stats["means", ], pch = 19, col = 6 )
  
  plot( x = lStats_sim1_cons$stats["sds", ], y = lStats_sim1_cons$stats["means", ], pch = 19 )
  points( x = lStats_sim2_cons$stats["sds", ], y = lStats_sim2_cons$stats["means", ], pch = 19, col = 2 )
  points( x = lStats_sim3_cons$stats["sds", ], y = lStats_sim3_cons$stats["means", ], pch = 19, col = 3 )
  points( x = lStats_sim4_cons$stats["sds", ], y = lStats_sim4_cons$stats["means", ], pch = 19, col = 4 )
  points( x = lStats_sim5_cons$stats["sds", ], y = lStats_sim5_cons$stats["means", ], pch = 19, col = 5 )
  points( x = lStats_sim6_cons$stats["sds", ], y = lStats_sim6_cons$stats["means", ], pch = 19, col = 6 )
  
  
  
  
  
  
  