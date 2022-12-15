
  
  ############################################################################
  ### RANDOM PORTFOLIOS - ALPHABET STRATEGY RANDOM PORTOLIOS - USA
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     28.07.2022
  # First version:    28.07.2022
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
  require(Backtest)
  require(RP)
  wd <- "R:/Asset_Management/Research_Projects/External_Research_Projects/GeomScale/Random_Portfolios/"
  source( paste0(wd, "Source/class_BacktestCustom.R") )
  source( paste0(wd, "Source/custom_functions.R") )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # PARAMETERS
  # --------------------------------------------------------------------------
  
  universe <- "usa"
  
  
  # --------------------------------------------------------------------------
  # LOAD DATA AND INITIALIZE REBALANCING ENVIRONMENTS
  # --------------------------------------------------------------------------
  
  BT0 <- loadBacktest( wd = paste0(wd, "Data/"),
                       name = universe )
  BT0$initRBE()
  
  
  # --------------------------------------------------------------------------
  # PREPARE BACKTEST OBJECT
  # --------------------------------------------------------------------------
  
  BT <- BacktestCustom$new()
  BT$spec <- BT0$spec
  BT$spec$selection_filter <- c("bm")
  BT$spec$cons_fun <- NULL
  BT$spec$GPS@covariance$method <- "duv"
  BT$spec$portfolio <- "score"
  BT$data <- BT0$data


  
  # --------------------------------------------------------------------------
  # BACKTEST EQW PORTFOLIO
  # --------------------------------------------------------------------------
  
  BT_eqw <- BT$copy()
  BT_eqw$spec$portfolio <- NULL
  BT_eqw$spec$GPS@covariance$method <- "duv"
  BT_eqw$runLoop()
  sim_eqw <- simPortfolio( X = BT_eqw$data$X_sim,
                           wghts = BT_eqw$output$weights, 
                           language = "C" )
  
  
  # --------------------------------------------------------------------------
  # BACKTEST Min Var PORTFOLIO
  # --------------------------------------------------------------------------

  # BT_mv <- BacktestPanel$new()
  # BT_mv$spec <- BT$spec
  # BT_mv$spec$portfolio <- "minvariance"
  # BT_mv$runLoop()
  # sim_mv <- simPortfolio( X = BT$data$X_sim,
  #                         wghts = BT_mv$output$weights, 
  #                         language = "C" )
  
  
  
  
  # --------------------------------------------------------------------------
  # LETTER SCORES
  # --------------------------------------------------------------------------
  
  Names <- BT0$data$name_mat
  letter_vec <- substr( aviationAlphabet(), 1, 1 )
  lScore <- list()
  names_vec <- toupper(as.vector(Names))
  for ( i in seq(along = letter_vec ) ) {
    score_vec <- unlist( lapply( names_vec, letterScore, preference = letter_vec[i] ) )
    Score <- timeSeries( matrix( score_vec,
                                 nrow = nrow(Names), ncol = ncol(Names),
                                 dimnames = list(rownames(Names), colnames(Names)) ),
                         time(Names) )
    lScore[[i]] <- Score
  }
  names(lScore) <- letter_vec
  
  headleft(Names)
  headleft(Score)    
  range(Score)
  
  
  
  
  # --------------------------------------------------------------------------
  # RANDOM PORTFOLIOS ON SCORE SUBSETS
  # --------------------------------------------------------------------------
  
  llStats <- list()
  llStats_inv <- list()
  lSim <- list()
  lSim_inv <- list()
  lPortf_scores <- list()
  lWeights <- list()
  lWeights_inv <- list()
  
  for ( i in seq(along = letter_vec) ) {
  # for ( i in c(1, 26) ) {
    
    BT1 <- BT$copy()
    BT1$spec$portfolio <- "score_rp"
    BT1$spec$sampling_method = "dirichlet"
    # BT1$spec$sampling_method = "moon"
    # BT1$spec$sampling_dist = "uniform"
    BT1$spec$sampling_dist = "capw"
    BT1$data$Score <- lScore[[i]]
    # BT1$spec$n_sim <- 10^3
    # debugonce( BT1$score_rpPortfolio )
    BT1$runLoop()  
    
    # Simulate
    simFUN <- function(wmat) 
    { 
      simPortfolio( X = BT0$data$X_sim, 
                    wghts = wmat,
                    language = "C" )
    }
    tmp <- lapply( BT1$output, FUN = simFUN )
    sim_tmp <- do.call( cbind, tmp )
    lStats <- descStats( sim_tmp, descStatsSpec(annualize = TRUE) )
    
    # Portfolio letter scores
    portf_scores <- do.call( cbind, portfolioScore( BT1$output, score_mat = Score ) )
    
    lWeights[[i]] <- BT1$output
    lSim[[i]] <- sim_tmp
    llStats[[i]] <- lStats
    lPortf_scores[[i]] <- portf_scores
    # lWeights_inv[[i]] <- BT1$output$weights_inv
    # lSim_inv[[i]] <- sim_inv_tmp
    # llStats_inv[[i]] <- lStats_inv
    
  }
  names(lWeights) <- letter_vec
  names(lSim) <- letter_vec
  names(llStats) <- letter_vec
  names(lPortf_scores) <- letter_vec
  # names(lWeights_inv) <- letter_vec
  # names(lSim_inv) <- letter_vec
  # names(llStats_inv) <- letter_vec
  
  
  
  
  
  # statsolz:::plot.stats( lStats, sortby = NULL )
  # mean(lStats$stats["means", ])
  # mean(lStats$stats["sds", ])
  # head(colnames(lStats$stats))
  # 
  # 
  # 
  # 
  # wmat <- lWeights[["Z"]][[1]]
  # a <- apply( wmat, 1, function(x) { sum(na.omit(x) > 0) } )
  # b <- apply( lScore[["Z"]][rownames(wmat), ], 1, function(x) { sum(na.omit(x) > 0) } )
  # cbind(a, b)
  
  
  
  
  # --------------------------------------------------------------------------
  # RANDOM PORTFOLIOS ON FIXED SUBSET
  # --------------------------------------------------------------------------
  
  samples <- rdirichlet( n = nrow(Score), alpha = rep(1, ncol(Score) ) )
  samples <- timeSeries( samples, time(Score) )
  samples_fev <- fevBias( samples, q = 5 )
  samples_fev[ samples_fev < 0.005 ] <- 0
  samples_fev <- timeSeries( samples_fev, time(Score) )
  colnames(samples_fev) <- colnames(Score)
  
  # weightsBarPlot( samples[1:100, ] )
  # weightsBarPlot( samples_fev[1:100, ] )
  
  apply(samples, 1, function(x) { sum(x > 0)})
  apply(samples_fev, 1, function(x) { sum(x > 0)})
  
  
  
  
  
  BT1 <- BT$copy()
  BT1$spec$portfolio <- "score_rp"
  BT1$spec$sampling_method = "dirichlet"
  # BT1$spec$sampling_method = "moon"
  # BT1$spec$sampling_dist = "uniform"
  BT1$spec$sampling_dist = "capw"
  BT1$data$Score <- samples_fev
  BT1$spec$n_sim <- 10^3
  # debugonce( BT1$score_rpPortfolio )
  BT1$runLoop()  
  
  # Simulate
  simFUN <- function(wmat) 
  { 
    simPortfolio( X = BT0$data$X_sim, 
                  wghts = wmat,
                  language = "C" )
  }
  tmp <- lapply( BT1$output, FUN = simFUN )
  sim_tmp <- do.call( cbind, tmp )
  lStats <- descStats( sim_tmp, descStatsSpec(annualize = TRUE) )
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # SAVE
  # --------------------------------------------------------------------------
  
  env <- new.env()
  env$lSim <- lSim
  env$lWeights <- lWeights
  env$llStats <- llStats
  env$lPortf_scores <- lPortf_scores
  # saveRDS( env, file = paste0(wd, "waRehouse/alphabet_strategy_rps.rds") )
  saveRDS( env, file = paste0(wd, "waRehouse/alphabet_strategy_rps_bmcentered.rds") )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # LOAD
  # --------------------------------------------------------------------------

  # env <- readRDS( file = paste0(wd, "waRehouse/alphabet_strategy_rps.rds") )
  env <- readRDS( file = paste0(wd, "waRehouse/alphabet_strategy_rps_bmcentered.rds") )
  lSim <- env$lSim
  lWeights <- env$lWeights
  llStats <- env$llStats
  lPortf_scores <- env$lPortf_scores
  
  
  
  
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
  
 
  
  
  
  
  
  
  
  # # Using LogConcDEAD
  # # install.packages("LogConcDEAD")
  # require(LogConcDEAD)
  # 
  # 
  # tic <- Sys.time()
  # A <- as.matrix( t( lStats$stats[c("sds", "means"), ] ) )
  # LCD <- mlelcd( x = A )
  # (toc <- Sys.time() - tic)
  # 
  # 
  # plot( LCD, addp = TRUE, uselog = TRUE, type = "ic", drawlabels = FALSE )
  # points( x = lStats$stats["sds", (ncol(lStats$stats)-27):(ncol(lStats$stats)-1)],
  #         y = lStats$stats["means", (ncol(lStats$stats)-27):(ncol(lStats$stats)-1)],
  #         col = "grey", cex = 2, pch = 19 )
  # points( x = lStats$stats["sds", 2],
  #         y = lStats$stats["means", 2], col = "blue", cex = 2, pch = 19 )
  # points( x = lStats$stats["sds", 1],
  #         y = lStats$stats["means", 1], col = "red", cex = 2, pch = 19 )
  
  
  
  
  
  
  
  
  
  
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
  
  
  
  
  
  
  
  
  
  
  
  
      
  