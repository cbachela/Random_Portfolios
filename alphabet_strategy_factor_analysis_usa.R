
  
  ############################################################################
  ### RANDOM PORTFOLIOS - ALPHABET STRATEGY - USA
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     10.06.2022
  # First version:    10.06.2022
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
  # LOAD
  # --------------------------------------------------------------------------
  
  env <- readRDS( file = paste0(wd, "waRehouse/alphabet_strategy.rds") )
  sim_eqw <- env$sim_eqw
  lSim <- env$lSim
  lSim_inv <- env$lSim_inv
  lWeights <- env$lWeights
  lWeights_inv <- env$lWeights_inv
  llStats <- env$llStats
  llStats_inv <- env$llStats_inv
  lPortf_scores <- env$lPortf_scores
  
  
  
  # --------------------------------------------------------------------------
  # SIMULATIONS OF LETTER-SCORE PORTFOLIOS
  # --------------------------------------------------------------------------
  
  # Full period
  end_date <- "2022-03-06"
  sim_sort <- do.call( cbind, lSim )
  sim_sort <- window( sim_sort, start(sim_sort), end_date )
  colnames(sim_sort) <- letter_vec
  sim_tmp <- na.omit( cbind( sim_sort, bm = BT0$data$X_bm, eqw = sim_eqw ) )
  lStats <- descStats( sim_tmp )
  
  sim_sort_inv <- do.call( cbind, lSim_inv )
  sim_sort_inv <- window( sim_sort_inv, start(sim_sort_inv), end_date )
  colnames(sim_sort_inv) <- paste0("|", letter_vec)
  sim_tmp_inv <- na.omit( cbind( sim_sort_inv, bm = BT0$data$X_bm, eqw = sim_eqw ) )
  
  # First decade
  sim_tmp_decade1 <- window( sim_tmp, start(sim_sort), "2011-01-01" )
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
  # FF REGRESSION
  # --------------------------------------------------------------------------
  
  path_ff <- "R:/Asset_Management/R/myRData/Factor/"
  env_ff <- readRDS( file = paste0(path_ff, "data.Rds") )
  FF3 <- env_ff$lFactor$`3F`$USA_daily
  FF5 <- env_ff$lFactor$`5F`$USA_daily
  
  dates <- intersect(rownames(FF5), rownames(sim_tmp))
  X_train <- FF5[dates, 1:5]
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
  
  
  for ( strat_name in colnames(sim_sort) ) {
    
    reg <- regression( Y_train = sim_sort[dates, strat_name], 
                       X_train = X_train, 
                       type = "ols" )
    lReg[[strat_name]] <- reg
    lCoeff[[strat_name]] <- reg$coeffmat
    
    reg_inv <- regression( Y_train = sim_sort_inv[dates, paste0("|", strat_name)], 
                           X_train = X_train, 
                           type = "ols" )
    lReg_inv[[strat_name]] <- reg_inv
    lCoeff_inv[[strat_name]] <- reg_inv$coeffmat
    
    reg_decade1 <- regression( Y_train = sim_sort[dates, strat_name], 
                               X_train = X_train_decade1, 
                               type = "ols" )
    lReg_decade1[[strat_name]] <- reg_decade1
    lCoeff_decade1[[strat_name]] <- reg_decade1$coeffmat
    
    reg_decade1_inv <- regression( Y_train = sim_sort_inv[dates, paste0("|", strat_name)], 
                                   X_train = X_train_decade1, 
                                   type = "ols" )
    lReg_decade1_inv[[strat_name]] <- reg_decade1_inv
    lCoeff_decade1_inv[[strat_name]] <- reg_decade1_inv$coeffmat
    
    reg_decade2 <- regression( Y_train = sim_sort[dates, strat_name], 
                               X_train = X_train_decade2, 
                               type = "ols" )
    lReg_decade2[[strat_name]] <- reg_decade2
    lCoeff_decade2[[strat_name]] <- reg_decade2$coeffmat
    
    reg_decade2_inv <- regression( Y_train = sim_sort_inv[dates, paste0("|", strat_name)], 
                                   X_train = X_train_decade2, 
                                   type = "ols" )
    lReg_decade2_inv[[strat_name]] <- reg_decade2_inv
    lCoeff_decade2_inv[[strat_name]] <- reg_decade2_inv$coeffmat
    
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
  
  
  
  
  # --------------------------------------------------------------------------
  # SAVE
  # --------------------------------------------------------------------------
  
  env <- new.env()
  env$alpha_score <- alpha_score
  env$beta_mat <- beta_mat
  env$beta_pval_mat <- beta_pval_mat
  env$p_values <- p_values
  saveRDS( env, file = paste0(wd, "waRehouse/alphabet_strategy_factor_analysis.rds") )
  
  
  
  
  ###
  
  barplot( alpha_score )
  
  barplot( sort(p_values) )
  abline( h = 0.05 )
  sum( p_values < 0.05 ) / length(p_values)
  sum( p_values < 0.1 ) / length(p_values)
  
  
  barplot( beta_pval_mat, beside = TRUE, col = 1:nrow(beta_pval_mat) )
  abline( h = 0.05 )
  
  
  ###
  
  
  
  
  
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
  
  
  
  
  
  
  
  
  
  
  
    
  
  

  
  
  
  
  
  
      
  