
  
  ############################################################################
  ### RANDOM PORTFOLIOS - FIND PONZI STRATEGY - USA
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
  BT$spec$selection_filter <- c(BT0$spec$selection_filter, "NA")
  BT$data <- BT0$data
  BT$data$upper_mat <- BT$data$upper_mat * 0 + 1
  BT$spec$fev_bias <- 2:10
  BT$spec$n_sim <- 10
  BT$spec$portfolio <- "dirichlet"
  BT$spec$sharpe_filter_th <- NULL
  
  
  # --------------------------------------------------------------------------
  # BACKTEST EQW PORTFOLIO
  # --------------------------------------------------------------------------
  
  BT_eqw <- BT$copy()
  BT_eqw$spec$GPS@covariance$method <- "duv"
  BT_eqw$spec$portfolio <- NULL
  BT_eqw$runLoop()
  # sim_eqw <- BT_eqw$simulate( language = "R" )
  sim_eqw <- simPortfolio( X = BT_eqw$data$X_sim, wghts = BT_eqw$output$weights, language = "R" )

  
  # --------------------------------------------------------------------------
  # BACKTEST DIRICHLET PORTFOLIOS
  # --------------------------------------------------------------------------
  
  BT1 <- BT$cop()
  BT1$runLoop()
    
  # Simulate
  lSim <- lapply( BT1$output, FUN = function(wmat) { simPortfolio( X = BT1$data$X_sim, 
                                                                   wghts = wmat, language = "R" ) } )
  sim <- do.call( cbind, lSim )
  
  
  
  letterScore <- function( ticker, preference = "A" ) 
  {
    ans <- 0
    if ( !is.na(ticker) ) {
      ans <- str_count( ticker, preference )
    }
    return( ans )
  }
  Names <- BT0$data$name_mat
  letter_vec <- substr( aviationAlphabet(), 1, 1)
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
  
  

  
  
  lPortf_score <- list()
  for ( i in 1:length(lScore) )  {
    FUN <- function(wmat) {
      esgScores( wmat = wmat, esg_mat = Score )
    }    
    Score <- lScore[[i]]
    lPortf_score[[i]] <- do.call( cbind, lapply( BT1$output, FUN = FUN ) )
  }
  names(lPortf_score) <- names(lScore)
  
  lportf_score <- lapply( lPortf_score, FUN = function(x) { apply(x, 2, mean) } )
  
  

  
  tmp <- unlist( lapply( lportf_score, FUN = mean ) ) 
  barplot(tmp)  
  
  
  tmp <- unlist( lapply( lportf_score, FUN = function(x) { cor(x, mu) } ) )
  idx <- which(tmp == max(tmp))
  tmp[idx]
  
  mu <- meanGeo( sim )
  plot( x = lportf_score[["K"]], y = mu )

  
  
  
  #####################
  
  
  w_array <- list2array( BT1$output )
  Score <- lScore[["K"]]
  # debugonce( simSortedBy )
  sim_A <- simSortedBy( sim = sim,
                        BT0 = BT0,
                        Score = Score,
                        w_array = w_array )
  
  colors <- fBasics::divPalette(n = ncol(sim), "RdYlGn" )
  plot( as.simTS(sim_A), col = colors )
  
  barplot( meanGeo(sim) )
  barplot( meanGeo(sim_A) )
  
  
  
  colors <- fBasics::divPalette(n = ncol(sim), "RdYlGn" )
  plot( as.simTS(na.omit(cbind(sim_A, BT0$data$X_bm))), col = c(colors, 1) )

  lStats <- descStats( sim_A )  
  plot.stats( lStats, sortby = NULL, col = colors )  
  
  plot_desc_stats( lStats )
  
  plot( lStats$stats["sds", ], lStats$stats["means", ], col = colors, pch = 19, cex = 2 )
  barplot( lStats$stats["sharpe", ] )  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # QUINTILE PORTFOLIOS
  # --------------------------------------------------------------------------
  
  llStats <- list()
  llSim <- list()
  llPortf_scores <- list()
  
  for ( i in seq(along = letter_vec) ) {
    
    # Run backtest
    BT2 <- BT$copy()
    BT2$spec$portfolio <- "quintile"
    BT2$spec$weighting_method <- "eqw"
    BT2$data$Score <- lScore[[i]]
    # debugonce( BT2$quantilePortfolio )
    # debognce( get_sequence_of_scores )
    BT2$runLoop()
    
    
    # Simulate
    lSim_deciles <- lapply( BT2$output[-1], FUN = function(wmat) { simPortfolio( X = BT0$data$X_sim, 
                                                                                 wghts = wmat,
                                                                                 language = "R" ) } )
    sim_deciles <- do.call( cbind, lSim_deciles )
    lStats <- descStats( na.omit(cbind(sim_deciles, BT0$data$X_bm)) )

    # Portfolio scores
    FUN <- function(wmat) {
      esgScores( wmat = wmat, esg_mat = Score )
    }    
    portf_scores <- do.call( cbind, lapply( BT2$output[-1], FUN = FUN ) )
    
    llSim[[i]] <- sim_deciles
    llStats[[i]] <- lStats
    llPortf_scores[[i]] <- portf_scores
    
  }
  
  
  i <- 1
  sim_deciles <- llSim[[i]]
  lStats <- llStats[[i]]
  portf_scores <- lPortf_scores[[i]]
  
  colors <- fBasics::divPalette( n = ncol(sim_deciles), "RdYlGn" )
  plot( as.simTS( na.omit(cbind(sim_deciles, BT0$data$X_bm)) ), col = c(colors, 1) )  
  
  barplot( lStats$stats["sharpe", ], col = c(colors, 1) )
  
  plot( x = lStats$stats["sds", ], y = lStats$stats["means", ], 
        col = c(colors, 1), pch = 19, cex = 2 )  
  
  plot( portf_scores, plot.type = "single" )
  
  boxplot( as.data.frame(portf_scores) )
  
  
  
  
  # --------------------------------------------------------------------------
  # SCORE PORTFOLIOS
  # --------------------------------------------------------------------------
  
  llStats_2 <- list()
  llSim_2 <- list()
  llPortf_scores_2 <- list()
  lWeights_2 <- list()
  
  for ( i in seq(along = letter_vec) ) {
    
    BT3 <- BT$copy()
    BT3$spec$portfolio <- "score"
    BT3$data$Score <- lScore[[i]]
    # debugonce( BT3$scorePortfolio )
    BT3$runLoop()  
    
    # Weights
    lWeights_2[[i]] <- BT3$output$weights
    
    # Simulate
    sim_tmp <-simPortfolio( X = BT0$data$X_sim, 
                            wghts = BT3$output$weights,
                            language = "R" )
    lStats <- descStats( sim_tmp )
    
    # Portfolio scores
    FUN <- function(wmat) {
      esgScores( wmat = wmat, esg_mat = Score )
    }    
    portf_scores <- do.call( cbind, lapply( BT3$output, FUN = FUN ) )
    
    llSim_2[[i]] <- sim_tmp
    llStats_2[[i]] <- lStats
    llPortf_scores_2[[i]] <- portf_scores
  
  }
  
  
  
  sim_sort <- do.call( cbind, llSim_2 )
  colnames(sim_sort) <- letter_vec
  sim_tmp <- na.omit( cbind( sim_sort, BT0$data$X_bm, sim_eqw ) )
  
  lStats <- descStats( sim_tmp )
  lStats
  
  colors <- fBasics::divPalette(n = ncol(sim_sort), "RdYlGn" )
  
  plot( as.simTS( sim_tmp ), col = c(colors, 1, 4) )
  plot( lStats$stats["sds", ], lStats$stats["means", ], pch = 19, cex = 2, col = c(colors, 1, 4) )
  barplot( lStats$stats["sharpe", ], col = c(colors, 1, 4) )
  
  
  
  weightsBarPlot( lWeights_2[["Z"]] )
  weightsBarPlot( lWeights_2[["A"]] )
  
  
  # Outperformance
  sim_outperf <- simOutperformance( x = sim_tmp[ ,"Z"], y = sim_tmp[ ,"bm"] )
  plot( as.simTS(sim_outperf) )
  descStats( sim_outperf )
  
  
  # Sector allocation
  sector_alloc <- groupWeights( wmat = lWeights[["Z"]], group_mat = BT0$data$sector_mat )
  weightsBarPlot( sector_alloc )
  boxplot( as.data.frame(sector_alloc) )
  
  
  # Sharpe test
  # X <- sim_tmp[ ,c("Z", "bm")]
  X <- sim_tmp
  ST <- sharpeTest( X = X, method = "HAC_PW" )
  ST
 
  
  # FF Regression
  path_ff <- "R:/Asset_Management/R/myRData/Factor/"
  env_ff <- readRDS( file = paste0(path_ff, "data.Rds") )
  FF3 <- env_ff$lFactor$`3F`$USA_daily
  FF5 <- env_ff$lFactor$`5F`$USA_daily
  
  # dates <- intersect(rownames(FF3), rownames(sim_tmp))
  # X_train <- FF3[dates, 1:3]
  dates <- intersect(rownames(FF5), rownames(sim_tmp))
  X_train <- FF5[dates, 1:5]
  lReg <- list()
  lCoeff <- list()
  for ( j in 1:ncol(sim_sort) ) {
    reg <- regression( Y_train = sim_sort[dates, j], 
                       X_train = X_train, 
                       type = "ols" )
    lReg[[j]] <- reg
    lCoeff[[j]] <- reg$coeffmat
  }  
  alpha_score <- unlist( lapply( lCoeff, FUN = function(x) { x[1, 1] } ) )
  beta_mat <- do.call( cbind, lapply( lCoeff, FUN = function(x) { x[-1, 1] } ) )
  beta_pval_mat <- do.call( cbind, lapply( lCoeff, FUN = function(x) { x[-1, 4] } ) )
  p_values <- unlist( lapply( lCoeff, FUN = function(x) { x[1, 4] } ) )
  alpha_mat <- cbind( alpha_score, p_values )
  
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
  
  
  barplot( beta_mat, beside = TRUE, col = 1:nrow(beta_mat) )
  
  
  
    
  
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
  
  
  
  
  
  
  
  
  
  
  
  
      
  