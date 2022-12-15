
  
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
  require(BBSolz)
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
  rebdates <- rownames(BT0$data$sel_mat)
  BT$spec$rebdates <- rebdates
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
                           fc = 0, 
                           vc = 0 )

  
  
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
  
  
  
  # A closer look at Z
  
  names_vec_unique <- unique(names_vec)
  idx_z <- which(grepl("z", tolower(names_vec_unique)))
  names_vec_unique[idx_z]
  
  Score <- lScore[["Z"]]
  idx <- which( Score >= 1, arr.ind = TRUE )
  Names[ ,unique(idx[ ,2])]
  
  series_id <- colnames(Names)[unique(idx[ ,2])]
  X_Z <- BT0$data$X_sim[ ,series_id]
  dim(X_Z)
  descStats(X_Z)
  
  plot( X_Z[ ,1:10] )
  Names[nrow(Names), series_id]
  
  mu <- apply( X_Z[isWeekday(time(X_Z)), ], 2, function(x) { meanGeo(na.omit(x), scalefactor = 252) } )
  mu <- setNames( as.numeric(mu), Names[nrow(Names), series_id] )
  mu
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # SCORE PORTFOLIOS
  # --------------------------------------------------------------------------
  
  llStats <- list()
  llStats_inv <- list()
  lSim <- list()
  lSim_inv <- list()
  lPortf_scores <- list()
  lWeights <- list()
  lWeights_inv <- list()
  fc <- 0
  vc <- 0
  
  for ( i in seq(along = letter_vec) ) {
    
    BT1 <- BT$copy()
    BT1$spec$portfolio <- "score"
    BT1$data$Score <- lScore[[i]]
    # debugonce( BT1$scorePortfolio )
    BT1$runLoop()  
    
    # Simulate
    sim_tmp <- simPortfolio( X = BT0$data$X_sim, 
                             wghts = BT1$output$weights,
                             fc = fc, 
                             vc = vc )
    lStats <- descStats( sim_tmp )
    sim_inv_tmp <- simPortfolio( X = BT0$data$X_sim, 
                                 wghts = BT1$output$weights_inv,
                                 fc = fc, 
                                 vc = vc )
    lStats_inv <- descStats( sim_inv_tmp )
    
    # Portfolio letter scores
    portf_scores <- do.call( cbind, portfolioScore( BT1$output, score_mat = Score ) )
    
    lWeights[[i]] <- BT1$output$weights
    lWeights_inv[[i]] <- BT1$output$weights_inv
    lSim[[i]] <- sim_tmp
    lSim_inv[[i]] <- sim_inv_tmp
    llStats[[i]] <- lStats
    llStats_inv[[i]] <- lStats_inv
    lPortf_scores[[i]] <- portf_scores
  
  }
  names(lWeights) <- names(lWeights_inv) <- letter_vec
  names(lSim) <- names(lSim_inv) <- letter_vec
  names(llStats) <- names(llStats_inv) <- letter_vec
  names(lPortf_scores) <- letter_vec
  
  
  
  # --------------------------------------------------------------------------
  # SAVE
  # --------------------------------------------------------------------------
  
  env <- new.env()
  env$sim_eqw <- sim_eqw
  env$lSim <- lSim
  env$lSim_inv <- lSim_inv
  env$lWeights <- lWeights
  env$lWeights_inv <- lWeights_inv
  env$llStats <- llStats
  env$llStats_inv <- llStats_inv
  env$lPortf_scores <- lPortf_scores
  env$lScore <- lScore
  saveRDS( env, file = paste0(wd, "waRehouse/alphabet_strategy.rds") )
  
  
  
  
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
  lScore <- env$lScore
  
  
  
  # --------------------------------------------------------------------------
  # SIMULATIONS OF LETTER-SCORE PORTFOLIOS
  # --------------------------------------------------------------------------
  
  # Full period
  end_date <- "2022-03-06"
  sim_sort <- do.call( cbind, lSim )
  sim_sort <- window( sim_sort, start(sim_sort), end_date )
  colnames(sim_sort) <- letter_vec
  sim_tmp <- na.omit( cbind( sim_sort, bm = BT0$data$X_bm, eqw = sim_eqw ) )
  sim_tmp <- sim_tmp[isWeekday(time(sim_tmp)), ]
  lStats <- descStats( sim_tmp )
  
  sim_sort_inv <- do.call( cbind, lSim_inv )
  sim_sort_inv <- window( sim_sort_inv, start(sim_sort_inv), end_date )
  colnames(sim_sort_inv) <- paste0("|", letter_vec)
  sim_tmp_inv <- na.omit( cbind( sim_sort_inv, bm = BT0$data$X_bm, eqw = sim_eqw ) )
  sim_tmp_inv <- sim_tmp_inv[isWeekday(time(sim_tmp_inv)), ]
  lStats_inv <- descStats( sim_tmp_inv )
  
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
  
  
  
  # Bear Bull Performance
  
  BBSObj <- BBSRC$new()
  BBSObj$setCtrl()
  BBSObj$data$X_level <- cumulated( BT0$data$X_bm, "discrete" )
  BBSObj$runRobust()
  BBSObj$data$X <- sim_tmp
  phase_stats <- BBSObj$phaseStats()
  barplot( t( cbind( bull = phase_stats$`states==1`["means", ],
                     bear = phase_stats$`states==-1`["means", ] ) ),
           beside = TRUE, col = c("green", "red") )
  
  
  
  
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
 
  
  
  
  