
  
  ############################################################################
  ### INDEX TRACKING - REPLICATE LOW SIZE PORTFOLIO WITH LARGE COMPANIES
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     02.02.2023
  # First version:    02.02.2023
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
  source("H:/R/github/Random_Portfolios/Source/class_BacktestLVA.R")
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
  standardizeScores <- function(scores)
  {
    idx_notNA <- which( !is.na(scores) )
    x <- scores[idx_notNA]
    # Z-score
    x_stdz <- (x - mean(x)) / sd(x)
    # Winsorize
    x_stdz[ which(x_stdz < -3) ] <- -3
    x_stdz[ which(x_stdz > 3) ] <- 3
    scores_stdz <- scores
    scores_stdz[ idx_notNA ] <- x_stdz
    return( scores_stdz )
  }
  
  
  
  
  # --------------------------------------------------------------------------
  # PARAMETERS
  # --------------------------------------------------------------------------
  
  # universe <- "usa"
  # ccy <- "USD"
  universe <- "aktien_ch"
  ccy <- "CHF"
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  end_date <- "2022-03-06"
  
  
  
  # --------------------------------------------------------------------------
  # BENCHMARK SERIES
  # --------------------------------------------------------------------------
  
  tickers <- setNames( c("SPI"),
                       c("SPI") )
  X_bm <- rodbcGetOLZDBReturns( assetName = tickers,
                                refCcy = ccy,
                                frqncy = "daily",
                                na.rm = "r" )
  X_bm <- X_bm[isWeekday(time(X_bm)), tickers]
  colnames(X_bm) <- names(tickers)
  
  
  
  # --------------------------------------------------------------------------
  # INITIALIZE BACKTEST OBJECT
  # --------------------------------------------------------------------------
  
  BT0 <- loadBacktest( wd = paste0(wd_data, "Data/"),
                       name = paste0(universe, "_fub") )
  # BT0 <- loadBacktest( wd = paste0(wd_data, "Data/"),
  #                      name = universe )
  # rebdates <- BT0$spec$rebdates
  # BT0$spec$rebdates <- rebdates[ seq(from = 1, to = length(rebdates), by = 2) ]
  # rebdates <- BT0$spec$rebdates
  # BT0$spec$rebdates <- rebdates[ rebdates > "2011-01-01" ]
  BT0$initRBE()
  
  
  
  
  # --------------------------------------------------------------------------
  # SIZE SCORES
  # --------------------------------------------------------------------------
  
  # Size score
  Score_size <- BT0$data$wmat_bm
  Score_size_log <- log(Score_size)
  Score_size_stdz <- timeSeries( t( apply( Score_size_log, 1, FUN = standardizeScores ) ),
                                 time(Score_size) )
  
  lScore <- list( size = Score_size,
                  size_log = Score_size_log,
                  size_stdz = Score_size_stdz )
  lFactor <- lapply( lScore, FUN = function(x) { na.omit(x, method = "z" ) } )
  
  
  
  
  # --------------------------------------------------------------------------
  # FACTOR SCORES BM
  # --------------------------------------------------------------------------
  
  lScore_bm <- list()
  for ( factor_name in names(lScore) ) {
    lScore_bm[[ factor_name ]] <- portfolioScore( wmat = BT0$data$wmat_bm, 
                                                  score_mat = na.omit( lScore[[ factor_name ]], method = "z" ) )
  }
  factor_scores_bm <- do.call( cbind, lScore_bm )
  
  
  
  # # --------------------------------------------------------------------------
  # # MINIMUM VARIANCE PORTFOLIO
  # # --------------------------------------------------------------------------
  # 
  # BT_minV <- BT_Fctrl$copy()
  # BT_minV$spec$portfolio <- NULL
  # BT_minV$spec$cons_fun <- "box"
  # BT_minV$spec$width <- BT0$spec$width
  # BT_minV$spec$selection_filter <- c("upperbound", "NA")
  # BT_minV$runLoop()
  
  
  
  
  # --------------------------------------------------------------------------
  # LOW-SIZE QUANTILE PORTFOLIO
  # --------------------------------------------------------------------------
  
  BT_size <- BacktestPanelFactor$new()
  BT_size$data <- BT0$data
  BT_size$spec <- BT0$spec
  BT_size$data$upper_mat <- BT0$data$upper_mat * 0 + 1
  BT_size$spec$portfolio <- "quintile"
  BT_size$spec$factor_name <- "size_stdz"
  BT_size$data$factor_scores_bm <- factor_scores_bm
  BT_size$data$lFactor <- lFactor
  BT_size$spec$cons_fun <- "box"
  BT_size$spec$width <- BT0$spec$width
  BT_size$spec$selection_filter <- c("upperbound", "NA")
  BT_size$spec$weighting_method <- "eqw"
  BT_size$runLoop()
  
  
  # Simulate
  lSim_size <- lapply( BT_size$output, FUN = function(wmat) { simPortfolio( X = BT_size$data$X_sim,
                                                                            wghts = wmat,
                                                                            fc = 0, vc = 0 ) } )
  sim_size <- do.call( cbind, lSim_size )
  

  lScore_port <- list()
  for ( i in 1:length(BT_size$output) ) {
    lScore_port[[ i ]] <- portfolioScore( wmat = BT_size$output[[i]],
                                          score_mat = na.omit( lScore[[ "size_stdz" ]], method = "z" ) )
  }
  names(lScore_port) <- names(BT_size$output)
  factor_scores_port <- do.call( cbind, lScore_port )
  
  
  
  
    
  # --------------------------------------------------------------------------
  # MINIMUM TRACKING ERROR PORTFOLIO
  # --------------------------------------------------------------------------
  
  
  BT_minTE <- BacktestPanelFactor$new()
  BT_minTE$data <- BT_size$data
  BT_minTE$spec <- BT_size$spec
  BT_minTE$spec$cons_fun <- "longshort"
  BT_minTE$spec$selection_filter <- c("upperbound", "NA")
  BT_minTE$spec$portfolio <- "minTE"
  BT_minTE$data$wmat_bm <- BT_size$data$wmat_bm[rownames(BT_size$output$q0.2), ] * 0
  BT_minTE$data$wmat_bm[ ,colnames(BT_size$output$q0.2)] <- BT_size$output$q0.2
  BT_minTE$data$wmat_bm[ is.na(BT_minTE$data$wmat_bm) ] <- 0
  # debugonce( BT_minTE$minTEPortfolio )
  BT_minTE$runLoop()
  sim_minTE <- BT_minTE$simulate( fc = 0, vc = 0 )
  
  factor_scores_minTE <- portfolioScore( wmat = BT_minTE$output$weights,
                                         score_mat = na.omit( lScore[[ "size_stdz" ]], method = "z" ) )
  
  
  
  # --------------------------------------------------------------------------
  # MINIMUM TRACKING ERROR PORTFOLIO WITH LARGE SIZE CONTROL
  # --------------------------------------------------------------------------
  
  
  BT_minTE_sizectrl <- BT_minTE$copy()
  BT_minTE_sizectrl$output <- list()
  factor_scores_tmp <- factor_scores_port[ ,"q0.6"]
  colnames(factor_scores_tmp) <- "size_stdz"
  BT_minTE_sizectrl$data$factor_scores_bm <- factor_scores_tmp
  BT_minTE_sizectrl$spec$cons_fun <- c("box", "controlbmex10")
  BT_minTE_sizectrl$spec$factor_names <- "size_stdz"
  BT_minTE_sizectrl$spec$factor_names_constraints <- BT_minTE_sizectrl$spec$factor_names
  BT_minTE_sizectrl$runLoop()
  sim_minTE_sizectrl <- BT_minTE_sizectrl$simulate( fc = 0, vc = 0 )
  
  factor_scores_minTE_sizectrl <- portfolioScore( wmat = BT_minTE_sizectrl$output$weights,
                                                  score_mat = na.omit( lScore[[ "size_stdz" ]], method = "z" ) )
  
  
  
  # --------------------------------------------------------------------------
  # ANALYZE SIMULATIONS
  # --------------------------------------------------------------------------
  
  
  bm_names <- c("SPI")
  sim <- na.omit( cbind( # X_bm[ ,bm_names], 
                         size_quantile = sim_size, 
                         minTE = sim_minTE, 
                         minTE_sizectrl = sim_minTE_sizectrl ) )
  
  
  plotSimTS( as.simTS( sim ) )
  
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD", "te")
  t( descStats( sim, descStatsSpec( bmName = "size_quantile.q0.2") )$stats[stats_fields, ] )
  
  cor(sim)
  
  
  
  factor_scores <- na.omit( cbind( bm = factor_scores_bm[ ,3],
                                   factor_scores_port,
                                   minTE = factor_scores_minTE,
                                   minTE_sizectrl = factor_scores_minTE_sizectrl ) )  
  plot( factor_scores, plot.type = "single" )
  colors <- 1:ncol(factor_scores)
  legend( "topright", colnames(factor_scores), lwd = 2, col = colors, text.col = colors, bty = "n")
  
  
  
  
  weightsBarPlot( BT0$data$wmat_bm )
  weightsBarPlot( BT_size$output$q0.2 )
  weightsBarPlot( BT_minTE$output$weights )
  weightsBarPlot( BT_minTE_sizectrl$output$weights )

  
  today <- "2022-07-20"
  tmp <- BT0$data$X[which(BT0$data$X[ ,"Date"]== today), c("Series_Id", "Name")]
  wghts <- getWeights( rebalenv( rebdate = today )$GPO )
  Names <- names(wghts)[ names(wghts) %in% tmp[ ,"Series_Id"] ]
  idx <- rev( order(wghts[Names]) )
  data.frame( tmp[idx, ], wghts[idx] )
    
  
  
  
  
  
  
  