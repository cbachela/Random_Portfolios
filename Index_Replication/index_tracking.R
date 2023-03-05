
  
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
  # INITIALIZE BACKTEST OBJECT
  # --------------------------------------------------------------------------
  
  BT0 <- loadBacktest( wd = paste0(wd_data, "Data/"),
                       name = "usa_maxw5" )
  rebdates <- BT0$spec$rebdates
  BT0$spec$rebdates <- rebdates[ seq(from = 1, to = length(rebdates), by = 2) ]
  BT0$initRBE()
  
 
  # Load Multifactor backtest object
  BT_Fctrl <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "multifac_fctrl" )
  BT_Fctrl$spec$GPS@covariance <- covCtrl( method = "ewma", 
                                           ellipsis = list(tau = 21) )
  BT_Fctrl$data$X_est <- BT_Fctrl$data$X_est_d[isWeekday(time(BT_Fctrl$data$X_est_d)), ]
  BT_Fctrl$spec$width <- 252 * 1
  
  
  
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
  
  BT_size <- BT_Fctrl$copy()
  BT_size$spec$portfolio <- "quintile"
  BT_size$spec$factor_name <- "size_stdz"
  BT_size$spec$cons_fun <- "box"
  BT_size$spec$width <- BT0$spec$width
  BT_size$spec$selection_filter <- c("upperbound", "NA")
  BT_size$spec$weighting_method <- "eqw"
  # debugonce( BT_size$quantilePortfolio )
  BT_size$runLoop()
  
  # Simulate
  lSim_size <- lapply( BT_size$output, FUN = function(wmat) { simPortfolio( X = BT_size$data$X_sim,
                                                                            wghts = wmat,
                                                                            fc = 0, vc = 0 ) } )
  sim_size <- do.call( cbind, lSim_size )
  
  
  names(BT_size$output)
  weightsBarPlot( BT_size$output$q0.2 )
  
  
  
    
  # --------------------------------------------------------------------------
  # MINIMUM TRACKING ERROR PORTFOLIO
  # --------------------------------------------------------------------------
  
  
  BT_minTE <- BT_Fctrl$copy()
  BT_minTE$spec$cons_fun <- "box"
  BT_minTE$spec$selection_filter <- c("upperbound", "NA")
  BT_minTE$spec$portfolio <- "minTE"
  BT_minTE$data$wmat_bm <- BT_Fctrl$data$wmat_bm[rownames(BT_size$output$q0.2), ] * 0
  BT_minTE$data$wmat_bm[ ,colnames(BT_size$output$q0.2)] <- BT_size$output$q0.2
  BT_minTE$data$wmat_bm[ is.na(BT_minTE$data$wmat_bm) ] <- 0
  # debugonce( BT_minTE$minTEPortfolio )
  BT_minTE$runLoop()
  sim_minTE <- BT_minTE$simulate( fc = 0, vc = 0 )
  
  
  
  # --------------------------------------------------------------------------
  # MINIMUM TRACKING ERROR PORTFOLIO WITH LARGE SIZE CONTROL
  # --------------------------------------------------------------------------
  
  
  BT_minTE_sizectrl <- BT_minTE$copy()
  BT_minTE_sizectrl$output <- list()
  BT_minTE_sizectrl$spec$cons_fun <- c("box", "controlbmex10")
  BT_minTE_sizectrl$spec$factor_names <- "size_stdz"
  BT_minTE_sizectrl$spec$factor_names_constraints <- BT_minTE_sizectrl$spec$factor_names
  BT_minTE_sizectrl$runLoop()
  sim_minTE_sizectrl <- BT_minTE_sizectrl$simulate( fc = 0, vc = 0 )
  
  
  
  
  # --------------------------------------------------------------------------
  # ANALYZE SIMULATIONS
  # --------------------------------------------------------------------------
  
  
  bm_names <- c("Capw", "Size")
  sim <- na.omit( cbind( X_fact[ ,bm_names], 
                         size_quantile = sim_size, 
                         minTE = sim_minTE, 
                         minTE_sizectrl = sim_minTE_sizectrl ) )
  
  
  plotSimTS( as.simTS( sim ) )
  
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD", "te")
  t( descStats( sim, descStatsSpec( bmName = "size_quantile.q0.2") )$stats[stats_fields, ] )
  t( descStats( sim, descStatsSpec( bmName = "Size") )$stats[stats_fields, ] )
  t( descStats( sim, descStatsSpec( bmName = "Capw") )$stats[stats_fields, ] )
  
  
  
  
  
  
  
  
  
  
  
  