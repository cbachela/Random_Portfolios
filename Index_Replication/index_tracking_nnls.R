
  
  ############################################################################
  ### INDEX TRACKING - NON-NEGATIVE LEAST SQUARES (NNLS)
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     06.03.2023
  # First version:    02.02.2023
  # --------------------------------------------------------------------------
  
  
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------

  # require(volesti)
  # require(stringr)
  # require(garcholz)
  # require(slolz)
  # require(simolz)
  # require(covolz)
  # require(BSS)
  # require(BBSolz)
  # require(Backtest)
  # require(RP)
  # wd_data <- "R:/Asset_Management/Research_Projects/External_Research_Projects/GeomScale/Random_Portfolios/"
  # wd <- "H:/R/github/Random_Portfolios/"
  # source( paste0(wd, "Source/class_BacktestCustom.R") )
  # source( paste0(wd, "Source/class_BacktestPanelFactor.R") )
  # source("H:/R/github/Random_Portfolios/Source/class_BacktestLVA.R")
  # source( paste0(wd, "Source/custom_functions.R") )

  
  require(garcholz)
  require(rugarch)
  require(DAARC)
  require(Backtest)
  
  wd_data <- "R:/Asset_Management/Research_Projects/External_Research_Projects/GeomScale/Random_Portfolios/"
  wd <- "H:/R/github/Random_Portfolios/"
  source( paste0(wd, "Source/class_BacktestCustom.R") )
  source( paste0(wd, "Source/class_BacktestIndexReplication.R") )
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # PARAMETERS
  # --------------------------------------------------------------------------
  
  universe <- "usa"
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  end_date <- "2022-03-06"
  fc <- 0
  vc <- 0
  
  

  # --------------------------------------------------------------------------
  # INITIALIZE BACKTEST OBJECT
  # --------------------------------------------------------------------------
  
  BT0 <- loadBacktest( wd = paste0(wd_data, "Data/"),
                       name = "usa_maxw5" )
  # rebdates <- BT0$spec$rebdates
  # BT0$spec$rebdates <- rebdates[ seq(from = 1, to = length(rebdates), by = 2) ]
  BT0$initRBE()
  

  rebdates <- BT0$spec$rebdates
  rebdates_w <- rownames(BT0$data$X_est)
  rebdates_w <- rebdates_w[ rebdates_w >= rebdates[1] ]
  
  
  
  
  # --------------------------------------------------------------------------
  # DEFINE INDEX SERIES TO BE REPLICATED
  # --------------------------------------------------------------------------

  X_bm <- BT0$data$X_bm
  
  
  SCNO <- DAARC:::Scenario$new()
  SCNO$data$X_bm <- X_bm
  # debugonce( SCNO$simUGagarchVola )
  # seed <- 7906
  seed <- 87654321
  prob <- SCNO$simUGagarchVola( seed = seed, alpha1 = 0.09, beta1 = 0.9  )
  SCNO$simulate( prob = prob)

  plotSimTS( as.simTS(SCNO$output$Z_eqw) )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # MINIMUM VARIANCE PORTFOLIO
  # --------------------------------------------------------------------------

  BT_minV <- BT0$copy()
  BT_minV$data$upper_mat <- BT0$data$upper_mat * 0 + 1
  BT_minV$spec$portfolio <- NULL
  BT_minV$spec$cons_fun <- "box"
  BT_minV$spec$selection_filter <- c("upperbound", "NA")
  BT_minV$runLoop()
  sim_minV <- BT_minV$simulate( fc = fc, vc = vc )
  
  
  
  # # --------------------------------------------------------------------------
  # # LOW-SIZE QUANTILE PORTFOLIO
  # # --------------------------------------------------------------------------
  # 
  # BT_size <- BT0$copy()
  # BT_size$spec$portfolio <- "quintile"
  # BT_size$spec$factor_name <- "size_stdz"
  # BT_size$spec$cons_fun <- "box"
  # BT_size$spec$width <- BT0$spec$width
  # BT_size$spec$selection_filter <- c("upperbound", "NA")
  # BT_size$spec$weighting_method <- "eqw"
  # # debugonce( BT_size$quantilePortfolio )
  # BT_size$runLoop()
  # 
  # # Simulate
  # lSim_size <- lapply( BT_size$output, FUN = function(wmat) { simPortfolio( X = BT_size$data$X_sim,
  #                                                                           wghts = wmat,
  #                                                                           fc = fc, vc = vc ) } )
  # sim_size <- do.call( cbind, lSim_size )

  
  
  
    
  # # --------------------------------------------------------------------------
  # # MINIMUM TRACKING ERROR PORTFOLIO TO LOW-SIZE PORTFOLIO
  # # --------------------------------------------------------------------------
  # 
  # BT_minTE <- BT0$copy()
  # BT_minTE$spec$cons_fun <- "box"
  # BT_minTE$spec$selection_filter <- c("upperbound", "NA")
  # BT_minTE$spec$portfolio <- "minTE"
  # BT_minTE$data$wmat_bm <- BT0$data$wmat_bm[rownames(BT_size$output$q0.2), ] * 0
  # BT_minTE$data$wmat_bm[ ,colnames(BT_size$output$q0.2)] <- BT_size$output$q0.2
  # BT_minTE$data$wmat_bm[ is.na(BT_minTE$data$wmat_bm) ] <- 0
  # # debugonce( BT_minTE$minTEPortfolio )
  # BT_minTE$runLoop()
  # sim_minTE <- BT_minTE$simulate( fc = fc, vc = vc )
  # 
  # 
  # 
  # # --------------------------------------------------------------------------
  # # MINIMUM TRACKING ERROR PORTFOLIO TO LOW-SIZE PORTFOLIO WITH LARGE SIZE CONTROL
  # # --------------------------------------------------------------------------
  # 
  # 
  # BT_minTE_sizectrl <- BT_minTE$copy()
  # BT_minTE_sizectrl$output <- list()
  # BT_minTE_sizectrl$spec$cons_fun <- c("box", "controlbmex10")
  # BT_minTE_sizectrl$spec$factor_names <- "size_stdz"
  # BT_minTE_sizectrl$spec$factor_names_constraints <- BT_minTE_sizectrl$spec$factor_names
  # BT_minTE_sizectrl$runLoop()
  # sim_minTE_sizectrl <- BT_minTE_sizectrl$simulate( fc = fc, vc = vc )
  
  
  
  
  # --------------------------------------------------------------------------
  # INDEX REPLICATING PORTFOLIO USING NNLS
  # --------------------------------------------------------------------------
  
  BT_nnls <- BacktestIndexReplication$new()
  BT_nnls$spec <- BT0$spec
  BT_nnls$data <- BT0$data
  BT_nnls$data$upper_mat <- BT0$data$upper_mat * 0 + 0.1
  BT_nnls$data$X_bm <- SCNO$output$Z_eqw
  # BT_nnls$spec$tau <- 52
  BT_nnls$spec$portfolio <- "nnls"
  # debugonce( BT_nnls$nnlsPortfolio )
  BT_nnls$runLoop()
  sim_nnls <- BT_nnls$simulate( fc = fc, vc = vc )
  
  
  
  
  # --------------------------------------------------------------------------
  # INDEX REPLICATING PORTFOLIO USING OLS
  # --------------------------------------------------------------------------
  
  BT_ols <- BacktestIndexReplication$new()
  BT_ols$spec <- BT0$spec
  BT_ols$data <- BT0$data
  BT_ols$data$upper_mat <- BT0$data$upper_mat * 0 + 0.1
  BT_ols$data$X_bm <- SCNO$output$Z_eqw
  BT_ols$spec$tau <- 52
  BT_ols$spec$portfolio <- "ols"
  # debugonce( BT_ols$olsPortfolio )
  BT_ols$runLoop()
  sim_ols <- BT_ols$simulate( fc = fc, vc = vc )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # INDEX REPLICATING PORTFOLIO USING NNLS - WEEKLY REBALANCING
  # --------------------------------------------------------------------------
  
  BT2 <- BacktestIndexReplication$new()
  BT2$data <- BT0$data
  BT2$spec <- BT0$spec
  dates_w <- rownames(BT0$data$X_est)
  BT2$spec$rebdates <- dates_w[dates_w >= BT0$spec$rebdates[1]]
  BT2$spec$portfolio <- "nnls"
  BT2$initRBE()
  BT2$runLoop()
  sim_bt2 <- BT2$simulate( fc = fc, vc = vc )
  
  
  
  # --------------------------------------------------------------------------
  # INDEX REPLICATING PORTFOLIO USING NNLS - WEEKLY REBALANCING, SHORT LOOKBACK
  # --------------------------------------------------------------------------
  
  BT3 <- BacktestIndexReplication$new()
  BT3$data <- BT0$data
  BT3$spec <- BT0$spec
  BT3$spec$width <- 2
  dates_w <- rownames(BT0$data$X_est)
  BT3$spec$rebdates <- dates_w[dates_w >= BT0$spec$rebdates[1]]
  BT3$spec$portfolio <- "nnls"
  BT3$initRBE()
  BT3$runLoop()
  sim_bt3 <- BT3$simulate( fc = fc, vc = vc )
  
  
  
  
  # --------------------------------------------------------------------------
  # INDEX REPLICATING PORTFOLIO USING ELSASTIC NET
  # --------------------------------------------------------------------------
  
  BT_elnet <- BacktestIndexReplication$new()
  BT_elnet$data <- BT0$data
  BT_elnet$spec <- BT0$spec
  BT_elnet$spec$rebdates <- rebdates_w
  BT_elnet$data$X_bm <- SCNO$output$Z_eqw
  BT_elnet$spec$tau <- 21 * 1
  BT_elnet$spec$elnet_alpha <- 1  # 1 == lasso (l1-penalty, induces sparsity), 0 = ridge (l2-penalty)
  BT_elnet$spec$portfolio <- "elnet"
  BT_elnet$initRBE()
  # debugonce( BT_elnet$elnetPortfolio )
  BT_elnet$runLoop()
  sim_elnet <- BT_elnet$simulate( fc = fc, vc = vc )
  
  
  
  BT_ridge <- BT_elnet$copy()
  BT_ridge$output <- list()
  BT_ridge$spec$tau <- NULL
  BT_ridge$spec$elnet_alpha <- 0  # 1 == lasso (l1-penalty, induces sparsity), 0 = ridge (l2-penalty)
  BT_ridge$spec$portfolio <- "elnet"
  BT_ridge$runLoop()
  sim_ridge <- BT_ridge$simulate( fc = fc, vc = vc )
  
  
  BT_lasso <- BT_elnet$copy()
  BT_lasso$output <- list()
  BT_lasso$spec$tau <- NULL
  BT_lasso$spec$elnet_alpha <- 1  # 1 == lasso (l1-penalty, induces sparsity), 0 = ridge (l2-penalty)
  BT_lasso$spec$portfolio <- "elnet"
  BT_lasso$runLoop()
  sim_lasso <- BT_lasso$simulate( fc = fc, vc = vc )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # WEIGHTS
  # --------------------------------------------------------------------------
  
  
  
  sector_mat_bm <- groupWeights( wmat = BT0$data$wmat_bm,
                                  group_mat = BT0$data$sector_mat )
  sector_mat_bt2 <- groupWeights( wmat = BT2$output$weights,
                                  group_mat = BT0$data$sector_mat )
  sector_mat_bt3 <- groupWeights( wmat = BT3$output$weights,
                                  group_mat = BT0$data$sector_mat )
  
  weightsBarPlot( sector_mat_bt2[BT0$spec$rebdates, ] )
  weightsBarPlot( sector_mat_bt3[BT0$spec$rebdates, ] )
  weightsBarPlot( sector_mat_bm[BT0$spec$rebdates, ] )
  
  weightsBarPlot( BT_nnls$output$weights )
  weightsBarPlot( BT_ols$output$weights )
  weightsBarPlot( BT_elnet$output$weights )
  
  
  
  
  # --------------------------------------------------------------------------
  # ANALYZE SIMULATIONS
  # --------------------------------------------------------------------------
  
  
  sim <- na.omit( cbind( bm = X_bm,
                         bm_synt = SCNO$output$Z_eqw,
                         minvar = sim_minV,
                         nnls = sim_nnls,
                         ols = sim_ols,
                         elnet = sim_elnet,
                         lasso = sim_lasso,
                         ridge = sim_ridge ) )
  
  sim <- na.omit( cbind( bm_synt = SCNO$output$Z_eqw,
                         elnet = sim_elnet ) )
  
  
  plotSimTS( as.simTS( sim ) )
  # plotSimTS( as.simTS( sim[rownames(sim) > "2003-01-01", ] ) )
  
  
  stats_fields <- c("cumret", "means", "sds", "dsds", "sharpe", "maxDD", "te", "ir")
  t( descStats( sim, descStatsSpec( bmName = "bm_synt") )$stats[stats_fields, ] )
  
  
  # Relative performance
  sim_rel <- timeSeries( apply( sim, 2, simOutperformance, y = as.numeric(sim[ ,"bm"]) ),
                         time(sim) )
  plotSimTS( as.simTS(sim_rel) )      
  
  stats_fields <- c("cumret", "means", "sds", "dsds", "sharpe", "maxDD")
  t( descStats( sim_rel, descStatsSpec( bmName = "bm") )$stats[stats_fields, ] )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # SAVE
  # --------------------------------------------------------------------------
  
  env <- new.env()
  env$sim_minV <- sim_minV
  env$sim_bt <- sim_bt
  
  saveRDS( object = env, file = paste0(wd_data, "waRehouse/index_tracking_nnls.rds") )
  
  
  
  
  
  
  
  
  
  
  
  