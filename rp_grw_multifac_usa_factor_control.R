
  
  ############################################################################
  ### RANDOM PORTFOLIOS - GEOMETRIC RANDOM WALK SAMPLING - FACTOR CONTROL - USA
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
  source( paste0(wd, "Source/class_BacktestPanelFactor.R") )
  source( paste0(wd, "Source/custom_functions.R") )

  
  
  
  
  # --------------------------------------------------------------------------
  # FUNCTIONS
  # --------------------------------------------------------------------------
  
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
  BT_qcapw$spec$name <- "qcapw"
  BT_qcapw$spec$portfolio <- "qcapw"
  BT_qcapw$spec$GPS@covariance$method <- "duv"
  BT_qcapw$runLoop()
  sim_qcapw <- BT_qcapw$simulate( fc = 0, vc = 0 )
  sim_qcapw <- sim_qcapw[isWeekday(time(sim_qcapw)), ]
  BT_qcapw$output$simulations <- sim_qcapw
  
  BT_qcapw$save( wd = paste0(wd_data, "waRehouse/"), 
                 name = BT_qcapw$spec$name,
                 without_data = TRUE )
  
  
  
  
  # --------------------------------------------------------------------------
  # BACKTEST QUASI EQW PORTFOLIO
  # --------------------------------------------------------------------------
  
  BT_qeqw <- BT_MF$copy()
  BT_qeqw$spec$name <- "qeqw"
  BT_qeqw$spec$portfolio <- NULL
  BT_qeqw$spec$GPS@covariance$method <- "duv"
  BT_qeqw$runLoop()
  sim_qeqw <- BT_qeqw$simulate( fc = 0, vc = 0 )
  sim_qeqw <- sim_qeqw[isWeekday(time(sim_qeqw)), ]
  BT_qeqw$output$simulations <- sim_qeqw
  
  BT_qeqw$save( wd = paste0(wd_data, "waRehouse/"),
                name = BT_qeqw$spec$name,
                without_data = TRUE )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # PRE-COMPUTE COMPANY FACTOR SCORES
  # --------------------------------------------------------------------------
  
  # Size score
  Score_size <- BT0$data$wmat_bm
  Score_size_log <- log(Score_size)
  Score_size_stdz <- timeSeries( t( apply( Score_size_log, 1, FUN = standardizeScores ) ),
                                 time(Score_size) )
  
  # Quality
  Score_qual <- BT0$data$quality_mat
  
  # Value 
  Score_val <- BT0$data$value_mat 
  
  
  # Compute Momentum scores
  BT_mom <- BT_MF$copy()
  BT_mom$spec$width <- 52
  BT_mom$spec$portfolio <- "momScore"
  BT_mom$runLoop()
  Score_mom <- BT_mom$output$scores
  Score_mom_stdz <- timeSeries( t( apply( Score_mom, 1, FUN = standardizeScores ) ),
                                 time(Score_mom) )
  
  # Compute Vola scores
  BT_vol <- BT_MF$copy()
  BT_vol$spec$width <- 52
  BT_vol$spec$portfolio <- "volScore"
  BT_vol$runLoop()
  Score_vol <- BT_vol$output$scores
  Score_vol_stdz <- timeSeries( t( apply( Score_vol, 1, FUN = standardizeScores ) ),
                                time(Score_vol) )
  
  
  lScore <- list( momentum = Score_mom,
                  momentum_stdz = Score_mom_stdz,
                  value = Score_val,
                  quality = Score_qual,
                  size = Score_size,
                  size_log = Score_size_log,
                  size_stdz = Score_size_stdz,
                  lowvola = Score_vol,
                  lowvola_stdz = Score_vol_stdz )
  lFactor <- lapply( lScore, FUN = function(x) { na.omit(x, method = "z" ) } )
  
  
  # --------------------------------------------------------------------------
  # FACTOR SCORES BM
  # --------------------------------------------------------------------------
  
  lScore_bm <- list()
  for ( factor_name in names(lScore) ) {
    lScore_bm[[ factor_name ]] <- portfolioScore( wmat = BT_MF$data$wmat_bm, 
                                                  score_mat = na.omit( lScore[[ factor_name ]], method = "z" ) )
  }
  factor_scores_bm <- do.call( cbind, lScore_bm )
  
  
  
  # --------------------------------------------------------------------------
  # BACKTEST FACTOR CONTROL
  # --------------------------------------------------------------------------
  
  BT_Fctrl <- BacktestPanelFactor$new()
  BT_Fctrl$spec <- BT_MF$spec
  BT_Fctrl$data <- BT_MF$data
  BT_Fctrl$data$factor_scores_bm <- factor_scores_bm
  BT_Fctrl$data$lFactor <- lFactor
  BT_Fctrl$spec$cons_fun <- c(BT_MF$spec$cons_fun, "controlbmex10")
  BT_Fctrl$spec$name <- "multifac_fctrl"

  BT_Fctrl$save( wd = paste0(wd_data, "waRehouse/"), 
                 name = BT_Fctrl$spec$name,
                 without_data = FALSE )
  
  # BT_Fctrl <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "multifac_fctrl" )
  
  
  
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
  BT2$output$simulations <- sim_2
  
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
  tmp <- lapply( BT2$output, FUN = simFUN )
  sim_2_cons <- do.call( cbind, tmp )
  sim_2_cons <- sim_2_cons[isWeekday(time(sim_2_cons)), ]
  colnames(sim_2_cons) <- paste0("rp_grw_uniform_cons_valctrl", 1:ncol(sim_2_cons))
  BT2_fev2cons$output$simulations <- sim_2_cons
  
  # Save
  BT2$save( wd = paste0(wd_data, "waRehouse/"), name = BT2$spec$name, without_data = TRUE )
  BT2_fev2cons$save( wd = paste0(wd_data, "waRehouse/"), name = BT2_fev2cons$spec$name, without_data = TRUE )
  
  
  

  
  # --------------------------------------------------------------------------
  # BACKTEST UNIFORM RP'S WITHIN CONSTRAINTS WITH QUALITY CONTROL
  # --------------------------------------------------------------------------
  
  BT3 <- BT_Fctrl$copy()
  BT3$spec$factor_names <- "quality"
  BT3$spec$factor_names_constraints <- BT3$spec$factor_names
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
  BT3$output$simulations <- sim_3
  
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
  BT3_fev2cons$output$simulations <- sim_3_cons
  
  # Save
  BT3$save( wd = paste0(wd_data, "waRehouse/"), name = BT3$spec$name, without_data = TRUE )
  BT3_fev2cons$save( wd = paste0(wd_data, "waRehouse/"), name = BT3_fev2cons$spec$name, without_data = TRUE )
  
 
  
  
  
  
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
  BT4$output$simulations <- sim_4
  
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
  BT4_fev2cons$output$simulations <- sim_4_cons
  
  # Save
  BT4$save( wd = paste0(wd_data, "waRehouse/"), name = BT4$spec$name, without_data = TRUE )
  BT4_fev2cons$save( wd = paste0(wd_data, "waRehouse/"), name = BT4_fev2cons$spec$name, without_data = TRUE )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BACKTEST UNIFORM RP'S WITHIN CONSTRAINTS WITH SIZE STANDARDIZED CONTROL
  # --------------------------------------------------------------------------
  
  BT5 <- BT_Fctrl$copy()
  BT5$spec$factor_names <- "size_stdz"
  BT5$spec$factor_names_constraints <- BT5$spec$factor_names
  BT5$spec$keep_all_rebdates <- FALSE
  BT5$spec$portfolio <- "grw"
  BT5$spec$n_sim <- 10^2 * 5
  BT5$spec$name <- "grw_multifac_usa_sizestdzctrl"
  BT5$output <- list()
  BT5$runLoop()
  
  
  # today <- "2002-03-06"
  # RBE <- rebalenv( rebdate = today )
  # Cons <- getConstraints(RBE$GPS)
  # Cons@linear$Amat[ ,1:5]
  # Cons@linear$rhs
  # portfolioScore( wmat = BT5$output[[4]], score_mat = BT5$data$lFactor[["size_stdz"]] )
  # 
  # 
  # selection <- getConstraints( Cons, "selection" )
  # wghts <- BT5$output[[100]][today, selection]
  # wghts <- setNames( as.numeric(wghts), names(wghts) )
  # w_fev <- fevBias( x = wghts, q = 6 )
  # 
  # 
  # # Remove some linear constraints for simplicity
  # Cons_alt <- Cons
  # idx <- c(1, 23, 24)
  # Cons_alt@linear$Amat <- Cons@linear$Amat[idx, ]
  # Cons_alt@linear$rhs <- Cons@linear$rhs[idx]
  # Cons_alt@linear$sense <- Cons@linear$sense[idx]
  # 
  # 
  # # debugonce( map2Constraints )
  # w_fev2con <- map2Constraints( w_unc = w_fev, 
  #                               # Constraints = Cons, 
  #                               Constraints = Cons_alt,
  #                               rm_superfluous = FALSE )
  # 
  # # debugonce( checkWeights )
  # checkWeights( wghts = wghts, spec = Cons )
  # checkWeights( wghts = w_fev, spec = Cons )
  # checkWeights( wghts = w_fev2con, spec = Cons, eps = 1e-06 )
  # 
  # cbind( wghts, w_fev, w_fev2con )
  # 
  # Cons_alt@linear$Amat %*% cbind(wghts, w_fev, w_fev2con)
  
  
  
  ###
  BT5 <- loadBacktest( wd = paste0(wd_data, "waRehouse/"), name = "grw_multifac_usa_sizestdzctrl" )
  BT5$data <- BT_Fctrl$data
  ###
  
  # Keep constraints
  lapply( BT5$spec$rebdates, FUN = BT5$gps )
  lCons <- lapply( BT5$spec$rebdates, FUN = function(today) { getConstraints(rebalenv(today)$GPS) } )
  names(lCons) <- BT5$spec$rebdates
  
  # Simulate
  tmp <- lapply( BT5$output[ grepl("grw", names(BT5$output)) ], FUN = simFUN )
  sim_5 <- do.call( cbind, tmp )
  sim_5 <- sim_5[isWeekday(time(sim_5)), ]
  colnames(sim_5) <- paste0("rp_grw_uniform_sizestdzctrl", 1:ncol(sim_5))
  BT5$output$simulations <- sim_5
  
  # Apply FEV-bias
  lWeights_ctrl <- BT5$output[ grepl("grw", names(BT5$output)) ]
  lWeights_fev_ctrl <- lapply( lWeights_ctrl, FUN = function(wmat) { fevBias( x = wmat, q = 6 ) } )
  lWeights_fev2_ctrl <- lapply( lWeights_ctrl, FUN = function(wmat) { fevBias( x = wmat, q = 2 ) } )
 
  # Project FEV-biased weights to boundary of feasible set
  lWeights_cons_ctrl <- lWeights_ctrl
  for ( i in 1:length(lWeights_cons_ctrl) ) {
    
    for ( today in rownames(lWeights_cons_ctrl[[i]]) ) {
      
      # RBE <- rebalenv( rebdate = today )
      w_unc <- lWeights_fev_ctrl[[i]][today, ]
      # idx_notNA <- which(!is.na(w_unc))
      idx_notNA <- which( !is.na( lWeights_ctrl[[i]][today, ] ) )
      w_unc <- w_unc[idx_notNA]
      w_unc <- setNames( as.numeric(w_unc), names(w_unc))
      # boxcon <- getConstraints( RBE$GPS, "bounds" )
      # upper_bounds <- boxcon$upper
      # # debugonce( map2boxcon )
      # lWeights_cons_ctrl[[i]][today, idx_notNA] <- map2boxcon( w_unc = w_unc,
      #                                                          upper = upper_bounds,
      #                                                          itermax = 10^3 )
      # Cons <- getConstraints( RBE$GPS )
      Cons <- lCons[[today]]
      selection <- getConstraints( Cons, "selection" )
      w_unc <- w_unc[ selection ]
      # debugonce( map2Constraints )
      w_con <- map2Constraints( w_unc = w_unc, 
                                Constraints = Cons, 
                                rm_superfluous = FALSE )
      lWeights_cons_ctrl[[i]][today, idx_notNA] <- w_con
    }
  }
  
  
  # --------------------------------------------------------------------------
  fev2cons <- function( BT )
  {
    # Apply FEV-bias
    lWeights <- BT$output[ grepl("grw", names(BT$output)) ]
    lWeights_fev <- lapply( lWeights, FUN = function(wmat) { fevBias( x = wmat, q = 6 ) } )
   
    # Project FEV-biased weights to boundary of feasible set
    lWeights_fev2cons <- lWeights
    for ( i in 1:length(lWeights_fev2cons) ) {
      for ( today in rownames(lWeights_fev2cons[[i]]) ) {
        
        RBE <- rebalenv( rebdate = today )
        w_unc <- lWeights_fev[[i]][today, ]
        # idx_notNA <- which(!is.na(w_unc))
        idx_notNA <- which( !is.na( lWeights[[i]][today, ] ) )
        w_unc <- w_unc[idx_notNA]
        w_unc <- setNames( as.numeric(w_unc), names(w_unc))
        # boxcon <- getConstraints( RBE$GPS, "bounds" )
        # upper_bounds <- boxcon$upper
        # # debugonce( map2boxcon )
        # lWeights_fev2cons[[i]][today, idx_notNA] <- map2boxcon( w_unc = w_unc,
        #                                                         upper = upper_bounds,
        #                                                         itermax = 10^3 )
        Cons <- getConstraints( RBE$GPS )
        # debugonce( map2Constraints )
        w_con <- map2Constraints( w_unc = w_unc, 
                                  Constraints = Cons, 
                                  rm_superfluous = FALSE )
        lWeights_fev2cons[[i]][today, idx_notNA] <- w_con
        
      }
    }
    return( lWeights_fev2cons )
  }
  
  
  lWeights_fev2cons <- fev2cons( BT = BT5 )
  
  tmp <- portfolioScore( wmat = lWeights_fev2cons[[1]], 
                         score_mat = BT5$data$lFactor[["size_stdz"]] )
  head(tmp)
  
  tmp2 <- portfolioScore( wmat = BT5$output[[2]], 
                          score_mat = BT5$data$lFactor[["size_stdz"]] )
  head(tmp2)
  

  fac_name <- "size_stdz"
  # lFacScore <- portfolioScore( wmat = lWeights_fev2cons,
  #                              score_mat = BT_Fctrl$data$lFactor[[fac_name]] )
  # lFacScore <- portfolioScore( wmat = lWeights_ctrl,
  #                              score_mat = BT_Fctrl$data$lFactor[[fac_name]] )
  lFacScore <- portfolioScore( wmat = lWeights_fev2_ctrl,
                               score_mat = BT_Fctrl$data$lFactor[[fac_name]] )
  fac_scores <- do.call( cbind, lFacScore )
  tmp <- cbind( BT_Fctrl$data$factor_scores_bm[ ,fac_name], 
                BT_Fctrl$data$factor_scores_bm[ ,fac_name] - 0.1,
                BT_Fctrl$data$factor_scores_bm[ ,fac_name] + 0.1 )
  dates <- intersect( rownames(fac_scores), rownames(tmp) )
  plot( fac_scores[dates, ], plot.type = "single", col = "grey" )
  lines( tmp[dates ,1], lwd = 2 )
  lines( tmp[dates ,2], col = 2, lwd = 2 )
  lines( tmp[dates ,3], col = 2, lwd = 2 )
  
  ####################
  
  
  lWeights_fev2cons <-  BT5$output[ grepl("grw", names(BT5$output)) ]
  for ( i in 1:length(lWeights_fev2cons) ) {
    for ( today in rownames(lWeights_fev2cons[[1]]) ) {
  
      Cons <- lCons[[today]]
      selection <- getConstraints( Cons, "selection" )
      w_unc <- rdirichlet( n = 1, alpha = rep(1, length(selection)) )
      w_unc <- setNames( as.numeric(w_unc), selection )
      w_unc_fev <- fevBias( x = w_unc, q = 6 )
      w_con <- map2Constraints( w_unc = w_unc_fev, 
                                Constraints = Cons, 
                                rm_superfluous = FALSE )
      lWeights_fev2cons[[i]][today, selection] <- w_con
    }
  }
  
  fac_name <- "size_stdz"
  lFacScore <- portfolioScore( wmat = lWeights_fev2cons,
                               score_mat = BT_Fctrl$data$lFactor[[fac_name]] )
  fac_scores <- do.call( cbind, lFacScore )
  tmp <- cbind( BT_Fctrl$data$factor_scores_bm[ ,fac_name], 
                BT_Fctrl$data$factor_scores_bm[ ,fac_name] - 0.1,
                BT_Fctrl$data$factor_scores_bm[ ,fac_name] + 0.1 )
  dates <- intersect( rownames(fac_scores), rownames(tmp) )
  plot( fac_scores[dates, ], plot.type = "single", col = "grey" )
  lines( tmp[dates ,1], lwd = 2 )
  lines( tmp[dates ,2], col = 2, lwd = 2 )
  lines( tmp[dates ,3], col = 2, lwd = 2 )
  
  
  
  
  
  
  BT5_fev2cons <- BT5$copy()
  BT5_fev2cons$spec$name <- paste0(BT5$spec$name, "_fev2cons")
  BT5_fev2cons$output <- fev2cons( BT = BT5 )
  
  # Simulate
  tmp <- lapply( BT5_fev2cons$output, FUN = simFUN )
  sim_5_cons <- do.call( cbind, tmp )
  sim_5_cons <- sim_5_cons[isWeekday(time(sim_5_cons)), ]
  colnames(sim_5_cons) <- paste0("rp_grw_uniform_cons_sizestdzctrl", 1:ncol(sim_5_cons))
  BT5_fev2cons$output$simulations <- sim_5_cons
  
  # Save
  BT5$save( wd = paste0(wd_data, "waRehouse/"), name = BT5$spec$name, without_data = TRUE )
  BT5_fev2cons$save( wd = paste0(wd_data, "waRehouse/"), name = BT5_fev2cons$spec$name, without_data = TRUE )
  
  
  
  
  # --------------------------------------------------------------------------
  # BACKTEST UNIFORM RP'S WITHIN CONSTRAINTS WITH VALUE, QUALITY, SIZE STANDARDIZED CONTROL
  # --------------------------------------------------------------------------
  
  BT6 <- BT_Fctrl$copy()
  BT6$spec$factor_names <- c("value", "quality", "size_stdz")
  BT6$spec$factor_names_constraints <- BT6$spec$factor_names
  BT6$spec$keep_all_rebdates <- FALSE
  BT6$spec$portfolio <- "grw"
  BT6$spec$n_sim <- 10^2 * 5
  BT6$spec$name <- "grw_multifac_usa_valqualsizectrl"
  BT6$output <- list()
  BT6$runLoop()
  
  # RBE <- rebalenv( rebdate = "2002-03-06" )
  # Cons <- getConstraints( RBE$GPS, "linear" )
  # rownames(Cons$Amat)
  
  # Simulate
  tmp <- lapply( BT6$output[ grepl("grw", names(BT6$output)) ], FUN = simFUN )
  sim_6 <- do.call( cbind, tmp )
  sim_6 <- sim_6[isWeekday(time(sim_6)), ]
  colnames(sim_6) <- paste0("rp_grw_uniform_valqualsizectrl", 1:ncol(sim_6))
  BT6$output$simulations <- sim_6
  
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
  BT6_fev2cons$output$simulations <- sim_6_cons
  
  # Save
  BT6$save( wd = paste0(wd_data, "waRehouse/"), 
            name = BT6$spec$name, 
            without_data = TRUE )
  BT6_fev2cons$save( wd = paste0(wd_data, "waRehouse/"), 
                     name = BT6_fev2cons$spec$name, 
                     without_data = TRUE )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BACKTEST UNIFORM RP'S WITHIN CONSTRAINTS WITH MOMENTUM CONTROL
  # --------------------------------------------------------------------------
  
  BT7 <- BT_Fctrl$copy()
  BT7$spec$factor_names <- c("momentum_stdz")
  BT7$spec$factor_names_constraints <- BT7$spec$factor_names
  BT7$spec$keep_all_rebdates <- FALSE
  BT7$spec$portfolio <- "grw"
  BT7$spec$n_sim <- 10^2 * 5
  BT7$spec$name <- "grw_multifac_usa_momctrl"
  BT7$output <- list()
  BT7$runLoop()
  
  # Simulate
  tmp <- lapply( BT7$output[ grepl("grw", names(BT7$output)) ], FUN = simFUN )
  sim_7 <- do.call( cbind, tmp )
  sim_7 <- sim_7[isWeekday(time(sim_7)), ]
  colnames(sim_7) <- paste0("rp_grw_uniform_momctrl", 1:ncol(sim_7))
  BT7$output$simulations <- sim_7
  
  # Apply FEV-bias
  lWeights_ctrl <- BT7$output[ grepl("grw", names(BT7$output)) ]
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
  BT7_fev2cons <- BT7$copy()
  BT7_fev2cons$spec$name <- paste0(BT7$spec$name, "_fev2cons")
  BT7_fev2cons$output <- lWeights_cons_ctrl
  
  # Simulate
  tmp <- lapply( BT7_fev2cons$output, FUN = simFUN )
  sim_7_cons <- do.call( cbind, tmp )
  sim_7_cons <- sim_7_cons[isWeekday(time(sim_7_cons)), ]
  colnames(sim_7_cons) <- paste0("rp_grw_uniform_cons_momctrl", 1:ncol(sim_7_cons))
  BT7_fev2cons$output$simulations <- sim_7_cons
  
  # Save
  BT7$save( wd = paste0(wd_data, "waRehouse/"), 
            name = BT7$spec$name, 
            without_data = TRUE )
  BT7_fev2cons$save( wd = paste0(wd_data, "waRehouse/"), 
                     name = BT7_fev2cons$spec$name, 
                     without_data = TRUE )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BACKTEST UNIFORM RP'S WITHIN CONSTRAINTS WITH 
  # VALUE, QUALITY, SIZE STANDARDIZED, MOMENTUM CONTROL
  # --------------------------------------------------------------------------
  
  BT8 <- BT_Fctrl$copy()
  BT8$spec$factor_names <- c("value", "quality", "size_stdz", "momentum_stdz")
  BT8$spec$factor_names_constraints <- BT8$spec$factor_names
  BT8$spec$keep_all_rebdates <- FALSE
  BT8$spec$portfolio <- "grw"
  BT8$spec$n_sim <- 10^2 * 5
  BT8$spec$name <- "grw_multifac_usa_valqualsizemomctrl"
  BT8$output <- list()
  BT8$runLoop()
  
  # RBE <- rebalenv( rebdate = "2002-03-06" )
  # Cons <- getConstraints( RBE$GPS, "linear" )
  # rownames(Cons$Amat)
  
  # Simulate
  tmp <- lapply( BT8$output[ grepl("grw", names(BT8$output)) ], FUN = simFUN )
  sim_8 <- do.call( cbind, tmp )
  sim_8 <- sim_8[isWeekday(time(sim_8)), ]
  colnames(sim_8) <- paste0("rp_grw_uniform_valqualsizemomctrl", 1:ncol(sim_8))
  BT8$output$simulations <- sim_8
  
  # Apply FEV-bias
  lWeights_ctrl <- BT8$output[ grepl("grw", names(BT8$output)) ]
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
  BT8_fev2cons <- BT8$copy()
  BT8_fev2cons$spec$name <- paste0(BT8$spec$name, "_fev2cons")
  BT8_fev2cons$output <- lWeights_cons_ctrl
  
  # Simulate
  tmp <- lapply( BT8_fev2cons$output, FUN = simFUN )
  sim_8_cons <- do.call( cbind, tmp )
  sim_8_cons <- sim_8_cons[isWeekday(time(sim_8_cons)), ]
  colnames(sim_8_cons) <- paste0("rp_grw_uniform_cons_valqualsizemomctrl", 1:ncol(sim_8_cons))
  BT8_fev2cons$output$simulations <- sim_8_cons
  
  # Save
  BT8$save( wd = paste0(wd_data, "waRehouse/"), 
            name = BT8$spec$name, 
            without_data = TRUE )
  BT8_fev2cons$save( wd = paste0(wd_data, "waRehouse/"), 
                     name = BT8_fev2cons$spec$name, 
                     without_data = TRUE )
  
  