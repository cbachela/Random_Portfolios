  
  
  ############################################################################
  ### GEOMSCALE LOW VOLA ANOMALY - PREPARE DATA
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     12.10.2021
  # First version:    12.10.2021
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # Require
  # --------------------------------------------------------------------------
  
  require(garcholz)
  require(simolz)
  require(BSS)
  require(covolz)
  require(Backtest)

  
  
  # --------------------------------------------------------------------------
  # Default parameters
  # --------------------------------------------------------------------------
  
  wd <- "R:/Asset_Management/Research_Projects/External_Research_Projects/GeomScale/Random_Portfolios/"
  wd_data <- paste0(wd, "Data/")
  wd_warehouse <- paste0(wd, "waRehouse/")
  # source( paste0(wd, "Source/custom_functions.R"), echo = TRUE )
  
  
  # --------------------------------------------------------------------------
  # Loov over universes
  # --------------------------------------------------------------------------
  
  universes <- c("usa", "dm", "europe_ex_ch", "europe", "eurozone", "aktien_ch_fub")
  
  for ( universe in universes ) {
    
    
    lSettings <- universe2lSettings(universe)
    width <- 260
    floatcon <- 1 / 4
    Covariance <- covCtrl(method = "qis")
    selection_filter <- c("db_flag")
    if ( universe == "aktien_ch_fub" ) {
      # cons_fun <- "box"
      ccy_sim <- "CHF"
    } else {
      # cons_fun <- c("box", "country", "sector", "reit", "esgWhenAvailable")
      ccy_sim <- "USD"
    }
    cons_fun <- "box"
    
    
    # --------------------------------------------------------------------------
    # Initialize Backtest object
    # --------------------------------------------------------------------------
    
    # Instantiate default Backtest class and prepare data for given universe
    BT0 <- BacktestPanel$new()
    
    # Default specifications for given universe
    BT0$setCtrl( universe = universe,
                 name = universe,
                 selection_filter = selection_filter,
                 cons_fun = cons_fun,
                 ccy_sim = ccy_sim,
                 width = width,
                 floatcon = floatcon,
                 wd = wd,
                 wd_warehouse = wd_warehouse,
                 clean_RBE = FALSE,
                 GPS = GPO::gps( Covariance = Covariance ) )
    
    
    # Default data for given universe
    BT0$inputData()
    
    # # Adaptive upper bounds
    # upper_mat_adaptive <- volume2BoxCon( Volume = BT0$data$Median,
    #                                      wmat_bm = BT0$data$wmat_bm,
    #                                      nav = 10^8 * 3,
    #                                      pct = 0.2,
    #                                      n_day = 5,
    #                                      bm_tol = 20,
    #                                      min_liq =  10^6 * 1.5,
    #                                      global_max = 0.03,
    #                                      global_min = 1e-04 )
    # BT0$data$upper_mat_adaptive <- upper_mat_adaptive
    
    
    # Save default backtest object
    BT0$save( wd = wd_data,
              name = BT0$spec$name,
              without_data = FALSE )
    
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  