  
  
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
  universe <- "usa"
  source( "H:/R/github/Random_Portfolios/Source/class_BacktestCustom.R" )
  
  
  # Instantiate default Backtest class and prepare data for given universe
  BT0 <- BacktestCustom$new()
  
  # Default specifications for given universe
  BT0$setCtrl( universe = universe,
               name = paste0(universe, "_maxw5"),
               selection_filter = "db_flag",
               # selection_filter = "bm",
               # cons_fun = c("box", "country", "sector", "reit", "esgWhenAvailable"),
               cons_fun = "box",
               width = 52 * 3,
               wd_data = paste0(wd, "Data/"),
               wd_warehouse = paste0(wd, "waRehouse/"),
               clean_RBE = FALSE )
  
  # Default data for given universe
  BT0$inputData()
  BT0$data$upper_mat <- BT0$data$upper_mat * 0 + 0.05
  BT0$spec$rebdates <- BT0$spec$rebdates[ which(BT0$spec$rebdates < tail(rownames(BT0$data$X_sim), 1) ) ]
  
  # Save default backtest object
  BT0$save( wd = BT0$spec$wd_data,
            name = BT0$spec$name,
            without_data = FALSE )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  