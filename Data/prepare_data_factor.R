  
  
  ############################################################################
  ### RANDOM PORTFOLIO - PREPARE DATA - FACTOR MATRICES
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     06.01.2023
  # First version:    06.01.2023
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # Require
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
  wd_data <- "R:/Asset_Management/Research_Projects/External_Research_Projects/GeomScale/Random_Portfolios/"
  wd <- "H:/R/github/Random_Portfolios/"
  source( paste0(wd, "Source/class_BacktestCustom.R") )
  source( paste0(wd, "Source/custom_functions.R") )

  
  
  # # --------------------------------------------------------------------------
  # # XXX
  # # --------------------------------------------------------------------------
  # 
  # universe <- "usa"
  # BT0 <- loadBacktest( wd = paste0(wd_data, "Data/"),
  #                      name = "usa_maxw5" )
  # 
  # 
  # # Compute Momentum scores
  # BT_mom <- BacktestCustom$new()
  # BT_mom$data <- BT0$data
  # BT_mom$spec <- BT0$spec
  # BT_mom$spec$width <- 13
  # BT_mom$spec$portfolio <- "momScore"
  # BT_mom$runLoop()
  # Score_mom <- BT_mom$output$scores
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Load factor matrices from factor Shiny tool
  # --------------------------------------------------------------------------

  wd <- "R:/Asset_Management/R/Shiny/Factor_Portfolios/"
  # global_data <- readRDS( file = paste0(wd, "Data/lOutput.rds") )
  env <- readRDS( file = paste0(wd, "Data/", universe, ".rds"))
  
  factor_names <- c("momentum", "value", "quality", "lowvola", "lowsize", "highsize")
  lFactor <- lapply( paste0(factor_names, "_mat"), FUN = get, envir = env$data_cons )
  # lFactor_stdz <- lapply( paste0(factor_names, "_stdz_mat"), FUN = get, envir = env$data_cons )
  
  
    
  
  
  
  
  
  