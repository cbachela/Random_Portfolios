
  
  ############################################################################
  ### GEOMSCALE - XXX - USA
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     17.05.2022
  # First version:    17.05.2022
  # --------------------------------------------------------------------------
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------

  require(garcholz)
  require(slolz)
  require(simolz)
  require(covolz)
  require(BSS)
  require(Backtest)
  require(RP)
  wd <- "R:/Asset_Management/Research_Projects/External_Research_Projects/GeomScale/Filters/"
  source( paste0(wd, "Source/class_BacktestCustom.R") )
  
  
  
  
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
  BT$spec$portfolio <- "dirichlet"
  
  
  scoreFUN <- function( rebdate = NULL, method = "xxx" )
  {
    eval( parse( text = "scoreFUN.", method, "(rebdate = ", rebdate, ")" ) )
  }
  
  scoreFUN.vola <- function( rebdate = NULL )
  {
    RBE <- rebalenv( rebdate = rebdate )
    X <- getData( RBE$GPS )
    sds <- apply( X, 2, sd )
    return( sds )
  }
  
  
  
  
  # --------------------------------------------------------------------------
  # RUN BACKTESTS
  # --------------------------------------------------------------------------
  
  BT1 <- BT$copy()
  # debugonce( BT1$accrejPortfolio )
  BT1$spec$portfolio <- "accrej"
  BT1$runLoop()
  
  
  
  
  
  
  
  
  
  
  
  BT2 <- BT$copy()
  BT2$spec$rebdates
  BT2$initRBE()
  BT2$spec$sharpe_filter_th <- 0.2
  # debugonce( BT2$setData )
  BT2$runLoop()
  
  BT2b <- BT$copy()
  BT2b$spec$rebdates
  BT2b$initRBE()
  BT2b$spec$sharpe_filter_th <- 0.8
  # debugonce( BT2$setData )
  BT2b$runLoop()
  
  
  
  
  BT3 <- BT$copy()
  BT3$spec$rebdates
  BT3$initRBE()
  BT3$spec$sharpe_filter_th <- 0.2
  BT3$spec$filter_by <- "mean"
  # debugonce( BT2$setData )
  BT3$runLoop()
  
 
  BT3b <- BT$copy()
  BT3b$spec$rebdates
  BT3b$initRBE()
  BT3b$spec$sharpe_filter_th <- 0.8
  BT3b$spec$filter_by <- "mean"
  # debugonce( BT2$setData )
  BT3b$runLoop()
  
  
  
  BT_tmp <- BT1
  lSim_tmp <- lapply( BT_tmp$output, FUN = function(wmat) { simPortfolio( X = BT_tmp$data$X_sim, wghts = wmat ) } )
  sim_bt1 <- do.call( cbindsimTS, lSim_tmp )
  colnames(sim_bt1) <- names(BT_tmp$output)
  stats_bt1 <- descStats(sim_bt1)$stats
  
  
  BT_tmp <- BT2
  lSim_tmp <- lapply( BT_tmp$output, FUN = function(wmat) { simPortfolio( X = BT_tmp$data$X_sim, wghts = wmat ) } )
  sim_bt2 <- do.call( cbindsimTS, lSim_tmp )
  colnames(sim_bt2) <- names(BT_tmp$output)
  stats_bt2 <- descStats(sim_bt2)$stats
  
  BT_tmp <- BT2b
  lSim_tmp <- lapply( BT_tmp$output, FUN = function(wmat) { simPortfolio( X = BT_tmp$data$X_sim, wghts = wmat ) } )
  sim_bt2b <- do.call( cbindsimTS, lSim_tmp )
  colnames(sim_bt2b) <- names(BT_tmp$output)
  stats_bt2b <- descStats(sim_bt2b)$stats
  
  
  BT_tmp <- BT3
  lSim_tmp <- lapply( BT_tmp$output, FUN = function(wmat) { simPortfolio( X = BT_tmp$data$X_sim, wghts = wmat ) } )
  sim_bt3 <- do.call( cbindsimTS, lSim_tmp )
  colnames(sim_bt3) <- names(BT_tmp$output)
  stats_bt3 <- descStats(sim_bt3)$stats
  
  BT_tmp <- BT3b
  lSim_tmp <- lapply( BT_tmp$output, FUN = function(wmat) { simPortfolio( X = BT_tmp$data$X_sim, wghts = wmat ) } )
  sim_bt3b <- do.call( cbindsimTS, lSim_tmp )
  colnames(sim_bt3b) <- names(BT_tmp$output)
  stats_bt3b <- descStats(sim_bt3b)$stats
  
  
  
  
  xlim <- range( c(stats_bt1["sds", ], stats_bt2["sds", ]) )
  ylim <- range( c(stats_bt1["means", ], stats_bt2["means", ]) )
  plot( x = stats_bt1["sds", ], y = stats_bt1["means", ], xlim = xlim, ylim = ylim )
  points( x = stats_bt2["sds", ], y = stats_bt2["means", ], col = 2 )
  points( x = stats_bt3["sds", ], y = stats_bt3["means", ], col = 3 )
  points( x = stats_bt2b["sds", ], y = stats_bt2b["means", ], col = 4 )
  points( x = stats_bt3b["sds", ], y = stats_bt3b["means", ], col = 5 )
  
  
  #// Note: excluding the worst might be a better strategy than including the best
  
  tmp <- cbind( stats_bt1["sharpe", ], 
                stats_bt2["sharpe", ],
                stats_bt2b["sharpe", ],
                stats_bt3["sharpe", ],
                stats_bt3b["sharpe", ] )
  boxplot( as.data.frame(tmp) )
  
  
  tmp <- cbind( stats_bt1["sharpe", -c(1:100)], 
                stats_bt2["sharpe", -c(1:100)], 
                stats_bt3["sharpe", -c(1:100)] )
  ldens <- apply( tmp, 2, density )
  plot.ldensity( ldens )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  