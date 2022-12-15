
  
  ############################################################################
  ### RANDOM PORTFOLIOS - THE CROSS-SECTION OF PORTFOLIO RETURNS - USA
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     22.06.2022
  # First version:    10.06.2022
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
  BT$data <- BT0$data
  BT$data$upper_mat <- BT$data$upper_mat * 0 + 1
  BT$spec$fev_bias <- 2:10
  BT$spec$n_sim <- 10
  BT$spec$portfolio <- "dirichlet"
  BT$spec$sharpe_filter_th <- NULL
  BT$spec$GPS@covariance$method <- "duv"
  
  
  BT1 <- BT$copy()
  # BT1$runLoop()
  
  
  
  # --------------------------------------------------------------------------
  # PREPARE DATA FOR PARTICULAR DATE
  # --------------------------------------------------------------------------
  
  # BT1$spec$rebdates
  # today <- "2008-08-27" 
  # tomorrow <- "2008-11-26"
  today <- "2005-03-02" 
  tomorrow <- "2005-06-01"
  BT1$gps( rebdate = today )
  BT1$gpo( rebdate = today )
  RBE <- rebalenv( rebdate = today )
  GPS <- RBE$GPS
  selection <- RBE$selection
  # X <- getData(GPS)
  # X_bm <- BT1$data$X_bm[rownames(X), ]
  
  dates <- rownames(BT1$data$X_sim)
  dates <- dates[ dates > today ]
  dates <- dates[ dates <= tomorrow ]
  X <- BT1$data$X_sim[dates, selection]
  X[is.na(X)] <- 0
  X_bm <- BT1$data$X_bm[dates, ]
  
  
  # Samples
  w_bm <- setNames( as.numeric(BT0$data$wmat_bm[today, selection]), selection )
  # w_bm <- w_bm / sum(w_bm)
  w_eqw <- rep(1/length(selection), length(selection))
  n_sim <- 10^4
  alpha <- rep(1, length(selection), length(selection))
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # RETURNS
  # --------------------------------------------------------------------------
  
  stats_bm <- descStats( X_bm )$stats
  stats <- descStats(X)$stats
  field <- "means"
  mu <- stats[field, ]

  sum( w_bm * mu )
  sum( w_eqw * mu )
  stats_bm[field, ]
  
  
  theta <- apply( samples, 1, function(w) { sum(w * mu) } )    
  
  plot(density(mu))  
  abline( v =  sum(w_bm * mu), col = "red" )
  abline( v = sum(w_eqw * mu), col = "black" )
  abline( v = stats_bm[field, ], col = "darkred" )
  lines( density(theta), col = "darkgrey" )

  
  
    
  
  
  
  
  # Dates where there is a change in the index composition
  
  IM <- BT0$data$Index_Member
  IM[is.na(IM)] <- 0
  review_dates <- rownames(IM)[1]
  for ( i in 2:nrow(IM) ) {
    a <- setNames( as.numeric(IM[i, ]), colnames(IM) )
    b <- setNames( as.numeric(IM[i-1, ]), colnames(IM) )
    if ( !all( a == b ) ) {
      review_dates <- c(review_dates, rownames(IM)[i])
    }
  }
  
  
  
  
  
  
  