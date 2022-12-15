
  
  ############################################################################
  ### RANDOM PORTFOLIOS - ALPHABET STRATEGIES - USA
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     13.11.2022
  # First version:    13.11.2022
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
  # FUNCTIONS
  # --------------------------------------------------------------------------
  
  simFUN <- function(wmat) 
  { 
    simPortfolio( X = BT0$data$X_sim, 
                  wghts = wmat,
                  fc = 0, 
                  vc = 0,
                  language = "C" )
  }
  
  
  
  
  # --------------------------------------------------------------------------
  # PARAMETERS
  # --------------------------------------------------------------------------
  
  universe <- "usa"
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  
  
  
  # --------------------------------------------------------------------------
  # LOAD INITIAL BACKTEST OBJECT
  # --------------------------------------------------------------------------
  
  BT0 <- loadBacktest( wd = paste0(wd, "Data/"),
                       name = "usa_maxw5" )
  BT0$initRBE()
  
  
  
  
  # --------------------------------------------------------------------------
  # ALPHABET STRATEGIES
  # --------------------------------------------------------------------------
  
  env <- readRDS( file = paste0(wd, "waRehouse/alphabet_strategy.rds") )
  lSim <- env$lSim
  X_alphabet <- do.call( cbind, lSim )
  
  
  
  # --------------------------------------------------------------------------
  # BENCHMARK SERIES
  # --------------------------------------------------------------------------
  
  tickers <- setNames( c("NDDUUS", "M1USMMT", "M1US000V", "M1USEV", "M1USQU", "M00IMVSO", "M1USSZT", "M1WODV"),
                       c("Capw", "Momentum", "Value", "Enhanced_Value", "Quality", "Min_Var", "Size", "Multifactor") )
  X_fact <- rodbcGetOLZDBReturns( assetName = tickers,
                                  refCcy = "USD",
                                  frqncy = "daily",
                                  na.rm = "r" )
  X_fact <- X_fact[isWeekday(time(X_fact)), tickers]
  colnames(X_fact) <- names(tickers)
  
  
  # Combine
  X_bm <- na.omit( cbind( X_fact, EQW = env$sim_eqw, X_alphabet ) )
  
  
  
  # --------------------------------------------------------------------------
  # RANDOM PORTFOLIOS - UNIFORM
  # --------------------------------------------------------------------------
  
  BT2 <- BacktestCustom$new()
  BT2$data <- BT0$data
  BT2$spec <- BT0$spec
  BT2$spec$portfolio <- "momentum_rp"
  BT2$spec$n_sim <- 10^2 * 5
  BT2$spec$sampling_dist = "uniform"
  # BT2$spec$m <- NULL # means m = length(selection)
  BT2$spec$m <- 30
  BT2$spec$shadow_dirichlet_transform <- FALSE
  BT2$spec$scl_by_capw <- FALSE
  BT2$runLoop()  
  
  # Simulate
  tmp <- lapply( BT2$output, FUN = simFUN )
  sim_2 <- do.call( cbind, tmp )
  colnames(sim_2) <- paste0("rp_uniform_", 1:ncol(sim_2))
  
  
 
  # --------------------------------------------------------------------------
  # RANDOM PORTFOLIOS - CAPW
  # --------------------------------------------------------------------------
  
  BT3 <- BacktestCustom$new()
  BT3$data <- BT0$data
  BT3$spec <- BT0$spec
  BT3$spec$portfolio <- "momentum_rp"
  BT3$spec$n_sim <- 10^2 * 5
  # BT3$spec$sampling_dist = "uniform"
  BT3$spec$sampling_dist = "capw"
  # BT3$spec$m <- NULL # means m = length(selection)
  BT3$spec$m <- 30
  # BT3$spec$th <- 0.1
  BT3$spec$shadow_dirichlet_transform <- FALSE
  BT3$spec$scl_by_capw <- TRUE # ~~~~~~~~~~~~~~~ 
  # BT3$spec$n_sim <- 10^3
  # debugonce( BT3$momentum_rpPortfolio )
  # debugonce( BT3$gps )
  # debugonce( BT3$selection.db_flag )
  BT3$runLoop()  
  
  # Simulate
  tmp <- lapply( BT3$output, FUN = simFUN )
  sim_3 <- do.call( cbind, tmp )
  colnames(sim_3) <- paste0("rp_capw", 1:ncol(sim_3))
  
  
 
  
  
  
  
  
  
  # Combine
  end_date <- "2022-03-06"
  strategy_names <- c("Capw", "EQW", colnames(X_alphabet))
  # sim_all <- na.omit( cbind( X_bm[ ,strategy_names], sim_2) )
  sim_all <- na.omit( cbind( X_bm[ ,strategy_names], sim_3) )
  # sim_all <- na.omit( cbind( X_bm[ ,strategy_names], sim_2, sim_3) )
  sim_all <- window( sim_all, start(sim_all), end_date )
  lStats <-  descStats( sim_all, descStatsSpec(annualize = TRUE) )
  # statsolz:::plot.stats( descStats( na.omit( cbind( X_bm[ ,strategy_names], X_alphabet) ) ), sortby = NULL )
  
  
  # First decade
  sim_all_decade1 <- window( sim_all, start(sim_all), "2011-01-01" )
  lStats_decade1 <-  descStats( sim_all_decade1, descStatsSpec(annualize = TRUE) )
  
  # Second decade
  sim_all_decade2 <- window( sim_all, "2011-01-01", end(sim_all) )
  lStats_decade2 <-  descStats( sim_all_decade2, descStatsSpec(annualize = TRUE) )
  
  # Bear Bull 
  BBSObj <- BBSRC$new()
  BBSObj$setCtrl()
  BBSObj$data$X_level <- cumulated( BT0$data$X_bm, "discrete" )
  BBSObj$runRobust()
  BBSObj$data$X <- sim_all
  phase_stats <- BBSObj$phaseStats()
  
  
  
  
  
  ###
  # Risk return plot
  
  df <- data.frame( x = lStats$stats["sds", ],
                    y = lStats$stats["means", ] )
  df <- data.frame( x = lStats_decade1$stats["sds", ],
                    y = lStats_decade1$stats["means", ] )
  df <- data.frame( x = lStats_decade2$stats["sds", ],
                    y = lStats_decade2$stats["means", ] )
  df <- data.frame( x = phase_stats$'states==1'["sds", ],
                    y = phase_stats$'states==1'["means", ] )
  df <- data.frame( x = phase_stats$'states==-1'["sds", ],
                    y = phase_stats$'states==-1'["means", ] )
  rownames(df) <- colnames(sim_all)
  h1 <- hist(df$x, breaks=250, plot=F)
  h2 <- hist(df$y, breaks=250, plot=F)
  top <- max(h1$counts, h2$counts) * 2
  lims_x <- range(c(df$x * 0.995, df$x * 1.005))
  lims_y <- range(df$y)
  k <- kde2d( df$x, df$y, n = 250, lims = c(lims_x, lims_y) )
  # margins
  oldpar <- par()
  par(mar=c(3,3,1,1))
  layout(matrix(c(2,0,1,3),2,2,byrow=TRUE),c(3,1), c(1,3))
  image(k)
  points( x = df$x, y = df$y, pch = 19, col = "orange", cex = 0.5 )
  tmp <- MASS::kde2d( x = df$x, y = df$y, n = 25, h = c(0.005, 0.02), lims = c(range(df$x), range(df$y)) ) 
  # contour( tmp, xlab = "Portfolio Standard Deviation", ylab = "Portfolio Return", 
  #          add = TRUE, drawlabels = FALSE, col = "red" )
  colors_bm <- c("black", "blue")
  colors <- c(colors_bm, rep("grey", ncol(X_alphabet)))
  for ( i in 1:length(colors) ) {
    points( x = df$x[i], y = df$y[i], pch = 19, cex = 2, col = colors[i] )
  }
  text( df[colnames(X_alphabet), "x"], 
        df[colnames(X_alphabet), "y"], 
        colnames(X_alphabet), font = 2, cex = 1, col = "black" )
  legend( "topright", strategy_names[1:2], pch = 19, col = colors_bm, text.col = colors_bm, bty = "n")
  par( mar=c(0,2,1,0) )
  barplot( h1$counts, axes=FALSE, ylim = c(0, top), space=0, col='grey' )
  par( mar=c(2,0,0.5,1) )
  barplot( h2$counts, axes=FALSE, xlim = c(0, top), space=0, col='grey', horiz=TRUE )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # ALPHA
  # --------------------------------------------------------------------------
  
  path_ff <- "R:/Asset_Management/R/myRData/Factor/"
  env_ff <- readRDS( file = paste0(path_ff, "data.Rds") )
  FF3 <- env_ff$lFactor$`3F`$USA_daily
  FF5 <- env_ff$lFactor$`5F`$USA_daily
 
  # Alpha's of naive random portfolios
  dates <- intersect( rownames(X_bm), rownames(FF5) )
  # dates <- intersect( dates, rownames(BBSObj$output$states)[ BBSObj$output$states == -1 ] )
  # dates <- intersect( dates, rownames(BBSObj$output$states)[ BBSObj$output$states == 1 ] )
  Y <- sim_2[dates, ]
  X_train <- FF5[dates, 1:5]
  lReg <- list()
  for ( j in 1:ncol(Y) ) {
    reg <- regression( Y_train = Y[ ,j], 
                       X_train = X_train, 
                       type = "ols" )
    lReg[[j]] <- reg$coeffmat
  }  
  alpha_score <- unlist( lapply( lReg, FUN = function(x) { x[1, 1] } ) )  
  t_score <- unlist( lapply( lReg, FUN = function(x) { x[1, 3] } ) )  
  p_score <- unlist( lapply( lReg, FUN = function(x) { x[1, 4] } ) )
  

  
  # Alpha's of simple random portfolios
  dates <- intersect( rownames(X_bm), rownames(FF5) )
  Y <- sim_3[dates, ]
  X_train <- FF5[dates, 1:5]
  lReg_simple <- list()
  for ( j in 1:ncol(Y) ) {
    reg <- regression( Y_train = Y[ ,j], 
                       X_train = X_train, 
                       type = "ols" )
    lReg_simple[[j]] <- reg$coeffmat
  }  
  alpha_score_simple <- unlist( lapply( lReg_simple, FUN = function(x) { x[1, 1] } ) )  
  t_score_simple <- unlist( lapply( lReg_simple, FUN = function(x) { x[1, 3] } ) )  
  p_score_simple <- unlist( lapply( lReg_simple, FUN = function(x) { x[1, 4] } ) )
  
  
  plot(density(t_score))  
  abline( v = 2 )
  plot(density(alpha_score))
  plot(density(p_score))  
  
  sum( p_score < 0.05 ) / length(p_score)
  sum( alpha_score < 0 ) / length(alpha_score)
  
  sum( p_score_simple < 0.05 ) / length(p_score_simple)
  sum( alpha_score_simple < 0 ) / length(alpha_score_simple)
  
  mean(alpha_score_simple)
  
  
  
  
  
  env <- readRDS( file = paste0(wd, "waRehouse/alphabet_strategy_factor_analysis.rds") )
  ldens <- list( alphabet_alphas = density(env$alpha_score), simple_rp = density(alpha_score_simple) )
  slolz:::plot.ldensity( ldens )
  
  
  
  
  
  
  
  
  

  
  # Relative performance
  FUN <- function(i) { simOutperformance( x = sim_all[ ,strategy_name], 
                                          y = sim_all[ ,i] ) }
  lSim_delta <- list()
  for ( strategy_name in colnames(X_bm) ) {
    lSim_delta_tmp <- lapply( (ncol(X_bm)+1):ncol(sim_all), FUN = FUN )
    lSim_delta[[strategy_name]] <- do.call( cbind, lSim_delta_tmp )
  }
  names(lSim_delta)
  
  # plot( log(cumulated(sim_delta, "discrete")), plot.type = "single", col = "orange" )
  # abline( h = 0, col = "grey" )
  
 
  
  
  
  # # Rolling relative performance
  # sim_delta <- lSim_delta[["Capw"]]
  # stats_roll <- applyRoll( Data = cumulated(sim_delta, "discrete"),
  #                          Width = 252 * 1,
  #                          By = 1,
  #                          FUN = function(x) { as.numeric(x[nrow(x), ]) / as.numeric(x[1, ]) - 1 } )
  # 
  # plot( stats_roll, plot.type = "single", col = "orange" )                                                    
  # abline( h = 0, col = "grey" )
  # 
  # 
  # perc_outperf <- apply( stats_roll, 1, function(x) { sum(x > 0) / length(x) } )
  # plot( perc_outperf )
  # abline(h = 0.5)
  
  
  
  # Market regimes
   
  BBSObj <- BBSRC$new()
  BBSObj$setCtrl()
  BBSObj$setData()
  BBSObj$runRobust()
  states <- BBSObj$output$states
  lPhase_stats <- list()
  for ( strategy_name in names(lSim_delta) ) {
    BBSObj$data$X <- lSim_delta[[strategy_name]]
    lPhase_stats[[strategy_name]] <- BBSObj$phaseStats()
  }
  
  # phase_stats <- lPhase_stats[[1]]
  # boxplot( as.data.frame( cbind( bull = phase_stats$`states==1`["means", ],
  #                                bear = phase_stats$`states==-1`["means", ] ) ), 
  #          beside = TRUE, col = c("green", "red") )
  
  
  tmp <- lapply( lPhase_stats, function(x) { cbind( bull = x$`states==1`["means", ],
                                                    bear = x$`states==-1`["means", ] ) } )
  tmp <- lapply( names(tmp), function(name) { x <- tmp[[name]]; colnames(x) <- paste0(name, "_", colnames(x)); return(x) } )
  tmp <- do.call( cbind, tmp )
  
  boxplot( as.data.frame( tmp[ ,3:18] ), beside = TRUE, col = c("green", "red") )
  abline( h = 0 )
  
  boxplot( as.data.frame( tmp ), beside = TRUE, col = c("green", "red") )
  abline( h = 0 )
  
 

  
  
  