
  
  ############################################################################
  ### RANDOM PORTFOLIOS - THE CROSS-SECTION OF FACTOR SCORES - USA
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     10.06.2022
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
  
  
  # Load factor matrices from factor Shiny tool
  wd <- "R:/Asset_Management/R/Shiny/Factor_Portfolios/"
  # global_data <- readRDS( file = paste0(wd, "Data/lOutput.rds") )
  env <- readRDS( file = paste0(wd, "Data/", universe, ".rds"))
  
  
  # --------------------------------------------------------------------------
  # PREPARE BACKTEST OBJECT
  # --------------------------------------------------------------------------
  
  BT <- BacktestCustom$new()
  BT$spec <- BT0$spec
  BT$spec$selection_filter <- c(BT0$spec$selection_filter, "NA")
  BT$data <- BT0$data
  BT$data$upper_mat <- BT$data$upper_mat * 0 + 1
  BT$spec$fev_bias <- 2:10
  BT$spec$n_sim <- 10
  BT$spec$portfolio <- "dirichlet"
  BT$spec$sharpe_filter_th <- NULL
  
  
  BT1 <- BT$copy()
  # BT1$runLoop()
  
  
  
  
  # --------------------------------------------------------------------------
  # PREPARE DATA FOR PARTICULAR DATE
  # --------------------------------------------------------------------------
  
  # today <- "2019-01-16"
  today <- "2008-02-27"
  BT1$gps( rebdate = today )
  BT1$gpo( rebdate = today )
  RBE <- rebalenv( rebdate = today )
  GPS <- RBE$GPS
  selection <- RBE$selection
  X <- getData(GPS)
  X_bm <- BT1$data$X_bm[rownames(X), ]
  
  
  # Samples
  n_sim <- 10^4
  w_bm <- setNames( as.numeric(BT0$data$wmat_bm[today, selection]), selection )
  w_bm <- w_bm / sum(w_bm)
  w_eqw <- rep(1/length(selection), length(selection))
  
  # Flat Dirichlet
  n <-  length(selection)
  alpha <- rep(1, n)
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  
  
  
  # m out of n bootstrap
  m <- 10
  prob <- rep(1/n, n)
  gridnodes <- matrix(0, nrow = n_sim, ncol = n)
  for ( i in 1:n_sim ) {
    idx <- sample( 1:n, size = m, replace = TRUE, prob = prob )
    tbl <- table(idx)
    gridnodes[i, as.numeric(names(tbl))] <- tbl / m
  }
  norm2_boot_m <- apply( gridnodes, 1, function(x) { sum(x^2) } )
 
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # THE DARTBOARD CONTEST
  # --------------------------------------------------------------------------
  
  source(file = "H:/R/notingit/Bootstrap/Source/custom_functions.R" )
  
  mu <- meanGeo( X[ ,names(w_bm)], scalefactor = 52 )
  n <- length(mu)
  m <- 4
  set.seed(1111)
  x <- mu[sample(1:n, size = m, replace = TRUE)]
  
  # Classical n out of n bootstrap
  M2Exact( z = mu, exp2norm2 = 2/n - 1/n^2 )
  
  # m out of n bootstrap
  mu2_moon <- M2Exact( z = mu, exp2norm2 = (n + m - 1) / (n * m) )
  
  # Same with corresponging lambda for BB 
  lambda <- mn2Lambda( n = n, m = m )
  mu2_lambda <- M2Exact( z = mu, exp2norm2 = (1 + lambda) / (1 + lambda * n) )
  
  # Classical bootstrap on the m subset
  mu2_m <- M2Exact( z = x, exp2norm2 = 2/m - 1/m^2 )
  
  mu2 <- c(x_bar = var(mu), moon = mu2_moon, lambda = mu2_lambda, m = mu2_m)
  sqrt(mu2)
  
  
  
  plot( density(mu) )
  abline( v = mean(mu), col = "grey" )
  abline( v = median(mu), col = "green" )
  abline( v = meanGeo(X_bm, scalefactor = 52) )
  abline( v = sum(w_bm * mu), col = 2 )
  
  sum( mu > meanGeo(X_bm, scalefactor = 52) ) / length(mu)
  sum( mu > mean(mu) ) / length(mu)
  sum( mu > median(mu) ) / length(mu)
  
  
  
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # ALPHABET SCORE
  # --------------------------------------------------------------------------

  Names <- setNames( as.character( BT0$data$name_mat[today, ] ), 
                     colnames(BT0$data$name_mat) )
  names_vec <- toupper(Names[names(w_bm)])
  letter_vec <- substr( aviationAlphabet(), 1, 1 )
  lScore_vec <- list()
  for ( i in seq(along = letter_vec ) ) {
    lScore_vec[[i]] <- unlist( lapply( names_vec, letterScore, preference = letter_vec[i] ) )
  }
  names(lScore_vec) <- letter_vec
  
  # Inverse score
  lScore_vec_inv <- list()
  for ( i in seq(along = letter_vec ) ) {
    lScore_vec_inv[[i]] <- as.numeric( lScore_vec[[i]] == 0 )
  }
  names(lScore_vec_inv) <- letter_vec
  
  
  
  # Letter "Z":
  score_vec <- lScore_vec[["Z"]]
  score_vec <- lScore_vec_inv[["Z"]]
  
  
  # Dirichlet
  theta <- apply( samples, 1, function(w) { sum(w * score_vec) } )
  # m out of n bootstrap
  theta_moonboot <- apply( gridnodes, 1, function(w) { sum(w * score_vec) } )
  
  plot(densFUN(score_vec), col = "darkgrey" )  
  abline( v =  sum(w_bm * score_vec), col = "blue" )
  abline( v = sum(w_eqw * score_vec) )
  lines( densFUN(theta) )
  lines( densFUN(theta_moonboot), col = 2 )
  
  
  ldens <- lapply( lScore_vec, FUN = densFUN )
  ldens <- lapply( lScore_vec_inv, FUN = densFUN )
  slolz:::plot.ldensity(ldens, fillin = FALSE, cex = 0.6, legend = "topright" )
  
  
  
  # Distribution of mean 
  mu <- meanGeo( X[ ,names(w_bm)], scalefactor = 52 )
  score_vec <- lScore_vec[["Z"]]
  w_portf <- score_vec / sum(score_vec)
  mu_portf <- meanGeo( timeSeries( X[ ,names(w_bm)] %*% w_portf, time(X) ) )
  m <- sum(score_vec > 0)
  m
  
  # moon mapped Dirichlet
  n <- length(mu)
  lambda <- mn2Lambda( n = n, m = m )
  alpha <- rep(lambda, n)
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  theta <- apply( samples, 1, function(w) { sum(w * mu) } )
  
  
  plot( density(theta) )
  abline( v = mu_portf, col = 1 )
  
  
  # Repeat analysis over entire 20y backtest
  # ...
  
  
  
  
  
  # --------------------------------------------------------------------------
  # MARKET BETA
  # --------------------------------------------------------------------------
  
  X_bm_tmp <- X_bm
  # X_bm_tmp <- timeSeries( X[ ,names(w_bm)] %*% w_bm, time(X) )
  
  sigma_sq_bm <- var(X_bm_tmp)
  beta_score <- apply( X[ ,names(w_bm)], 2, function(x) { cov(x, X_bm_tmp) / sigma_sq_bm } )
  sim <- timeSeries( X[ ,names(w_bm)] %*% cbind(t(samples), w_bm), time(X) )
  sim_moonboot <- timeSeries( X[ ,names(w_bm)] %*% cbind(t(gridnodes), w_bm), time(X) )
  
  theta <- apply( sim, 2, function(y) { cov(y, X_bm_tmp) / sigma_sq_bm } )
  theta_moonboot <- apply( sim_moonboot, 2, function(y) { cov(y, X_bm_tmp) / sigma_sq_bm } )
  
  
  plot( density(beta_score), col = "darkgrey" )  
  abline( v = mean(beta_score), col = "darkgrey" )
  abline( v =  sum(w_bm * beta_score), col = "blue" )
  abline( v =  theta[length(theta)], col = "blue" )
  abline( v = sum(w_eqw * beta_score) )
  lines( density(theta) )
  lines( density(theta_moonboot), col = 2 )
  
  
  
  
  # --------------------------------------------------------------------------
  # SIZE
  # --------------------------------------------------------------------------
  
  size_score <- w_bm
  
  sum(w_bm * size_score)
  sum(w_eqw * size_score)
  

  theta <- apply( samples, 1, function(w) { sum(w * size_score) } )
  # m out of n bootstrap
  theta_moonboot <- apply( gridnodes, 1, function(w) { sum(w * size_score) } )
  
  plot(density(size_score), col = "darkgrey" )  
  abline( v =  sum(w_bm * size_score), col = "blue" )
  abline( v = sum(w_eqw * size_score) )
  lines( density(theta) )
  lines( density(theta_moonboot), col = 2 )
  
  
  lX <- list( score = size_score,
              bb = theta,
              moonboot = theta_moonboot )
  ldens <- lapply( lX, density )
  colors <- c("darkgrey", "black", "red")
  slolz:::plot.ldensity( ldens, fillin = FALSE, col = colors, legen = "topright" )
  abline( v =  sum(w_bm * size_score), col = "blue" )
  abline( v = sum(w_eqw * size_score) )
  
  
  
  # Power law?
  market_cap <- 10^10
  size <- size_score * market_cap
  dens <- density( size )
  idx <- intersect( which(dens$y > 0), which(dens$x > 0) )
  plot( y = log(dens$y[idx]), x = log(dens$x[idx]) )
  plot( y = cumsum(dens$y[idx] / sum(dens$y[idx])), x = dens$x[idx] )
  plot( y = log( cumsum(dens$y[idx] / sum(dens$y[idx])) ), x = log( dens$x[idx]) )
  
  
  
  ldens_adj <- lapply( ldens, FUN = function(dens) 
  { 
    idx <- intersect( which(dens$y > 0), which(dens$x > 0) )
    ans <- list(y = log( cumsum(dens$y[idx] / sum(dens$y[idx])) ), x = log(dens$x[idx]) ) 
  } )
  slolz:::plot.ldensity( ldens_adj, fillin = FALSE, col = colors, legen = "topright" )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # VALUE
  # --------------------------------------------------------------------------
  
  value_score <- setNames( as.numeric(env$data_cons$value_mat[today, selection]), selection )
  value_score[is.na(value_score)] <- 0
  
  theta <- apply( samples, 1, function(w) { sum(w * value_score) } )  
  # m out of n bootstrap
  theta_moonboot <- apply( gridnodes, 1, function(w) { sum(w * value_score) } )
  
  
  plot(density(value_score), col = "darkgrey" )  
  abline( v =  sum(w_bm * value_score), col = "blue" )
  abline( v = sum(w_eqw * value_score) )
  lines( density(theta) )
  lines( density(theta_moonboot), col = 2 )
  
  
  
  
  # --------------------------------------------------------------------------
  # QUALITY
  # --------------------------------------------------------------------------
  
  quality_score <- setNames( as.numeric(env$data_cons$quality_mat[today, selection]), selection )
  quality_score[is.na(quality_score)] <- 0
  
  # Names <- intersect(names(w_bm), names(quality_score))
  # sum(w_bm[Names] * quality_score[Names])
  # sum(w_eqw[Names] * quality_score[Names])
  
  sum(w_bm * quality_score)
  sum(w_eqw * quality_score)
  
  
  theta <- apply( samples, 1, function(w) { sum(w * quality_score) } )  
  # m out of n bootstrap
  theta_moonboot <- apply( gridnodes, 1, function(w) { sum(w * quality_score) } )
  
  
  plot(density(quality_score), col = "darkgrey" )  
  abline( v =  sum(w_bm * quality_score), col = "blue" )
  abline( v = sum(w_eqw * quality_score) )
  lines( density(theta) )
  lines( density(theta_moonboot), col = 2 )
  
  
  
  # --------------------------------------------------------------------------
  # MOMENTUM
  # --------------------------------------------------------------------------
  
  mom_score <- apply( X, 2, mean )

  sum( w_bm * mom_score )
  sum( w_eqw * mom_score )
  
  theta <- apply( samples, 1, function(w) { sum(w * mom_score) } )    
  
  plot(density(mom_score), col = "darkgrey" )  
  abline( v =  sum(w_bm * mom_score), col = "blue" )
  abline( v = sum(w_eqw * mom_score) )
  lines( density(theta) )

  
    
  # --------------------------------------------------------------------------
  # VOLA
  # --------------------------------------------------------------------------
  
  vola_score <- apply( X, 2, sd )^2
  
  sum( w_bm * vola_score )
  sum( w_eqw * vola_score )
  
  theta <- apply( samples, 1, function(w) { sum(w * vola_score) } )    

  plot(density(vola_score), col = "darkgrey" )  
  abline( v =  sum(w_bm * vola_score), col = "blue" )
  abline( v = sum(w_eqw * vola_score) )
  lines( density(theta) )
  
  
  lX <- list( cross_section = vola_score,
              random_portfolio = theta )
  ldens <- lapply( lX, density )
  colors <- c("black", "blue")
  slolz:::plot.ldensity( ldens, col = colors, legen = "topright" )
  abline( v =  sum(w_bm * vola_score), col = "red", lwd = 2 )

  
  
  # Power law?
  dens <- density( vola_score )
  idx <- intersect( which(dens$y > 0), which(dens$x > 0) )
  plot( y = log(dens$y[idx]), x = log(dens$x[idx]) )
  plot( y = cumsum(dens$y[idx] / sum(dens$y[idx])), x = dens$x[idx] )
  plot( y = log( cumsum(dens$y[idx] / sum(dens$y[idx])) ), x = log( dens$x[idx]) )
  
  
  
  ldens_adj <- lapply( ldens, FUN = function(dens) 
  { 
    idx <- intersect( which(dens$y > 0), which(dens$x > 0) )
    ans <- list(y = log( cumsum(dens$y[idx] / sum(dens$y[idx])) ), x = log(dens$x[idx]) ) 
  } )
  slolz:::plot.ldensity( ldens_adj, fillin = FALSE, col = colors, legen = "topright" )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # SHARPE
  # --------------------------------------------------------------------------
  
  sharpe_score <- descStats(X)$stats["sharpe", ]
  
  sum( w_bm * mom_score )
  sum( w_eqw * mom_score )
  
  theta <- apply( samples, 1, function(w) { sum(w * sharpe_score) } )    
  
  plot(density(sharpe_score))  
  abline( v =  sum(w_bm * sharpe_score), col = "blue" )
  abline( v = sum(w_eqw * sharpe_score), col = "grey" )
  lines( density(theta), col = "grey" )
  
  
  
  
  # --------------------------------------------------------------------------
  # ALPHA
  # --------------------------------------------------------------------------
  
  path_ff <- "R:/Asset_Management/R/myRData/Factor/"
  env_ff <- readRDS( file = paste0(path_ff, "data.Rds") )
  
  ls(env_ff)  
  ls(env_ff$lFactor)
  
  # Alpha's of single stocks
  FF3 <- env_ff$lFactor$`3F`$USA_daily
  dates <- intersect(rownames(FF3), rownames(X))
  X_train <- FF3[dates, 1:3]
  lReg <- list()
  for ( j in 1:ncol(X) ) {
    reg <- regression( Y_train = X[dates, j], 
                       X_train = X_train, 
                       type = "ols" )
    lReg[[j]] <- reg$coeffmat
  }  
  alpha_score <- unlist( lapply( lReg, FUN = function(x) { x[1, 1] } ) )  
  
  theta <- apply( samples, 1, function(w) { sum(w * alpha_score) } )   
  
  plot(density(alpha_score))  
  abline( v =  sum(w_bm * alpha_score) )
  abline( v = sum(w_eqw * alpha_score), col = "grey" )
  lines( density(theta), col = "grey" )
  
  
  
  # Alpha's of random portfolios
  Y <- timeSeries( X[dates, ] %*% t(samples), dates )
  lRegY <- list()
  for ( j in 1:nrow(Y) ) {
    reg <- regression( Y_train = Y[ ,j], 
                       X_train = X_train, 
                       type = "ols" )
    lRegY[[j]] <- reg$coeffmat
  }  
  alpha_score_y <- unlist( lapply( lRegY, FUN = function(x) { x[1, 1] } ) )  
  
  
  
  plot(density(alpha_score))  
  abline( v =  sum(w_bm * alpha_score) )
  abline( v = sum(w_eqw * alpha_score), col = "grey" )
  lines( density(theta), col = "grey" )
  lines( density(alpha_score_y), col = "red" )
  
  
  
  
  t_score <- unlist( lapply( lReg, FUN = function(x) { x[1, 3] } ) )  
  t_score_y <- unlist( lapply( lRegY, FUN = function(x) { x[1, 3] } ) )  
  
  plot(density(t_score))  
  lines( density(t_score_y), col = "red" )
  abline( v = 2 )
    
  
  
  p_score <- unlist( lapply( lReg, FUN = function(x) { x[1, 4] } ) )
  p_score_y <- unlist( lapply( lRegY, FUN = function(x) { x[1, 4] } ) )
  
  plot(density(p_score))  
  lines( density(p_score_y), col = "red" )
  abline( v = 0.05 )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  