
  
  ############################################################################
  ### RANDOM PORTFOLIOS - CLASS BACKTEST INDEX REPLICATION
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard   01.03.2023
  # First version:    01.03.2023
  # --------------------------------------------------------------------------

  
  # BacktestIndexReplication.setData
  # BacktestIndexReplication.gps
  # BacktestIndexReplication.gpo
  # BacktestIndexReplication.nnlsPortfolio
  # BacktestIndexReplication.olsPortfolio
  # BacktestIndexReplication.lsPortfolio
  # BacktestIndexReplication.selection.100
  
  
  

  require(Backtest)
  
  
  
  
  # --------------------------------------------------------------------------
  BacktestIndexReplication <- setRefClass( Class = "BacktestIndexReplication", 
                                           contains = "BacktestPanel",
                                           methods = list() )
  
 
  
  # --------------------------------------------------------------------------
  BacktestIndexReplication.gps <- function( rebdate = NULL )
  {
    .self$setGPS(rebdate = rebdate)
    .self$setSelection(rebdate = rebdate)
    .self$selection.100(rebdate = rebdate) # ~~~~~~~~~~~~~~~~~~~~~~  
    .self$setData(rebdate = rebdate)
    .self$setCovariance(rebdate = rebdate)
    .self$setConstraints(rebdate = rebdate)
    return(TRUE)
  }
  BacktestIndexReplication$methods( gps = BacktestIndexReplication.gps )
  
    
  
  # --------------------------------------------------------------------------
  BacktestIndexReplication.gpo <- function( rebdate = NULL )
  {
    portfolio_fun <- spec$portfolio
    
    if ( !is.null(portfolio_fun) ) {
      
      txt <- paste0(".self$", portfolio_fun, "Portfolio( rebdate = '", rebdate, "')" )
      eval( parse( text = txt ) )
    
    } else {
      
      tic <- Sys.time()
      # .self$initialWeights( rebdate = rebdate )
      # .self$floatUpperBounds( rebdate = rebdate )
      RBE <- rebalenv( rebdate = rebdate )
      RBE$GPO <- GPO::gpo( GPS = RBE$GPS )
      toc <- Sys.time() - tic
      wghts <- matrix( getWeights(RBE$GPO), nrow = 1, 
                       dimnames = list(rebdate, names(getWeights(RBE$GPO))) )
      tictoc <- matrix(toc, dimnames = list(rebdate, "tictoc"))
      
      # Append output
      .self$appendOutput( list( weights = as.timeSeries(wghts),
                                runtime = as.timeSeries(tictoc) ) )
      
      # Clean rebalancing environment
      if ( isTRUE(spec$clean_RBE) && exists("RBE") ) {
        cleanEnvir( keep = "rebdate", envir = RBE )
      }
    }
    
    return( TRUE )
  }
  BacktestIndexReplication$methods( gpo = BacktestIndexReplication.gpo )
  
  
  
  # --------------------------------------------------------------------------
  BacktestIndexReplication.nnlsPortfolio <- function( rebdate = NULL )
  {
    
    RBE <- rebalenv( rebdate = rebdate )
    GPS <- RBE$GPS
    X_w <- getData( GPS )
    X_d <- window(data$X_sim[ ,colnames(X_w)], start(X_w), end(X_w))  # Use daily returns  -- ccy?
    X_d <- X_d[isWeekday(time(X_d)), ]
    X_tmp <- t( na.omit( t(X_d) ) )
    X <- timeSeries( X_tmp, rownames(X_d) )
    dates <- intersect( rownames(data$X_bm), rownames(X) )
    y <- log(1 + data$X_bm[dates, ])
    X <- log(1 + X[dates, ])
    colnames(y) <- "RV"
    colnames(X) <- paste0("EV_", colnames(X))
    
    
    # # Check
    # w <- rep(1/ncol(X), ncol(X))
    # ans <- y * NA
    # ans2 <- ans
    # for ( today in dates ) {
    #   a <- as.numeric(X[today, ])
    #   ans[today, ] <- (t(a) %*% w - y[today, ])^2
    #   A <- a %*% t(a)
    #   ans2[today, ] <- t(w) %*% A %*% w - t(2*y[today, ]*a) %*% w + y[today, ]^2
    # }
    # sum(ans)
    # sum(ans2)  # same same
    # sum( (X %*% w - y)^2 ) # same same
    # t(w) %*% (t(X) %*% X) %*% w - t(y) %*% X %*% w - t(w) %*% t(X) %*% y - t(y) %*% y  # different ?
    # t(w) %*% (t(X) %*% X) %*% w - 2 * t(y) %*% X %*% w - t(y) %*% y
    
    # Weights
    tau <- spec$tau
    wt <- y * 0 + 1
    if ( !is.null(tau) ) {
      lambda <- exp(-log(2) / tau)
      i <- (0:(length(dates)-1))
      wt_tmp <- lambda^i
      wt[ ,1] <- rev( wt_tmp / sum(wt_tmp) * length(wt_tmp) )
    } 
    
    A_mat <- 0
    a_vec <- 0
    for ( today in dates ) {
      A_mat <- A_mat + as.numeric(wt[today, ]) * as.numeric(X[today, ]) %*% t(as.numeric(X[today, ]))
      a_vec <- a_vec + as.numeric(wt[today, ]) * 2 * as.numeric(X[today, ]) * as.numeric(y[today, ])
    }
    # Check
    # t(w) %*% A_mat %*% w - a_vec %*% w + sum(y^2) # same same
    # t(w) %*% A_mat %*% w; t(w) %*% (t(X) %*% X) %*% w
    # a_vec %*% w; 2 * t(y) %*% X %*% w
    # t(y) %*% y; sum(y^2)
    # cov(X) * length(dates)
    
    
    GPS@solver$portfolio <- "meanvariance"
    GPS@solver$obj_quad <- A_mat
    GPS@solver$obj_lin <- a_vec
    GPS@solver$obj_constant <- sum(y^2)
    
    GPO <- GPO::gpo( GPS = GPS )
    RBE$GPS <- GPS
    RBE$GPO <- GPO
    wghts <- getWeights(GPO)
    wghts <- timeSeries( matrix(wghts, nrow = 1, 
                                dimnames = list(rebdate, names(wghts)) ),
                         rebdate )
    lWeights <- list( weights = wghts )
    .self$appendOutput( lWeights )
    
    y_prime <- timeSeries( matrix(X %*% getWeights(GPO), ncol = 1), dates )
    Y <- cbind(y, prime = y_prime)
    plotSimTS( as.simTS(Y) )
    
    
    return( TRUE )
  }
  BacktestIndexReplication$methods( nnlsPortfolio = BacktestIndexReplication.nnlsPortfolio )
  
  
  
  
  # --------------------------------------------------------------------------
  BacktestIndexReplication.olsPortfolio <- function( rebdate = NULL )
  {
    
    RBE <- rebalenv( rebdate = rebdate )
    GPS <- RBE$GPS
    X_w <- getData( GPS )
    X_d <- window(data$X_sim[ ,colnames(X_w)], start(X_w), end(X_w))  # Use daily returns  -- ccy?
    X_d <- X_d[isWeekday(time(X_d)), ]
    X_tmp <- t( na.omit( t(X_d) ) )
    X <- timeSeries( X_tmp, rownames(X_d) )
    dates <- intersect( rownames(data$X_bm), rownames(X) )
    y <- log(1 + data$X_bm[dates, ])
    X <- log(1 + X[dates, ])
    colnames(y) <- "RV"
    colnames(X) <- paste0("EV_", colnames(X))
    
    # ###
    # y <- cumulated(y, "discrete")
    # X <- cumulated(X, "discrete")
    # ###
    
    # OLS Regression 
    DF <- as.data.frame(cbind(RV = y, X))
    f <- as.formula( paste0("RV ~ 0 + ", paste0(colnames(X), collapse = " + ")) )
    reg <- lm( formula = f, data = DF )
    coeff <- summary(reg)$coefficients
    wghts <- setNames( as.numeric(tail(coeff[ ,1], ncol(X))), colnames(X_tmp) )

    wghts <- timeSeries( matrix(wghts, nrow = 1, 
                                dimnames = list(rebdate, names(wghts)) ),
                         rebdate )
    lWeights <- list( weights = wghts )
    .self$appendOutput( lWeights )
    
    y_prime <- timeSeries( matrix(X %*% as.numeric(wghts), ncol = 1), dates )
    Y <- cbind(y, prime = y_prime)
    # plotSimTS( as.simTS(Y) )
    
    return( TRUE )
  }
  BacktestIndexReplication$methods( olsPortfolio = BacktestIndexReplication.olsPortfolio )
  
  
  
  
  
  # --------------------------------------------------------------------------
  BacktestIndexReplication.lsPortfolio <- function( rebdate = NULL )
  {
    
    RBE <- rebalenv( rebdate = rebdate )
    GPS <- RBE$GPS
    X_w <- getData( GPS )
    X_d <- window(data$X_sim[ ,colnames(X_w)], start(X_w), end(X_w))  # Use daily returns  -- ccy?
    X_d <- X_d[isWeekday(time(X_d)), ]
    X_tmp <- t( na.omit( t(X_d) ) )
    X <- timeSeries( X_tmp, rownames(X_d) )
    dates <- intersect( rownames(data$X_bm), rownames(X) )
    y <- log(1 + data$X_bm[dates, ])
    X <- log(1 + X[dates, ])
    colnames(y) <- "RV"
    colnames(X) <- paste0("EV_", colnames(X))
    
  
    
    # Regression
    # reg <- regression( Y_train = y,
    #                    X_train = X,
    #                    type = "elnet" )
    # beta_elnet <- reg$coeffmat[-1, 1]
    # wghts <- setNames( beta_elnet, colnames(X) )
    
    # DF <- as.data.frame(cbind(RV = y, X))
    # colnames(DF) <- c("RV", colnames(X))
    # f <- as.formula( paste0("RV ~ ", paste0(colnames(X), collapse = " + ")) )
    # reg <- lm( formula = f, data = DF )
    # coeff <- summary(reg)$coefficients
    # wghts <- coeff[ ,1]
    # y_prime <- timeSeries( X %*% wghts[-1], time(X) )
    # 
    # f2 <- as.formula( paste0("RV ~ 0 + ", paste0(colnames(X), collapse = " + ")) )
    # reg2 <- lm( formula = f2, data = DF )
    # coeff2 <- summary(reg2)$coefficients
    # wghts2 <- coeff2[ ,1]
    # y_prime_nointercept <- timeSeries( X %*% wghts2, time(X) )
    # 
    # y_prime3 <- timeSeries( cbind(X[ ,1] * 0 + 1, X) %*% wghts, time(X) )
    #   
    #   
    # Y <- cbind(y, prime = y_prime, prime_noint = y_prime_nointercept, prime3 = y_prime3)
    # plotSimTS( as.simTS(Y) )
    
    
    # DF <- as.data.frame(cbind(RV = y, X))
    # f <- as.formula( paste0("RV ~ 0 + ", paste0(colnames(X), collapse = " + ")) )
    # reg <- lm( formula = f, data = DF )
    # coeff <- summary(reg)$coefficients
    # wghts <- setNames( as.numeric(coeff[ ,1]), colnames(X_tmp) )
    
    
    # Elastic Net regression
    elnet_alpha <- 1
    # elnet_weights <- rep(1/nrow(X), nrow(X))
    tau <- 52
    lambda <- exp(-log(2) / tau)
    n <- dim(X)[1]
    i <- (0:(n-1))
    wt <- lambda^i
    elnet_weights <- rev( length(wt) * wt / sum(wt) )
    reg <- glmnet::cv.glmnet( y = as.matrix(y), 
                              x = as.matrix(X), 
                              # weights = elnet_weights,
                              intercept = FALSE,
                              alpha = elnet_alpha, 
                              family = "gaussian", 
                              type.measure = "mse" )
    y_prime <- timeSeries( predict( reg,
                                    s = reg$lambda.min,
                                    newx = as.matrix(X) ),
                           rownames(X) )
    coeffmat <- coef(reg, reg$lambda.min)
    wghts <- setNames( as.numeric(coeffmat[-1, 1]), colnames(X_tmp) )
    
    Y <- cbind(y, prime = y_prime)
    plotSimTS(as.simTS(Y))
    
    
    # Scale weights to sum to one
    if ( mean(wghts) > 0 ) {
      wghts <- wghts / sum(wghts)
    }
    
    # beta <- pseudoinverse( t(X) %*% X ) %*% t(X) %*% y
    # wghts <- setNames( beta, colnames(X) )
    
   
    wghts <- timeSeries( matrix(wghts, nrow = 1, 
                                dimnames = list(rebdate, names(wghts)) ),
                         rebdate )
    lWeights <- list( weights = wghts )
    .self$appendOutput( lWeights )
    
    return( TRUE )
  }
  BacktestIndexReplication$methods( lsPortfolio = BacktestIndexReplication.lsPortfolio )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  BacktestIndexReplication.elnetPortfolio <- function( rebdate = NULL )
  {
    
    RBE <- rebalenv( rebdate = rebdate )
    GPS <- RBE$GPS
    X_w <- getData( GPS )
    X_d <- window(data$X_sim[ ,colnames(X_w)], start(X_w), end(X_w))  # Use daily returns  -- ccy?
    X_d <- X_d[isWeekday(time(X_d)), ]
    X_tmp <- t( na.omit( t(X_d) ) )
    X <- timeSeries( X_tmp, rownames(X_d) )
    dates <- intersect( rownames(data$X_bm), rownames(X) )
    y <- log(1 + data$X_bm[dates, ])
    X <- log(1 + X[dates, ])
    colnames(y) <- "RV"
    colnames(X) <- paste0("EV_", colnames(X))
    
    ###
    y <- cumulated(y)
    X <- cumulated(X)
    ###
    
    # Elastic Net regression
    elnet_alpha <- spec$elnet_alpha
    # Weights
    tau <- spec$tau
    wt <- y * 0 + 1
    if ( !is.null(tau) ) {
      lambda <- exp(-log(2) / tau)
      i <- (0:(length(dates)-1))
      wt_tmp <- lambda^i
      wt[ ,1] <- rev( wt_tmp / sum(wt_tmp) * length(wt_tmp) )
    } 
    reg <- glmnet::cv.glmnet( y = as.matrix(y), 
                              x = as.matrix(X), 
                              weights = wt,
                              intercept = FALSE,
                              alpha = elnet_alpha, 
                              family = "gaussian", 
                              type.measure = "mse" )
    y_prime <- timeSeries( predict( reg,
                                    s = reg$lambda.min,
                                    newx = as.matrix(X) ),
                           rownames(X) )
    coeffmat <- coef(reg, reg$lambda.min)
    wghts <- setNames( as.numeric(coeffmat[-1, 1]), colnames(X_tmp) )
    
    # Y <- exp( cbind(y, prime = y_prime) ) - 1
    # plotSimTS(as.simTS(Y))
    Y <- cbind(y, y_prime)
    plot( Y, plot.type = "single" )
    
    
    # Scale weights to sum to one
    if ( mean(wghts) > 0 ) {
      wghts <- wghts / sum(wghts)
    }
    
    
    wghts <- timeSeries( matrix(wghts, nrow = 1, 
                                dimnames = list(rebdate, names(wghts)) ),
                         rebdate )
    lWeights <- list( weights = wghts )
    .self$appendOutput( lWeights )
    
    return( TRUE )
    
    
    
    
  }
  BacktestIndexReplication$methods( elnetPortfolio = BacktestIndexReplication.elnetPortfolio )
  
  
  
  
  # --------------------------------------------------------------------------
  BacktestIndexReplication.selection.100 <- function( rebdate = NULL )
  {
    RBE <-  rebalenv( rebdate = rebdate )
    RBE$selection <- RBE$selection[ 1:min(100, length(RBE$selection)) ]
    return( TRUE )
  }
  BacktestIndexReplication$methods( selection.100 = BacktestIndexReplication.selection.100 )
  
  
  