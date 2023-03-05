
  
  ############################################################################
  ### RANDOM PORTFOLIOS - CLASS BACKTEST INDEX REPLICATION
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard   01.03.2023
  # First version:    01.03.2023
  # --------------------------------------------------------------------------

  
  # BacktestIndexReplication.setData
  # BacktestIndexReplication.gpo
  # BacktestIndexReplication.nnlsPortfolio
  # BacktestIndexReplication.lsPortfolio
  
  
  

  require(Backtest)
  
  
  
  
  # --------------------------------------------------------------------------
  BacktestIndexReplication <- setRefClass( Class = "BacktestIndexReplication", 
                                           contains = "BacktestPanel",
                                           methods = list() )
  
 
  
  
  
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
    X <- window(data$X_sim[ ,colnames(X_w)], start(X_w), end(X_w))  # Use daily returns  -- ccy?
    X <- X[isWeekday(time(X)), ]
    dates <- rownames(X)
    y <- data$X_bm[dates, ]
    ###
    # y <- abs(y)
    ###
    
    
    # Check
    w <- rep(1/ncol(X), ncol(X))
    ans <- y * NA
    ans2 <- ans
    for ( today in dates ) {
      a <- as.numeric(X[today, ])
      ans[today, ] <- (t(a) %*% w - y[today, ])^2
      A <- a %*% t(a)
      ans2[today, ] <- t(w) %*% A %*% w - t(2*y[today, ]*a) %*% w + y[today, ]^2
    }
    sum(ans)
    sum(ans2)  # same same
    sum( (X %*% w - y)^2 ) # same same
    t(w) %*% (t(X) %*% X) %*% w - t(y) %*% X %*% w - t(w) %*% t(X) %*% y - t(y) %*% y  # different ?
    t(w) %*% (t(X) %*% X) %*% w - 2 * t(y) %*% X %*% w - t(y) %*% y
    
    
    
    
    A_mat <- 0
    a_vec <- 0
    for ( today in dates ) {
      A_mat <- A_mat + as.numeric(X[today, ]) %*% t(as.numeric(X[today, ]))
      a_vec <- a_vec + 2 * as.numeric(X[today, ]) * as.numeric(y[today, ])
    }
    # Check
    # t(w) %*% A_mat %*% w - a_vec %*% w + sum(y^2) # same same
    # t(w) %*% A_mat %*% w; t(w) %*% (t(X) %*% X) %*% w
    # a_vec %*% w; 2 * t(y) %*% X %*% w
    # t(y) %*% y; sum(y^2)
    
    
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
    
    return( TRUE )
  }
  BacktestIndexReplication$methods( nnlsPortfolio = BacktestIndexReplication.nnlsPortfolio )
  
  
  # --------------------------------------------------------------------------
  BacktestIndexReplication.lsPortfolio <- function( rebdate = NULL )
  {
    
    RBE <- rebalenv( rebdate = rebdate )
    GPS <- RBE$GPS
    X <- getData( GPS )
    ###
    X <- scale(X, TRUE, FALSE)
    ###
    dates <- rownames(X)
    y <- data$X_bm[dates, ]
    y <- abs(y)
    colnames(y) <- "RV"
    # colnames(X) <- paste0("EV", 1:ncol(X))

    
    # Regression
    reg <- regression( Y_train = y,
                       X_train = X,
                       type = "elnet" )
    beta_elnet <- reg$coeffmat[-1, 1]
    wghts <- setNames( beta_elnet, colnames(X) )
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
  