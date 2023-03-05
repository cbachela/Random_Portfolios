
  
  ############################################################################
  ### GEOMSCALE LOW VOLATILITY ANAOMALY - CUSTOM FUNCTIONS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard   19.05.2021
  # First version:    19.05.2021
  # --------------------------------------------------------------------------

  
  # BacktestLVA.gps
  # BacktestLVA.gpo
  # BacktestLVA.setData
  # BacktestLVA.decilePortfolio
  # BacktestLVA.quintilePortfolio
  # BacktestLVA.quantilePortfolio
  # BacktestLVA.qcquantilePortfolio
  # BacktestLVA.mvPortfolio
  # BacktestLVA.selection.subset
  # BacktestLVA.selection.size
  # BacktestLVA.selection.sector
  
  
 
  

  require(Backtest)
  
  
  
  # --------------------------------------------------------------------------
  BacktestLVA <- setRefClass( Class = "BacktestLVA", 
                              contains = "BacktestPanel",
                              methods = list() )
  
  # --------------------------------------------------------------------------
  BacktestLVA.gps <- function( rebdate = NULL )
  {
    .self$setGPS(rebdate = rebdate)
    .self$setSelection(rebdate = rebdate)
    .self$setData(rebdate = rebdate)
    .self$setCovariance(rebdate = rebdate)
    .self$setConstraints(rebdate = rebdate)
    return(TRUE)
  }
  BacktestLVA$methods( gps = BacktestLVA.gps )
  
  
  # --------------------------------------------------------------------------
  BacktestLVA.selection.size <- function( rebdate = NULL )
  {
    
    # Sort by size
    if ( !spec$sort_by == "none" ) {
      w_bm <- setNames( as.numeric(data$wmat_bm[rebdate, selection]), selection )
      if ( spec$sort_by == "large" ) {
        # Select largest stocks
        id <- names(rev(sort(w_bm)))
      } else if ( spec$sort_by == "small" ) {
        # Select smallest stocks
        id <- names(sort(w_bm)) 
      } else {
        stop( "sort_by argument matches no option." )
      }
      n_stocks <- ifelse( is.null(spec$n_stocks), 70, spec$n_stocks )
      if ( n_stocks < 1 ) {
        n_stocks <- floor(length(id) * n_stocks)
        ###
        n_max <- ifelse( is.null(spec$n_max), 150, spec$n_max)
        n_stocks <- min( n_stocks, n_max )   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        ###
      }
      id <- id[1:n_stocks]
      X_est <- X_est[ ,id]  
    }
    ###
  }  
  
  
  # --------------------------------------------------------------------------
  BacktestLVA.gpo <- function( rebdate = NULL )
  {
    txt <- paste0(".self$", spec$portfolio, "Portfolio( rebdate = '", rebdate, "')" )
    eval( parse( text = txt ) )
    return( TRUE )
  }
  BacktestLVA$methods( gpo = BacktestLVA.gpo )
  
  
  # --------------------------------------------------------------------------
  BacktestLVA.setData <- function( rebdate = NULL ) 
  {
    RBE <- rebalenv(rebdate = rebdate)
    if (is.null(RBE$selection)) {
      .self$setSelection(rebdate = rebdate)
    }
    selection <- RBE$selection
    lridx <- which(rownames(data$X_est) == rebdate)
    if (spec$width == 0 | spec$width > lridx) {
      fridx <- 1
    } else {
      fridx <- max(1, lridx - spec$width)
    }
    X_est <- data$X_est[fridx:lridx, selection]
    
    ###
    # Sort by size
    if ( !spec$sort_by == "none" ) {
      w_bm <- setNames( as.numeric(data$wmat_bm[rebdate, selection]), selection )
      if ( spec$sort_by == "large" ) {
        # Select largest stocks
        id <- names(rev(sort(w_bm)))
      } else if ( spec$sort_by == "small" ) {
        # Select smallest stocks
        id <- names(sort(w_bm)) 
      } else {
        stop( "sort_by argument matches no option." )
      }
      n_stocks <- ifelse( is.null(spec$n_stocks), 70, spec$n_stocks )
      if ( n_stocks < 1 ) {
        n_stocks <- floor(length(id) * n_stocks)
        ###
        n_max <- ifelse( is.null(spec$n_max), 150, spec$n_max)
        n_stocks <- min( n_stocks, n_max)   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        ###
      }
      id <- id[1:n_stocks]
      X_est <- X_est[ ,id]  
    }
    ###
    
    setData(RBE$GPS) <- X_est
    # .self$transformData(rebdate = rebdate)
    
    return( TRUE )
  }
  BacktestLVA$methods( setData = BacktestLVA.setData )
  
  
  # --------------------------------------------------------------------------
  BacktestLVA.decilePortfolio <- function( rebdate = NULL )
  {
    
    spec$q_vec <<- seq(from = 0.1, to = 1, by = 0.1)
    .self$quantilePortfolio( rebdate = rebdate )
   
    return( TRUE )
  }
  BacktestLVA$methods( decilePortfolio = BacktestLVA.decilePortfolio )
  
  
  # --------------------------------------------------------------------------
  BacktestLVA.quintilePortfolio <- function( rebdate = NULL )
  {
    
    spec$q_vec <<- seq(from = 0.2, to = 1, by = 0.2)
    .self$quantilePortfolio( rebdate = rebdate )
    
    return( TRUE )
  }
  BacktestLVA$methods( quintilePortfolio = BacktestLVA.quintilePortfolio )
  
  
  # --------------------------------------------------------------------------
  BacktestLVA.quantilePortfolio <- function( rebdate = NULL )
  {
    
    q_vec <- spec$q_vec
    RBE <- rebalenv( rebdate = rebdate )
    X <- getData(RBE$GPS)
    covmat <- covariance( RBE$GPS )
    wmat_bm <- data$wmat_bm
    weighting_method <- spec$weighting_method
    
    
    
    
    # # Sort by volatility
    # # sds <- apply( X, 2, sd )
    # covmat <- covariance( RBE$GPS )
    # sds <- sqrt(diag(covmat))
    # th <- quantile(sds, q_vec)
    # lID <- list()
    # for ( i in seq(along = th) ) {
    #   if ( i == 1 ) {
    #     lID[[i]] <- names(sds)[ which(sds <= th[i]) ]
    #   } else {
    #     lID[[i]] <- names(sds)[ intersect( which(sds > th[i-1]), which(sds <= th[i]) ) ]
    #   }
    # }
    # # lSD <- lapply( lID, FUN = function(id) { sds[id] } )
    # weightingFUN <- function( id, weighting_method = c("eqw", "capw") )
    # { 
    #   if ( weighting_method == "eqw" ) {
    #     wghts <- setNames( rep(1/length(id), length(id)), id ) 
    #   } else {
    #     wghts_bm <- setNames( as.numeric(wmat_bm[rebdate, id]), id )
    #     wghts <- wghts_bm / sum(wghts_bm)
    #   }
    # }
    # lW <- lapply( lID, FUN = weightingFUN, weighting_method = weighting_method )
    # FUN <- function(w) { timeSeries( matrix( w, nrow = 1, dimnames = list(NULL, names(w)) ), rebdate ) }
    # lWeights <- lapply( lW, FUN = FUN )
    # names(lWeights) <- paste0("q", q_vec)
    
    
    ### 
    n_volas <- length(q_vec)
    sigmas <- get_sequence_of_volatilities_CB( sigma = covmat,
                                               m = n_volas,
                                               ignore_cov = FALSE )   #// assumes eqw
    lIdx <- attr(sigmas, "lIdx")
    weightingFUN <- function( idx, weighting_method = c("eqw", "capw") )
    { 
      id <- colnames(covmat)[idx]
      if ( weighting_method == "eqw" ) {
        wghts <- setNames( rep(1/length(id), length(id)), id ) 
      } else {
        wghts_bm <- setNames( as.numeric(wmat_bm[rebdate, id]), id )
        wghts <- wghts_bm / sum(wghts_bm)
      }
    }
    lW2 <- lapply( lIdx, FUN = weightingFUN, weighting_method = weighting_method )
    FUN <- function(w) { timeSeries( matrix( w, nrow = 1, dimnames = list(NULL, names(w)) ), rebdate ) }
    lWeights <- lapply( lW2, FUN = FUN )
    names(lWeights) <- paste0("q", q_vec)
    
    
    # sigmas_check <- sigmas * NA
    # for ( i in 1:length(lW2) ) {
    #   w <- lW2[[i]]
    #   sigmas_check[i] <- t(w) %*% covmat[names(w), names(w)] %*% w
    # }
    # cbind( sigmas, sigmas_check )
    
    
    sigmas_ts <- timeSeries( matrix(sigmas, nrow = 1, 
                                    dimnames = list( rebdate, 
                                                     paste0("q", q_vec) ) ) )
    .self$appendOutput( list( sigmas = sigmas_ts ) )
    .self$appendOutput( lWeights )
  
    return( TRUE )
  }
  BacktestLVA$methods( quantilePortfolio = BacktestLVA.quantilePortfolio )
  
  
  
  # --------------------------------------------------------------------------
  BacktestLVA.qcquantilePortfolio <- function( rebdate = NULL )
  {
    
    RBE <- rebalenv( rebdate = rebdate )
    GPS <- RBE$GPS
    GPS@solver$recursion <- FALSE
    GPS@constraints@nonlinear <- list()
    X <- getData(GPS)
    q_vec <- spec$q_vec
    if ( is.null(q_vec) ) {
      q_vec <- seq(from = 0.2, to = 1, by = 0.2)
    }
    
    # Minimum variance portfolio
    gpo_mv <- minvariancePortfolio( GPS = GPS )
    w_mv <- getWeights(gpo_mv)
    
    # Max variance-score portfolio
    covmat <- covariance(GPS)
    gps_maxvar <- GPS
    gps_maxvar@solver$obj_lin <- (-1) * diag(covmat)
    gps_maxvar@solver$obj_quad <- covmat * 0
    gpo_maxvar <- linquadPortfolio( GPS = gps_maxvar )
    w_maxvar <- getWeights(gpo_maxvar)
    
    
    sigma_min <- as.numeric( w_mv %*% covmat %*% w_mv )
    sigma_max <- as.numeric( w_maxvar %*% covmat %*% w_maxvar )
    sigma_vec <- seq(from = sigma_min, to = sigma_max, length.out = 10)
    sigma_th <- quantile( sigma_vec, q_vec )
    
    GPS <- gps_maxvar
    # GPS@solver$portfolio <- "minturnover"
    GPS@solver$progtype <- "QCP"
    # setInitialWeights(GPS) <- w_mv * 0 + rep(1/length(w_mv), length(w_mv))
    
    lWeights <- list()
    for ( i in seq(along = sigma_th) ) {
      
      # Add quadratic constraint
      addConstraint(GPS) <- varianceConstraint( rhs = sigma_th[i],
                                                sense = "<=",
                                                Qmat = covmat )
      # GPO <- GPO::gpo( GPS = GPS )
      # debugonce( gpp.qcp.gurobi )
      GPO <- linquadPortfolio( GPS = GPS )
      w <- getWeights(GPO)
      lWeights[[i]] <- timeSeries( matrix( w, nrow = 1, dimnames = list(NULL, names(w)) ), rebdate )
    }
    names(lWeights) <- paste0("q", q_vec)
    
    .self$appendOutput( lWeights )
    return( TRUE )
  }
  BacktestLVA$methods( qcquantilePortfolio = BacktestLVA.qcquantilePortfolio )
  
  
  # --------------------------------------------------------------------------
  BacktestLVA.mvPortfolio <- function( rebdate = NULL )
  {
    tic <- Sys.time()
    .self$initialWeights(rebdate = rebdate)
    .self$floatUpperBounds(rebdate = rebdate)
    RBE <- rebalenv(rebdate = rebdate)
    RBE$GPO <- GPO::gpo(GPS = RBE$GPS)
    toc <- Sys.time() - tic
    wghts <- matrix(getWeights(RBE$GPO), nrow = 1, dimnames = list(rebdate, 
                                                                   names(getWeights(RBE$GPO))))
    tictoc <- matrix(toc, dimnames = list(rebdate, "tictoc"))
    .self$appendOutput(list(weights = as.timeSeries(wghts), runtime = as.timeSeries(tictoc)))
    if (isTRUE(spec$clean_RBE) && exists("RBE")) {
      cleanEnvir(keep = "rebdate", envir = RBE)
    }
    return(TRUE)
  }
  BacktestLVA$methods( mvPortfolio = BacktestLVA.mvPortfolio )
  
  
  # --------------------------------------------------------------------------
  BacktestLVA.selection.subset <- function( rebdate = NULL )
  {
    quantile_subset <- spec$quantile_subset
    selection <- data$lSelection[[rebdate]]
    idx1 <- floor( length(selection) * (quantile_subset - 0.2) ) + 1
    idx2 <- floor( length(selection) * quantile_subset )
    if ( idx1 == 2 ) { idx1 <- 1 }
    selection <- selection[idx1:idx2]
    return( selection )
  }
  BacktestLVA$methods( selection.subset = BacktestLVA.selection.subset )
  
  
  
  # --------------------------------------------------------------------------
  BacktestLVA.selection.sector <- function( rebdate = NULL )
  {
    quantile_subset <- spec$quantile_subset
    selection <- data$lSelection[[rebdate]]
    idx1 <- floor( length(selection) * (quantile_subset - 0.2) ) + 1
    idx2 <- floor( length(selection) * quantile_subset )
    if ( idx1 == 2 ) { idx1 <- 1 }
    selection <- selection[idx1:idx2]
    return( selection )
  }
  BacktestLVA$methods( selection.subset = BacktestLVA.selection.subset )
  
  
  
  