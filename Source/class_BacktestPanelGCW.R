  
  
  ############################################################################
  ### GEOMSCALE LOW VOLATILITY ANAOMALY - CUSTOM BACKTEST CLASS - GREAT CYCLE WALK
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard   19.05.2021
  # First version:    19.05.2021
  # --------------------------------------------------------------------------
  
  
  
  # BacktestPanelGCW.setData
  # BacktestPanelGCW.gpo
  # BacktestPanelGCW.quantileSimplexPortfolio
  # BacktestPanelGCW.gcwPortfolio
  # BacktestPanelGCW.dirichletPortfolio  
  # BacktestPanelGCW.selection.subset
  
  
  require(Backtest)
  
  
  
  # --------------------------------------------------------------------------
  BacktestPanelGCW <- setRefClass( Class = "BacktestPanelGCW", 
                                   contains = "BacktestPanel",
                                   methods = list() )
  
  
  # # --------------------------------------------------------------------------
  # BacktestPanelGCW.gps <- function (rebdate = NULL) 
  # {
  #   .self$setGPS(rebdate = rebdate)
  #   .self$setSelection(rebdate = rebdate)
  #   .self$setData(rebdate = rebdate)
  #   .self$setCovariance(rebdate = rebdate)
  #   .self$setConstraints(rebdate = rebdate)
  #   return(TRUE)
  # }
  # BacktestPanelGCW$methods( gps = BacktestPanelGCW.gps )
  
 
  # --------------------------------------------------------------------------
  BacktestPanelGCW.setData <- function( rebdate = NULL ) 
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
  BacktestPanelGCW$methods( setData = BacktestPanelGCW.setData )
  
  
  # --------------------------------------------------------------------------
  BacktestPanelGCW.gpo <- function( rebdate = NULL )
  {
    txt <- paste0(".self$", spec$portfolio, "Portfolio( rebdate = '", rebdate, "')" )
    eval( parse( text = txt ) )
    return( TRUE )
  }
  BacktestPanelGCW$methods( gpo = BacktestPanelGCW.gpo )
  
  
  # --------------------------------------------------------------------------
  BacktestPanelGCW.quantileSimplexPortfolio <- function( rebdate = NULL )
  {
    
    q_vec <- spec$q_vec
    q_quant <- spec$quant
    n_sim <- ifelse( is.null(spec$n_sim), 10^2, spec$n_sim )
    RBE <- rebalenv( rebdate = rebdate )
    X <- getData(RBE$GPS)
    
    if ( q_quant <= length(q_vec) ) {
      
      # Sort by volatility
      sds <- apply( X, 2, sd )
      th <- quantile(sds, q_vec)
      lID <- list()
      for ( i in seq(along = th) ) {
        if ( i == 1 ) {
          lID[[i]] <- names(sds)[ which(sds <= th[i]) ]
        } else {
          lID[[i]] <- names(sds)[ intersect( which(sds > th[i-1]), which(sds <= th[i]) ) ]
        }
      }
      
      # Sample from quantile-buckets
      id <- lID[[ spec$quant ]]
      S <- Simplex$new( d = length(id) )
      samples <- S$runif( n_sim )
      colnames(samples) <- id
      
    } else {
      
      # Sample from entire dataset
      id <- colnames(X)
      S <- Simplex$new( d = length(id) )
      samples <- S$runif( n_sim )
      colnames(samples) <- id
    }
    
    # FEV-Bias
    samples_fev <- fevBias( samples, q = 10 )
    
    
    FUN <- function(w) { timeSeries( matrix( w, nrow = 1, dimnames = list(NULL, names(w)) ), rebdate ) }
    lWeights <- list()
    lWeights_fev <- list()
    for ( i in 1:nrow(samples)) {
      lWeights[[i]] <- FUN( w = samples[i, ] )
      lWeights_fev[[i]] <- FUN( w = samples_fev[i, ] )
    }
    names(lWeights) <- paste0("gcw", 1:length(lWeights))
    names(lWeights_fev) <- paste0("gcw_fev", 1:length(lWeights))
    
    .self$appendOutput( lWeights )
    .self$appendOutput( lWeights_fev )
    
    return( TRUE )
  }
  BacktestPanelGCW$methods( quantileSimplexPortfolio = BacktestPanelGCW.quantileSimplexPortfolio )
  
  
  
  # --------------------------------------------------------------------------
  BacktestPanelGCW.gcwPortfolio <- function( rebdate = NULL )
  {
    
    n_volas <- spec$n_volas
    vola_idx <- spec$vola_idx
    n_sim <- ifelse( is.null(spec$n_sim), 10^2, spec$n_sim )
    
    # Load data and compute covariance
    RBE <- rebalenv( rebdate = rebdate )
    X <- getData(RBE$GPS)
    covmat <- covariance(RBE$GPS)
    
    # Variance target
    sigmas <- get_sequence_of_volatilities_CB( sigma = covmat, 
                                               m = n_volas,
                                               ignore_cov = FALSE )
    sigma_target <- sigmas[ vola_idx ]
    
    # Sampling
    # parameters <- get_parameters( ncol(X)-1 )
    lSamples <- try( sample_ptfs_constant_volatility( sigma = covmat, 
                                                      c = sigma_target, 
                                                      M = n_sim,
                                                      ignore_smallest_components = FALSE ) )
    if ( !inherits(lSamples, "try-error") ) {
        
      output$n_components <<- c(output$n_components, length(lSamples[[1]]))
      
      samples <- t( do.call( cbind, lSamples[[1]] ) )
      if ( nrow(samples) > n_sim ) {
        idx <- sample( x = 1:nrow(samples), size = n_sim )
        samples <- samples[idx, ]
      }
      colnames(samples) <- colnames(X)
      
      # Cast samples to weights
      FUN <- function(w) { timeSeries( matrix( w, nrow = 1, 
                                               dimnames = list(NULL, names(w)) ), 
                                       rebdate ) }
      lWeights <- list()
      for ( i in 1:nrow(samples)) {
        lWeights[[i]] <- FUN( w = samples[i, ] )
      }
      names(lWeights) <- paste0("gcw", 1:length(lWeights))
      
      # Append weights
      .self$appendOutput( lWeights )
      
    } else {
      if ( isTRUE(spec$skip_error) ) {
        out <- list( attr(samples, "condition") )
        names(out) <- rebdate
        output$error_log <<- c(output$error_log, out)
        return( TRUE )
      } else {
        return( FALSE )
      }
    }
  }
  BacktestPanelGCW$methods( gcwPortfolio = BacktestPanelGCW.gcwPortfolio )
  
  
  
  # --------------------------------------------------------------------------
  BacktestPanelGCW.dirichletPortfolio <- function( rebdate = NULL )
  {

    n_sim <- ifelse( is.null(spec$n_sim), 10^2, spec$n_sim )
    k <- ifelse( is.null(spec$k), 0, spec$k )
    b_inverse_vola <- spec$b_inverse_vola
    RBE <- rebalenv( rebdate = rebdate )
    X <- getData(RBE$GPS)
    sds <- apply( X, 2, sd )
    
    
    # Sample from Dirichlet distribution with parameter proportional to asset volas
    if ( isTRUE(b_inverse_vola) ) {
      alpha <- (1 / sds^k) / sum(1 / sds^k)
    } else {
      alpha <- sds^k / sum(sds^k)
    }
    samples <- rdirichlet( n = n_sim, alpha = alpha * length(alpha) )
    colnames(samples) <- colnames(X)
  
    # FEV-Bias
    samples_fev <- fevBias( samples, q = 10 )
    
    FUN <- function(w) { timeSeries( matrix( w, nrow = 1, dimnames = list(NULL, names(w)) ), rebdate ) }
    lWeights <- list()
    lWeights_fev <- list()
    for ( k in 1:nrow(samples)) {
      lWeights[[k]] <- FUN( w = samples[k, ] )
      lWeights_fev[[k]] <- FUN( w = samples_fev[k, ] )
    }
    names(lWeights) <- paste0("dirichlet", 1:length(lWeights))
    names(lWeights_fev) <- paste0("dirichlet_fev", 1:length(lWeights))
    
    .self$appendOutput( lWeights )
    .self$appendOutput( lWeights_fev )
  }
  BacktestPanelGCW$methods( dirichletPortfolio = BacktestPanelGCW.dirichletPortfolio )
  
  
  
  # --------------------------------------------------------------------------
  BacktestPanelGCW.selection.subset <- function( rebdate = NULL )
  {
    quantile_subset <- spec$quantile_subset
    selection <- data$lSelection[[rebdate]]
    idx1 <- floor( length(selection) * (quantile_subset - 0.2) ) + 1
    idx2 <- floor( length(selection) * quantile_subset )
    if ( idx1 == 2 ) { idx1 <- 1 }
    selection <- selection[idx1:idx2]
    return( selection )
  }
  BacktestPanelGCW$methods( selection.subset = BacktestPanelGCW.selection.subset )
  
  
