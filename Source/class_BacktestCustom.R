
  
  ############################################################################
  ### GEOMSCALE LOW VOLATILITY ANAOMALY - CUSTOM FUNCTIONS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard   19.05.2021
  # First version:    19.05.2021
  # --------------------------------------------------------------------------

  
  # BacktestCustom.setData
  # BacktestCustom.gpo
  # BacktestCustom.dirichletPortfolio
  # BacktestCustom.accrejPortfolio
  # BacktestCustom.decilePortfolio
  # BacktestCustom.quintilePortfolio
  # BacktestCustom.quantilePortfolio
  # BacktestCustom.scorePortfolio
  # BacktestCustom.score_rpPortfolio
  # BacktestCustom.momentum_rpPortfolio
  
  
  
 
  

  require(Backtest)
  
  
  
  # --------------------------------------------------------------------------
  BacktestCustom <- setRefClass( Class = "BacktestCustom", 
                                 contains = "BacktestPanel",
                                 methods = list() )
  
 
  
  # --------------------------------------------------------------------------
  BacktestCustom.setData <- function( rebdate = NULL ) 
  {
    
    filter_th <- spec$sharpe_filter_th
    filter_by <- ifelse( is.null(spec$filter_by), "sharpe", spec$filter_by )
    
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
    
    # Sort by Sharpe
    if ( !is.null(filter_th) ) {
      if ( filter_by == "sharpe" ) {
        # mu <- meanGeo( tail(X_est, 52), scalefactor = 1 )
        # mu <- exp( apply( log(1 + tail(X_est, 52)), 2, mean ) ) - 1
        mu <- apply( tail(X_est, 52), 2, mean )
        sds <- apply( X_est, 2, sd )
        score <- mu / sds
      } else if ( filter_by == "mean" ) {
        score <- apply( tail(X_est, 52), 2, mean )
      }
      
      id <- colnames(X_est)[ which( score > quantile(score, filter_th) ) ]
      X_est <- X_est[ ,id]  
    }

    setData(RBE$GPS) <- X_est
    # .self$transformData(rebdate = rebdate)
    
    return( TRUE )
  }
  BacktestCustom$methods( setData = BacktestCustom.setData )
  
  
  # --------------------------------------------------------------------------
  BacktestCustom.gpo <- function( rebdate = NULL )
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
  BacktestCustom$methods( gpo = BacktestCustom.gpo )
  
  
  
  # --------------------------------------------------------------------------
  BacktestCustom.dirichletPortfolio <- function( rebdate = NULL )
  {
    
    RBE <- rebalenv( rebdate = rebdate )
    X <- getData( RBE$GPS )
    n_sim <- ifelse( is.null(spec$n_sim), 10^2, spec$n_sim )
    fev_bias <- spec$fev_bias
    if ( is.null(fev_bias) ) {
      fev_bias <- 10
    }
   
    # Sample uniformly from Dirichlet distribution 
    alpha <- rep(1, ncol(X))
    samples <- rdirichlet( n = n_sim, alpha = alpha )
    colnames(samples) <- colnames(X)
    
    # FEV-Bias
    samples <- samples[rep(1:nrow(samples), (length(fev_bias)+1)), ]
    for ( i in seq(along = fev_bias) ) {
      idx <- ((n_sim*i)+1):(n_sim*(i+1))
      samples[idx, ]  <- fevBias( samples[1:n_sim, ], q = fev_bias[i] )
    }
    
    FUN <- function(w) { timeSeries( matrix( w, nrow = 1, dimnames = list(NULL, names(w)) ),
                                     rebdate ) }
    lWeights <- list()
    for ( k in 1:nrow(samples)) {
      lWeights[[k]] <- FUN( w = samples[k, ] )
    }
    names(lWeights) <- paste0("dirichlet", 1:length(lWeights))
   
    .self$appendOutput( lWeights )
  }
  BacktestCustom$methods( dirichletPortfolio = BacktestCustom.dirichletPortfolio )
  
  
  
  
  # --------------------------------------------------------------------------
  BacktestCustom.accrejPortfolio <- function( rebdate = NULL )
  {
    
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    n_sim <- ifelse( is.null(spec$n_sim), 10^2, spec$n_sim )
    fev_bias <- spec$fev_bias
    if ( is.null(fev_bias) ) {
      fev_bias <- 10
    }
    th <- as.numeric( spec$th[rebdate, ] )
    Score <- data$Score
    scores_vec <- setNames( as.numeric( Score[rebdate, selection] ), selection )
    scores_vec[is.null(scores_vec)] <- 0
    
    # Sample uniformly from Dirichlet distribution 
    alpha <- rep(1, length(selection))
    samples <- rdirichlet( n = n_sim, alpha = alpha )
    colnames(samples) <- selection
    
    # FEV-Bias
    samples <- samples[rep(1:nrow(samples), (length(fev_bias)+1)), ]
    for ( i in seq(along = fev_bias) ) {
      idx <- ((n_sim*i)+1):(n_sim*(i+1))
      samples[idx, ]  <- fevBias( samples[1:n_sim, ], q = fev_bias[i] )
    }
    
    # Acceptance - rejection step
    portfolio_scores <- apply( samples, 1, function(w) { t(w) %*% scores_vec }  )
    acc_idx <- matrix( portfolio_scores > th, nrow = 1, dimnames = list(rebdate, NULL) )
    acc_idx <- timeSeries( acc_idx, rebdate )
    lIdx <- list( acc_idx = acc_idx )
    
    FUN <- function(w) { timeSeries( matrix( w, nrow = 1, 
                                             dimnames = list(NULL, paste0("dirichlet", 1:length(selection))) ),
                                     rebdate ) }
    lWeights <- list()
    for ( k in 1:nrow(samples)) {
      lWeights[[k]] <- FUN( w = samples[k, ] )
    }
    names(lWeights) <- paste0("dirichlet", 1:length(lWeights))
    
    .self$appendOutput( lWeights )
    .self$appendOutput( lIdx )
  }
  BacktestCustom$methods( accrejPortfolio = BacktestCustom.accrejPortfolio )
  
    
  
  
  # --------------------------------------------------------------------------
  BacktestCustom.decilePortfolio <- function( rebdate = NULL )
  {
    
    spec$q_vec <<- seq(from = 0.1, to = 1, by = 0.1)
    .self$quantilePortfolio( rebdate = rebdate )
    
    return( TRUE )
  }
  BacktestCustom$methods( decilePortfolio = BacktestCustom.decilePortfolio )
  
  
  # --------------------------------------------------------------------------
  BacktestCustom.quintilePortfolio <- function( rebdate = NULL )
  {
    
    spec$q_vec <<- seq(from = 0.2, to = 1, by = 0.2)
    .self$quantilePortfolio( rebdate = rebdate )
    
    return( TRUE )
  }
  BacktestCustom$methods( quintilePortfolio = BacktestCustom.quintilePortfolio )
  
  
  # --------------------------------------------------------------------------
  BacktestCustom.quantilePortfolio <- function( rebdate = NULL )
  {
    
    q_vec <- spec$q_vec
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    Score <- data$Score
    wmat_bm <- data$wmat_bm
    weighting_method <- spec$weighting_method
    
    scores_vec <- setNames( as.numeric( Score[rebdate, selection] ), selection )
    scores <- get_sequence_of_scores( scores = scores_vec, m = length(q_vec) )
    lIdx <- attr(scores, "lIdx")
    weightingFUN <- function( idx, weighting_method = c("eqw", "capw") )
    { 
      id <- selection[idx]
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
    
    scores_ts <- timeSeries( matrix(scores, nrow = 1, 
                                    dimnames = list( rebdate, 
                                                     paste0("q", q_vec) ) ) )
    .self$appendOutput( list( scores = scores_ts ) )
    .self$appendOutput( lWeights )
    
    return( TRUE )
  }
  BacktestCustom$methods( quantilePortfolio = BacktestCustom.quantilePortfolio )
  
  
  
  # --------------------------------------------------------------------------
  BacktestCustom.scorePortfolio <- function( rebdate = NULL )
  {
    
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    Score <- data$Score
    scores_vec <- setNames( as.numeric( Score[rebdate, selection] ), selection )
    scores_vec[is.null(scores_vec)] <- 0
    wghts <- scores_vec / sum(scores_vec)
    wghts <- timeSeries( matrix(wghts, nrow = 1, dimnames = list(rebdate, selection) ),
                         rebdate )
    eps <- 1e-10
    wghts_inv <- (1 / (wghts + eps) ) / sum(1 / (wghts + eps))
    lWeights <- list( weights = wghts,
                      weights_inv = wghts_inv )
    .self$appendOutput( lWeights )
    return( TRUE )
  }
  BacktestCustom$methods( scorePortfolio = BacktestCustom.scorePortfolio )
  
  
  
  
  # --------------------------------------------------------------------------
  BacktestCustom.score_rpPortfolio <- function( rebdate = NULL )
  {
    #
    .self$scorePortfolio( rebdate = rebdate )
    #
    # Parameters
    RBE <- rebalenv( rebdate = rebdate )
    n_sim <- ifelse( is.null(spec$n_sim), 10^2, spec$n_sim )
    Names <- paste0("dirichlet", 1:n_sim)
    sampling_method <- ifelse( is.null(spec$sampling_method), "moon", spec$sampling_method )
    sampling_dist <- ifelse( is.null(spec$sampling_dist), "uniform", spec$sampling_dist )
    
    selection <- RBE$selection
    Score <- data$Score
    score_vec <- setNames( as.numeric( Score[rebdate, selection] ), selection )
    score_vec[is.na(score_vec)] <- 0
    m <- sum( score_vec > 0 )
    n <- length(selection)
    
    if ( sampling_dist == "uniform" ) {
      alpha <- rep(1/n, n)
    } else if ( sampling_dist == "capw" ) {
      .self$benchmarkWeights( rebdate = rebdate )
      w_bm <- getBenchmarkWeights(RBE$GPS)
      alpha <- w_bm / sum(w_bm)
    } else {
      stop("invalid sampling_dist.")
    }
    
    if ( sampling_method == "dirichlet" ) {
      
      # Sample uniformly from Dirichlet distribution
      # lambda <- mn2Lambda( n = n, m = m )
      # alpha_scl <- alpha * n * lambda    #// Not sure if this makes sense if alpha is asymmetric.
      alpha_scl <- alpha * (m - 1)         #// Not sure if this makes sense if alpha is asymmetric.
      samples <- rdirichlet( n = n_sim, alpha = alpha_scl )
      # # Rescale mini-positions
      # tol <- max( 0.003, 100 / m )
      # samples <- t( apply( samples, 1, function(x) { x[x < tol] <- 0; x <- x / sum(x) } ) )
      colnames(samples) <- selection
      
    } else if (sampling_method == "moon" ) {
      
      # m-out-of-n bootstrap
      samples <- matrix(0, nrow = n_sim, ncol = n, dimnames = list(NULL, selection))
      for ( i in 1:nrow(samples) ) {
        idx <- sample( 1:n, size = m, replace = TRUE, prob = alpha )
        tbl <- table(idx)
        samples[i, as.numeric(names(tbl))] <- tbl / m
      }
    } else {
      stop("invalid sampling_method.")
    }
    
    # Add capw to samples
    if ( sampling_dist == "capw" ) {
      samples <- rbind( alpha, samples )
      colnames(samples) <- selection
    }

    # Load initial weights and only append new random weights if there are is a change in the selection 
    # compared to the old allocation (otherwise buy and hold).
    idx <- which(spec$rebdates == rebdate)
    if (idx > 1) {
      lw_init <- list()
      for ( i in seq(along = Names) ) {
        yesterday <- max( rownames( output[[Names[[i]]]] )[ rownames( output[[Names[[i]]]] ) < rebdate ] )
        w_old <- output[[Names[i]]][yesterday, ]
        w_old <- w_old[ ,which(w_old > 0)]
        lw_init[[i]] <- na.omit( setNames( as.numeric(w_old), names(w_old)) )
      }
      tmp <- unlist( lapply( lw_init, FUN = function(w){ !all(names(w) %in% selection) } ) )
      not_ok <- which(tmp)
    } else {
      not_ok <- 1:n_sim
    }
    
    FUN <- function(w) { timeSeries( matrix( w, nrow = 1, dimnames = list(NULL, names(w)) ),
                                     rebdate ) }
    lWeights <- lapply( 1:nrow(samples), FUN = function(k) { FUN( samples[k, ] ) } )
    names(lWeights) <- Names
    if ( length(not_ok) > 0 ) {
      .self$appendOutput( lWeights[ not_ok ] )
    }
    
    return( TRUE )
  }
  BacktestCustom$methods( score_rpPortfolio = BacktestCustom.score_rpPortfolio )
  
  
  
  
  # --------------------------------------------------------------------------
  BacktestCustom.momentum_rpPortfolio <- function( rebdate = NULL )
  {
   
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    n <- length(selection)
    n_sim <- ifelse( is.null(spec$n_sim), 10^2, spec$n_sim )
    m <- ifelse( is.null(spec$m), length(selection), spec$m )
    th <- ifelse( is.null(spec$th), 0.05, spec$th )
    Names <- paste0("dirichlet", 1:n_sim)
    shadow_dirichlet_transform <- spec$shadow_dirichlet_transform
    scl_by_capw <- spec$scl_by_capw
    # sampling_method <- ifelse( is.null(spec$sampling_method), "moon", spec$sampling_method )
    sampling_dist <- ifelse( is.null(spec$sampling_dist), "uniform", spec$sampling_dist )

    # Get benchmark weights
    .self$benchmarkWeights( rebdate = rebdate )
    w_bm <- getBenchmarkWeights(RBE$GPS)
    w_bm <- w_bm[selection]
    
    # Sampling probabilities
    if ( sampling_dist == "uniform" ) {
      alpha <- rep(1/n, n)
    } else if ( sampling_dist == "capw" ) {
      alpha <- w_bm / sum(w_bm)
    } else {
      stop("invalid sampling_dist.")
    }

    # m-out-of-n bootstrap
    samples <- matrix(0, nrow = n_sim, ncol = n, dimnames = list(NULL, selection))
    for ( i in 1:nrow(samples) ) {
      idx <- sample( 1:n, size = m, replace = TRUE, prob = alpha )
      tbl <- table(idx)
      # samples[i, as.numeric(names(tbl))] <- tbl / m
      # Scale by cap-weight
      id <- as.numeric(names(tbl))
      w <- tbl / m
      if ( isTRUE(scl_by_capw) ) {
        alpha_scl <- w_bm[id] / sum(w_bm[id])
        w <- (w * alpha_scl) / sum( w * alpha_scl)
      }
      samples[i, id] <- w
    }
    
    # Shadow Dirichlet transformation
    if ( isTRUE(shadow_dirichlet_transform) ) {
      M <- diag( rep( th, n ) )
      row_sum <- apply( M, 2, sum )
      idx_row_sum <- which( 1 - row_sum > 0 )
      for ( i in seq(along = idx_row_sum) ) {
        M[-idx_row_sum[i], idx_row_sum[i]] <- (1 - row_sum[i]) / (n-1)
      }
      samples <- t( apply( samples, 1, function(w) { M %*% w } ) )
      colnames(samples) <- selection
    }
  
    
    # Load initial weights and only append new random weights if there are is a change in the selection
    # compared to the old allocation (otherwise buy and hold).
    idx <- which(spec$rebdates == rebdate)
    if (idx > 1) {
      lw_init <- list()
      for ( i in seq(along = Names) ) {
        yesterday <- max( rownames( output[[Names[[i]]]] )[ rownames( output[[Names[[i]]]] ) < rebdate ] )
        w_old <- output[[Names[i]]][yesterday, ]
        w_old <- w_old[ ,which(w_old > 0)]
        lw_init[[i]] <- na.omit( setNames( as.numeric(w_old), names(w_old)) )
      }
      tmp <- unlist( lapply( lw_init, FUN = function(w){ !all(names(w) %in% selection) } ) )
      not_ok <- which(tmp)
    } else {
      not_ok <- 1:n_sim
    }
    # not_ok <- 1:n_sim
    
    FUN <- function(w) { timeSeries( matrix( w, nrow = 1, dimnames = list(NULL, names(w)) ),
                                     rebdate ) }
    lWeights <- lapply( 1:nrow(samples), FUN = function(k) { FUN( samples[k, ] ) } )
    names(lWeights) <- Names
    if ( length(not_ok) > 0 ) {
      .self$appendOutput( lWeights[ not_ok ] )
    }
    
    return( TRUE )
  }
  BacktestCustom$methods( momentum_rpPortfolio = BacktestCustom.momentum_rpPortfolio )
  
  
  
  

  
  # --------------------------------------------------------------------------
  BacktestCustom.getInitialWeights <- function( rebdate  = NULL )
  {
    w_init <- NULL
    RBE <- rebalenv(rebdate = rebdate)
    idx <- which(spec$rebdates == rebdate)
    if (idx > 1) {
      dates <- rownames(output$weights)
      if (!is.null(dates)) {
        yesterday <- tail(dates[dates < rebdate], 1)
        old_names <- colnames(output$weights[yesterday, which(output$weights[yesterday, 
        ] != 0)])
        if (length(old_names) == 0) {
          old_names <- NULL
        }
      }
      else {
        old_names <- NULL
      }
      if (is.null(old_names)) {
        if (isTRUE(spec$verbose)) {
          warning("No weights from last rebalancing available that could be floated.")
        }
      }
      else {
        if (any(old_names %in% colnames(data$X_sim) == FALSE)) {
          stop("X_sim does not contain all asset names from last rebalancing.")
        }
        w_float <- floatingWeights(X = data$X_sim[, old_names], 
                                   w = output$weights[yesterday, old_names], startDate = yesterday, 
                                   endDate = timeDate(rebdate), rescale = TRUE)
        new_names <- getConstraints(RBE$GPS, "selection")
        same_names <- intersect(old_names, new_names)
        w_init <- setNames(numeric(length(new_names)), new_names)
        w_init[same_names] <- w_float[nrow(w_float), same_names]
      }
    }
    return( w_init )
  }
  BacktestCustom$methods( getInitialWeights = BacktestCustom.getInitialWeights )
  
  
            
            