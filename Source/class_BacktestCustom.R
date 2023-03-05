
  
  ############################################################################
  ### RANDOM PORTFOLIOS - CUSTOM BACKTEST CLASS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard   19.05.2021
  # First version:    19.05.2021
  # --------------------------------------------------------------------------

  
  # BacktestCustom.setData
  # BacktestCustom.gpo
  # BacktestCustom.constraints.vola
  # BacktestCustom.dirichletPortfolio
  # BacktestCustom.accrejPortfolio
  # BacktestCustom.decilePortfolio
  # BacktestCustom.quintilePortfolio
  # BacktestCustom.quantilePortfolio
  # BacktestCustom.scorePortfolio
  # BacktestCustom.score_rpPortfolio
  # BacktestCustom.momentum_rpPortfolio
  # BacktestCustom.maxMomPortfolio
  # BacktestCustom.momScorePortfolio
  # BacktestCustom.volScorePortfolio
  # BacktestCustom.grwPortfolio
  # BacktestCustom.maxRCPortfolio
  # BacktestCustom.MSCIMultifactorConstraints
  
  
 
  

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
  BacktestCustom.constraints.vola <- function( rebdate = NULL )
  {
    vola_multiple <- ifelse( is.null(spec$vola_multiple), 1, spec$vola_multiple)
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    X <- data$X_est_d[which(rownames(data$X_est_d) <= rebdate), selection]
    X <- tail(X, 252)
    X[is.na(X)] <- 0
    sds <- apply( X, 2, sd )
    .self$benchmarkWeights( rebdate = rebdate )
    bm_weights <- getBenchmarkWeights( RBE$GPS )
    bm_weights <- bm_weights[ selection ] / sum( bm_weights[ selection ] )
    bm_score <- sum( bm_weights * sds ) * vola_multiple
    ans <- linearConstraint( name = "vola",
                             sense = "<=",
                             rhs = bm_score,
                             Amat = matrix( sds, nrow = 1, dimnames = list("vola", names(sds)) ) ) 
    return( ans )
  }
  BacktestCustom$methods( constraints.vola = BacktestCustom.constraints.vola )
  
    
  
  
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
    keep_all_rebdates <- spec$keep_all_rebdates
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
  
    # Load initial weights and only append new random weights if there are is a change 
    # in the selection compared to the old allocation (otherwise buy and hold).
    not_ok <- 1:n_sim
    if ( !isTRUE(keep_all_rebdates) ) {
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
      } 
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
  
  
            
  # --------------------------------------------------------------------------
  BacktestCustom.maxMomPortfolio <- function( rebdate  = NULL )
  {
    RBE <- rebalenv(rebdate = rebdate)
    GPS <- RBE$GPS
    X_w <- getData( GPS )
    X <- window( data$X_est_d[ ,colnames(X_w)], start(X_w), end(X_w) )
    X[is.na(X)] <- 0
    # mom <- exp( apply( log(1 + X_w[-1, ]), 2, mean ) ) - 1
    mom <- exp( apply( X[-1, ], 2, mean ) ) - 1 # estimation data is already continuous returns
    GPS@solver$obj_lin <- -mom
    GPS@solver$utility <- list(riskaversion = 0)
    GPS@solver$portfolio <- "meanvariance"
    GPO <- GPO::gpo(GPS)
    wghts <- GPO::getWeights(GPO)
    wghts <- timeSeries( matrix( wghts, nrow = 1, dimnames = list(NULL, names(wghts)) ),
                         rebdate )
    .self$appendOutput( value = list( weights = wghts ) )
    return( TRUE )
  }
  BacktestCustom$methods( maxMomPortfolio = BacktestCustom.maxMomPortfolio )
  
  # --------------------------------------------------------------------------
  BacktestCustom.momScorePortfolio <- function( rebdate  = NULL )
  {
    RBE <- rebalenv(rebdate = rebdate)
    X_w <- getData( RBE$GPS )
    X <- window( data$X_est_d[ ,colnames(X_w)], start(X_w), end(X_w) )
    X[is.na(X)] <- 0
    # mom <- momentum( Data = X, spec = momCtrl( method = "cumretEwma" ) )
    # mom <- exp( apply( log(1 + X[-1, ]), 2, mean ) ) - 1
    mom <- exp( apply( X[-1, ], 2, mean ) ) - 1   # estimation data is already continuous returns
    mom <- timeSeries( matrix( mom, nrow = 1, dimnames = list(NULL, names(mom)) ),
                       rebdate )
    .self$appendOutput( value = list( scores = mom ) )
    return( TRUE )
  }
  BacktestCustom$methods( momScorePortfolio = BacktestCustom.momScorePortfolio )
  
  # --------------------------------------------------------------------------
  BacktestCustom.volScorePortfolio <- function( rebdate  = NULL )
  {
    RBE <- rebalenv(rebdate = rebdate)
    X_w <- getData( RBE$GPS )
    X <- window( data$X_est_d[ ,colnames(X_w)], start(X_w), end(X_w) )
    X[is.na(X)] <- 0
    sds <- setNames( as.numeric( apply( X, 2, sd ) ), colnames(X) )
    sds <- timeSeries( matrix( sds, nrow = 1, dimnames = list(NULL, names(sds)) ),
                       rebdate )
    .self$appendOutput( value = list( scores = sds ) )
  }
  BacktestCustom$methods( volScorePortfolio = BacktestCustom.volScorePortfolio )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  BacktestCustom.grwPortfolio <- function( rebdate  = NULL )
  {
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    n_sim <- ifelse( is.null(spec$n_sim), 10^2, spec$n_sim )
    Names <- paste0("grw", 1:n_sim)
    keep_all_rebdates <- spec$keep_all_rebdates
    
    # Extract constraints
    Cons <- getConstraints( RBE$GPS )
    P <- lincon2Polytope( object = Cons )
    idx_eq <- which(P$sense == "=")
    A <- P$A[-idx_eq, ]
    b <- P$b[-idx_eq]
    Aeq <- matrix( P$A[idx_eq, ], nrow = length(idx_eq), ncol = ncol(P$A), byrow = TRUE )
    beq <- P$b[idx_eq]
    
    # Sample with volesti
    # pre_proc_list = preprocess_with_quadprog(A, b, Aeq, beq)
    # debugonce( preprocess_with_quadprog )
    result_list <- samples_uniform_portfolios( A = A, b = b, Aeq = Aeq, beq = beq, ess = n_sim )
    # samples <- t(result_list$samples)
    samples <- t(result_list$random_portfolios)
    colnames(samples) <- selection
    samples <- samples[ sample(1:nrow(samples))[1:n_sim], ]
    
    # Load initial weights and only append new random weights if there are is a change 
    # in the selection compared to the old allocation (otherwise buy and hold).
    not_ok <- 1:n_sim
    if ( !isTRUE(keep_all_rebdates) ) {
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
      } 
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
  BacktestCustom$methods( grwPortfolio = BacktestCustom.grwPortfolio )
  
  
  
  
  # --------------------------------------------------------------------------
  BacktestCustom.maxRCPortfolio <- function( rebdate  = NULL )
  {
    # Maximize random centroid portfolio 
    
    # Parameters
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    n_sim <- ifelse( is.null(spec$n_sim), 10^2, spec$n_sim )
    Names <- paste0("maxRC", 1:n_sim)
    keep_all_rebdates <- spec$keep_all_rebdates
    n_pos <- spec$n_pos
    
    # Compute centroid
    z <- ffv::centroid( n = length(selection) )
    
    # Shift centroid so that more values are negative 
    # (which means that less stocks will be in the portfolio if n_pos < 1/max_weight)
    if ( !is.null(n_pos) ) {
      z <- z - z[n_pos]
    }
    
    # Prepare linear optimization problem
    GPS <- RBE$GPS
    GPS@solver$progtype = "LP"
    GPS@solver$verbose <- FALSE
    GPS@solver$obj_lin <- -z
    GPS@solver$obj_quad <- NULL
    GPP <- gpp( GPS = GPS )
    
    # Sampling
    samples <- matrix( NA, nrow = n_sim, ncol = length(selection), 
                       dimnames = list( Names, selection) )
    for ( i in 1:n_sim ) {
      mu <- sample(z)
      GPP@model$model$obj <- mu
      GPO <- gpsolve(GPP)
      samples[i, ] <- getWeights(GPO)
    }
    
    # Load initial weights and only append new random weights if there are is a change 
    # in the selection compared to the old allocation (otherwise buy and hold).
    not_ok <- 1:n_sim
    if ( !isTRUE(keep_all_rebdates) ) {
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
      } 
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
  BacktestCustom$methods( maxRCPortfolio = BacktestCustom.maxRCPortfolio )
  
  
  
  
  
  # --------------------------------------------------------------------------
  BacktestCustom.MSCIMultifactorConstraints <- function( rebdate  = NULL )
  {
    
    wmat_bm <- data$wmat_bm
    if ( is.null(wmat_bm) ) { 
      stop( "Matrix of benchmark weights is missing.")
    }
    country_mat <- data$country_mat
    if ( is.null(country_mat) ) { 
      stop( "Matrix of country belongingness is missing.") 
    }
    sector_mat <- data$sector_mat
    if ( is.null(sector_mat) ) { 
      stop( "Matrix of sectort belongingness is missing.") 
    }
    
    # Keep raw data
    if ( !is.null(data$upper_mat) ) data$upper_mat_raw <<- data$upper_mat
    if ( !is.null(data$country_UB) ) data$country_UB_raw <<- data$country_UB
    if ( !is.null(data$sector_UB) ) data$sector_UB_raw <<- data$sector_UB
    
    # Upper and lower bounds
    # upper_mat <- wmat_bm + 0.02
    upper_mat <- wmat_bm + 0.05
    upper_mat[ upper_mat > 1 ] <- 1
    # lower_mat <- wmat_bm - 0.02
    lower_mat <- wmat_bm - 0.05
    lower_mat[ lower_mat < 0 ] <- 0
    if ( !is.null( data$upper_mat) ) {
      Names <- colnames(wmat_bm)
      idx_na <- which( is.na( data$upper_mat[ ,Names] ) )
      upper_mat[ idx_na ] <- NA
      lower_mat[ idx_na ] <- NA
    }
  
    # Sector constraints
    wmat_bm_sector <- groupWeights( wmat = wmat_bm, 
                                    group_mat = sector_mat )
    sector_UB <- wmat_bm_sector + 0.1
    # sector_UB <- wmat_bm_sector + 0.05
    sector_UB[ sector_UB > 1 ] <- 1
    sector_LB <- wmat_bm_sector - 0.1
    # sector_LB <- wmat_bm_sector - 0.05
    sector_LB[ sector_LB < 0 ] <- 0
    
    # Country constraints
    wmat_bm_country <- groupWeights( wmat = wmat_bm, 
                                     group_mat = country_mat )
    country_UB <- wmat_bm_country + 0.05
    country_UB[ country_UB > 1 ] <- 1
    country_LB <- wmat_bm_country - 0.05
    country_LB[ country_LB < 0 ] <- 0
    
    # Correction for countries which have weight < 2.5% in BM
    idx <- which( wmat_bm_country < 0.025 )
    if ( length(idx) > 0 ) {
      country_UB[ idx ] <- wmat_bm_country[ idx ] * 3
    }
    
    # Attach
    data$upper_mat <<- upper_mat
    data$lower_mat <<- lower_mat
    data$sector_LB <<- sector_LB
    data$sector_UB <<- sector_UB
    data$country_LB <<- country_LB
    data$country_UB <<- country_UB
    
    return( TRUE )
  }
  BacktestCustom$methods( MSCIMultifactorConstraints = BacktestCustom.MSCIMultifactorConstraints )
  
  
  
  
  # --------------------------------------------------------------------------
  BacktestCustom.selection.best_momentum <- function( rebdate  = NULL )
  {
    
  }
  
            