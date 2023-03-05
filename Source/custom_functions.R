  
  
  ############################################################################
  ### GEOMSCALE - RANDOM PORTFOLIOS - CUSTOM FUNCTIONS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard   21.06.2022
  # First version:    21.06.2022
  # --------------------------------------------------------------------------
  
 
  # letterScore
  # simSortedBy
  # get_sequence_of_scores
  # mn2Lambda
  # lincon2Polytope
  
  
  # --------------------------------------------------------------------------
  letterScore <- function( ticker, preference = "A" ) 
  {
    ans <- 0
    if ( !is.na(ticker) ) {
      ans <- stringr::str_count( ticker, preference )
    }
    return( ans )
  }
  
  
  # --------------------------------------------------------------------------
  simSortedBy <- function( sim = NULL, 
                           BT0 = NULL, 
                           Score = NULL,
                           w_array = NULL,
                           insample = FALSE )
  {
    Score[is.na(Score)] <- 0
    sim_tmp <- sim
    sim_new <- sim_tmp * NA
    id <- dimnames(w_array)[[2]]
    # rebdates <- BT0$spec$rebdates
    rebdates <- intersect( rownames(w_array), rownames(Score) )
    rebdates_d <- rownames(sim_tmp)
    # rebdates[1] <- rebdates_d[1]
    sim_scores <- as.timeSeries( matrix(NA, nrow = length(rebdates), ncol = ncol(sim), 
                                        dimnames = list(rebdates, colnames(sim)) ), rebdates )
    ordering_mat <- sim_scores
    
    FUN <- function(i)
    {
      wghts <- as.numeric( w_array[yesterday, Names, i] )
      wghts[is.na(wghts)] <- 0
      ans <- as.numeric( wghts %*% scores )
      return( ans )
    }
    
    for ( i in 2:(length(rebdates) - 1) ) {
      
      today <- rebdates[i]
      yesterday <- rebdates[i-1]
      Names <- intersect( colnames(w_array[ ,,1]), colnames(Score) )
      idx1 <- which(rebdates_d == as.character(as.Date(rebdates[i-1]) + 1))
      idx2 <- which(rebdates_d == rebdates[i])

      if ( isTRUE(insample) ) {
        
        # In-sample
        scores <- setNames( as.numeric( Score[today, Names] ), Names )
        score_vec <- unlist( lapply( 1:dim(w_array)[[3]], FUN = FUN ) )
        ordering <- order(score_vec)
        ordering_mat[today, ] <- ordering
      
        X_tmp <- sim_tmp[idx1:idx2, ]
        # mu <- exp( apply( log(1 + X_tmp), 2, mean ) ) - 1
        # ordering <- order(mu)
        # sds <- apply( X_tmp, 2, sd )
        # ordering <- order(sds)
        sim_new[idx1:idx2, ] <- X_tmp[ ,ordering]
        sim_scores[today, ] <- score_vec[ordering]
        
      } else {
        
        # Out-of-sample
        scores <- setNames( as.numeric( Score[yesterday, Names] ), Names )
        score_vec <- unlist( lapply( 1:dim(w_array)[[3]], FUN = FUN ) )
        ordering <- order(score_vec)
        ordering_mat[yesterday, ] <- ordering
        
        X_tmp <- sim_tmp[idx1:idx2, ]
        sim_new[idx1:idx2, ] <- X_tmp[ ,ordering]
        sim_scores[yesterday, ] <- score_vec[ordering]
        
        # X_tmp <- window( BT0$data$X_sim[ ,id], rebdates[i-1], rebdates[i] )
        # X_tmp[is.na(X_tmp)] <- 0
        # lW <- lapply( 1:dim(w_array)[[3]], FUN = function(j) { w_array[today, ,j] } )
        # sim <- do.call( cbind, lapply( lW, function(w) { w[is.na(w)] <- 0; X_tmp %*% w } ) )
        # 
        # 
        # idx3 <- which(rebdates_d == rebdates[i+1])
        # idx2 <- which(rebdates_d == rebdates[i])
        # idx1 <- which(rebdates_d == as.character(as.Date(rebdates[i-1]) + 1))
        # if ( i > 2 ) {
        #   idx1 <- idx1 + 1
        # }
        # X_tmp_1_2 <- sim_tmp[idx1:idx2, ]
        # X_tmp_2_3 <- sim_tmp[(idx2+1):idx3, ]
        # sim_new[(idx2+1):idx3, ] <- sim_tmp[(idx2+1):idx3, ordering]
        # if ( i == 2 ) {
        #   sim_new[idx1:idx2, ] <- sim_tmp[idx1:idx2, ordering]
        # }
      }
    }
    sim_new <- na.omit(sim_new)
    attr(sim_new, "scores") <- sim_scores
    attr(sim_new, "ordering") <- ordering_mat
    
    return( sim_new )
  }
  
  
  
  # --------------------------------------------------------------------------
  get_sequence_of_scores <- function( scores, m )
  {
    # Sort by increasing scores
    ordering <- order(scores)
    
    # Estimate the variance of an equally-weighted portfolio within 
    # each vola-quantile
    portfolio_scores <- numeric(m)
    lIdx <- list()
    for ( i in 1:m ) {
      k <- ceiling( length(ordering) / m )   # ~~~~~~~~~~~
      if ( i < m ) {
        idx <- ordering[(i*k-k+1):(i*k)]
      } else {
        idx <- ordering[(i*k-k+1):length(ordering)]
      }
      lIdx[[i]] <- idx
      wghts <- rep(1/length(idx), length(idx))
      portfolio_scores[i] <- t(wghts) %*% scores[idx]
    }
    attr(portfolio_scores, "lIdx") <- lIdx
    
    return( portfolio_scores )
  }
  
  
  
  # --------------------------------------------------------------------------
  # Function finds lambda such that the variance of the Bayesian bootstrap
  # based on Dir(\alpha \lambda) corresponds to the variance of the m out of n bootstrap.
  # --------------------------------------------------------------------------
  mn2Lambda <- function(n, m) { (1-m)/(1-n) + (m-1)/(n*(1-n)) }
  
  
  
  # --------------------------------------------------------------------------
  # Part of the AAA package
  lincon2Polytope <- function(object)
  {
    
    selection <- getConstraints(object, "selection")
    boxcon <- getConstraints(object, "bounds")
    lincon <- getConstraints(object, "linear")
    if ( length(selection) < length(boxcon$upper) ) {
      selection <- names(boxcon$upper)
    }
    
    n <- length(selection)
    stopifnot( length(boxcon$lower) == n )
    stopifnot( length(boxcon$upper) == n )
    
    # Bounds
    A <- rbind(diag(n) * (-1),
               diag(n) )
    colnames(A) <- selection
    b <- c(boxcon$lower * (-1), boxcon$upper)
    sense <- rep("<=", n * 2)
    
    P <- list(A = A, 
              b = b, 
              sense = sense)
    A = b = sense = NULL
    
    # Linear constraints
    if ( !is.null(lincon$Amat) ) {
      
      A <- lincon$Amat
      rhs <- lincon$rhs
      sense <- lincon$sense
      
      idx_geq <- which(sense == ">=")
      if ( length(idx_geq) > 0 ) {
        A[idx_geq, ] <- A[idx_geq, ] * (-1)
        rhs[idx_geq] <- rhs[idx_geq] * (-1)
        sense[idx_geq] <- "<="
      }
      
      P$A <- rbind(P$A, A)
      P$b <- c(P$b, rhs)
      P$sense <- c(P$sense, sense)
    }
    
    colnames(P$A) <- selection
    
    return( P )
  }
  
  