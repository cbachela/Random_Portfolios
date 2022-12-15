  
  
  ############################################################################
  ### GEOMSCALE - RANDOM PORTFOLIOS - CUSTOM FUNCTIONS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard   21.06.2022
  # First version:    21.06.2022
  # --------------------------------------------------------------------------
  
 
  # letterScore
  # simSortedByMom
  # get_sequence_of_scores
  # mn2Lambda
  
  
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
                           w_array = NULL )
  {
    
    sim_tmp <- sim
    sim_new <- sim_tmp * NA
    id <- dimnames(w_array)[[2]]
    rebdates_d <- rownames(sim_tmp)
    # rebdates <- BT0$spec$rebdates
    rebdates <- rownames(w_array)
    rebdates[1] <- rebdates_d[1]
    
    for ( i in 2:(length(rebdates) - 1) ) {
      
      # Out-of-sample
      X_tmp <- window( BT0$data$X_sim[ ,id], rebdates[i-1], rebdates[i] )
      X_tmp[is.na(X_tmp)] <- 0
      today <- rebdates[i]
      lW <- lapply( 1:dim(w_array)[[3]], FUN = function(j) { w_array[today, ,j] } )
      sim <- do.call( cbind, lapply( lW, function(w) { w[is.na(w)] <- 0; X_tmp %*% w } ) )
      
      Names <- intersect( colnames(w_array[ ,,1]), colnames(Score) )
      scores <- as.numeric( Score[today, Names] )
      FUN <- function(i)
      {
        wghts <- as.numeric( w_array[today, Names, i] )
        wghts[is.na(wghts)] <- 0
        ans <- as.numeric( wghts %*% scores )
        return( ans )
      }
      score_vec <- unlist( lapply( 1:ncol(sim), FUN = FUN ) )
      ordering <- order(score_vec)
      idx3 <- which(rebdates_d == rebdates[i+1])
      idx2 <- which(rebdates_d == rebdates[i])
      idx1 <- which(rebdates_d == rebdates[i-1])
      if ( i > 2 ) {
        idx1 <- idx1 + 1
      }
      X_tmp_1_2 <- sim_tmp[idx1:idx2, ]
      X_tmp_2_3 <- sim_tmp[(idx2+1):idx3, ]
      sim_new[(idx2+1):idx3, ] <- sim_tmp[(idx2+1):idx3, ordering]
      if ( i == 2 ) {
        sim_new[idx1:idx2, ] <- sim_tmp[idx1:idx2, ordering]
      }
      # # In-sample
      # X_tmp <- sim_tmp[idx1:idx2, ]
      # mu <- exp( apply( log(1 + X_tmp), 2, mean ) ) - 1
      # eps <- rnorm(length(mu), 1, k)
      # mu <- mu * eps
      # ordering <- order(mu)
      # sim_new[idx1:idx2, ] <- X_tmp[ ,ordering]
    }
    sim_new <- na.omit(sim_new)
    
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
  
  
  
  