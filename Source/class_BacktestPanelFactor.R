
  
  ############################################################################
  ### CLASS BacktestPanelFactorGRW
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard   06.01.2023
  # First version:    06.01.2023
  # --------------------------------------------------------------------------

  
  # BacktestPanelFactor
  # BacktestPanelFactor.selection.factor
  # # BacktestPanelFactor.setSolver
  # BacktestPanelFactor.constraints.equalbmex
  # BacktestPanelFactor.constraints.equalbmexcap
  # BacktestPanelFactor.constraints.controlbmex
  # BacktestPanelFactor.constraints.controlbmex25
  # BacktestPanelFactor.constraints.controlbmex10
  # BacktestPanelFactor.constraints.controlbmex10ub
  # BacktestPanelFactor.constraints.controlbmex10lb
  # BacktestPanelFactor.constraints.controlbmex00
  # BacktestPanelFactor.constraints.controlbmex00ub
  # BacktestPanelFactor.constraints.controlbmex00lb
  # BacktestPanelFactor.factorScores
  # BacktestPanelFactor.centroidBuckets
  # BacktestPanelFactor.standardizeFactorScores
  # BacktestPanelFactor.quantilePortfolio
  # BacktestPanelFactor.decilePortfolio
  # BacktestPanelFactor.quintilePortfolio
  # BacktestPanelFactor.minTEPortfolio
  # BacktestPanelFactor.maxExposurePortfolio
  # BacktestPanelFactor.constraints.longshort 
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Backtest Class Definition
  # --------------------------------------------------------------------------
  
  BacktestPanelFactor <- setRefClass( Class = "BacktestPanelFactor", 
                                      contains = "BacktestCustom",
                                      methods = list() )
  
  
  
  # --------------------------------------------------------------------------
  # Main Methods
  # --------------------------------------------------------------------------

  # # --------------------------------------------------------------------------
  # BacktestPanelFactor.gps <- function(rebdate = NULL)
  # {
  #   # Copy the initial GPS object from the specs slot
  #   # and assign it to the rebalancing environment.
  #   .self$setGPS( rebdate = rebdate )
  #   # Modify the selection in the rebalancing environment
  #   .self$setSelection( rebdate = rebdate )
  #   # Modify the GPS object in the rebalancing environment
  #   .self$setData( rebdate = rebdate )
  #   .self$setCovariance( rebdate = rebdate )
  #   .self$setSolver( rebdate = rebdate ) # Computes weighted factor score vector
  #   .self$setConstraints( rebdate = rebdate )
  #   return( TRUE )
  # }
  # BacktestPanelFactor$methods( gps = BacktestPanelFactor.gps )
  
  
  
  
  # --------------------------------------------------------------------------
  BacktestPanelFactor.selection.factor <- function(rebdate = NULL)
  {
    factor_names_selection <- spec$factor_names_selection
    stopifnot( is.character(factor_names_selection) )
    p_sel <- spec$p_sel
    if ( is.null(p_sel) ) p_sel <- 1
    if ( length(p_sel) < factor_names_selection ) {
      p_sel <- rep( p_sel, length(factor_names_selection) )
    }
    if ( !identical( names(p_sel), factor_names_selection) ) {
      p_sel <- setNames( p_sel, factor_names_selection )
    }
        
    selection <- colnames(data$X_est)
    
    for ( factor_name in factor_names_selection ) {
      
      factor_mat <- data$lFactor[[ factor_name ]]
      
      if ( !is.null(factor_mat) ) {
        
        selection <- intersect( selection, colnames(factor_mat) )
        idx <- which( rownames(factor_mat) == rebdate )
        fac_vec <- setNames( as.numeric(factor_mat[idx, selection]), selection )
        names_fac <- names(fac_vec[ which(!is.na(fac_vec)) ])
        selection <- intersect( selection, names_fac )
        
        # Filter out a given quantile of the factor scores
        if ( p_sel < 1 ) {
          
          fac_vec <- fac_vec[selection]
          names_fac_vec <- names( which( fac_vec < quantile( fac_vec, p_sel[factor_name] ) ) )
          selection <- intersect( selection, names_fac_vec )
        }
      }
    }
    return(selection)
  }
  BacktestPanelFactor$methods( selection.factor = BacktestPanelFactor.selection.factor )
  
  
  
  # # --------------------------------------------------------------------------
  # BacktestPanelFactor.setSolver <- function( rebdate = NULL )
  # {
  #   # Add weighted factor scores to GPS in rebalancing environment
  #   RBE <- rebalenv( rebdate = rebdate )
  #   selection <- RBE$selection
  #   Solver <- getSolver(RBE$GPS)
  #   # Add factor timing to GPS in rebalancing environment
  #   factor_list <- .self$factorScores( rebdate = rebdate )
  #   weighted_scores <- factor_list$weighted_scores
  #   Solver$obj_lin <- setNames(as.numeric(weighted_scores * (-1)), names(weighted_scores))
  #   setSolver(RBE$GPS) <- Solver
  #   return( TRUE )
  # }
  # BacktestPanelFactor$methods( setSolver = BacktestPanelFactor.setSolver )
  
  
  
  # --------------------------------------------------------------------------
  BacktestPanelFactor.constraints.equalbmex <- function(rebdate = NULL) 
  {
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    stopifnot( is.character(spec$factor_names_constraints) )
    factor_names <- spec$factor_names_constraints
    factor_list <- .self$factorScores( rebdate = rebdate )
    factor_mat <- do.call(rbind, lapply(factor_list, as.numeric))
    factor_scores_bm <- data$factor_scores_bm[rebdate, factor_names]
    
    Amat <- matrix(factor_mat[factor_names, ], nrow = length(factor_names))
    rownames(Amat) <- factor_names
    colnames(Amat) <- selection
    rhs <- as.numeric(factor_scores_bm)
    sense <- rep("=", nrow(Amat))
    ans <- linearConstraint( rhs = rhs,
                             sense = sense,
                             Amat = Amat )
    return(ans)
  }
  BacktestPanelFactor$methods( constraints.equalbmex = BacktestPanelFactor.constraints.equalbmex )
  
  
  
  # --------------------------------------------------------------------------
  BacktestPanelFactor.constraints.equalbmexcap <- function(rebdate = NULL) 
  {
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    stopifnot( is.character(spec$factor_names_constraints) )
    factor_names <- spec$factor_names_constraints
    factor_list <- .self$factorScores( rebdate = rebdate )
    factor_mat <- do.call(rbind, lapply(factor_list, as.numeric))
    factor_scores_bm <- data$factor_scores_bm[rebdate, factor_names]
    
    Amat <- matrix(factor_mat[factor_names, ], nrow = length(factor_names))
    rownames(Amat) <- factor_names
    colnames(Amat) <- selection
    factor_scores_bm[factor_scores_bm < -1.25] <- -1.25
    factor_scores_bm[factor_scores_bm > 1.25] <- 1.25
    rhs <- as.numeric(factor_scores_bm)
    sense <- rep("=", nrow(Amat))
    ans <- linearConstraint(rhs = rhs,
                            sense = sense,
                            Amat = Amat)
    return(ans)
  }
  BacktestPanelFactor$methods( constraints.equalbmexcap = BacktestPanelFactor.constraints.equalbmexcap )
  
  
  
  # --------------------------------------------------------------------------
  BacktestPanelFactor.constraints.controlbmex <- function(rebdate = NULL) 
  {
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    stopifnot( is.character(spec$factor_names_constraints) )
    factor_names <- spec$factor_names_constraints
    factor_list <- .self$factorScores( rebdate = rebdate )
    factor_mat <- do.call(rbind, lapply(factor_list, as.numeric))
    factor_scores_bm <- data$factor_scores_bm[rebdate, factor_names]
    
    Amat <- matrix(factor_mat[factor_names, ], nrow = length(factor_names))
    Amat <- rbind(Amat, Amat)
    rownames(Amat) <- c(paste0(factor_names, "_ub"), paste0(factor_names, "_lb"))
    colnames(Amat) <- selection
    factor_scores_bm[factor_scores_bm < -1.25] <- -1.25
    factor_scores_bm[factor_scores_bm > 1.25] <- 1.25
    rhs <- c(as.numeric(factor_scores_bm) + 0.1,
             as.numeric(factor_scores_bm) - 0.1)
    sense <- c(rep("<=", length(factor_names)), rep(">=", length(factor_names)))
    ans <- linearConstraint(rhs = rhs,
                            sense = sense,
                            Amat = Amat)
    return(ans)
  }
  BacktestPanelFactor$methods( constraints.controlbmex = BacktestPanelFactor.constraints.controlbmex )
  
 
  
  # --------------------------------------------------------------------------
  BacktestPanelFactor.constraints.controlbmex25 <- function(rebdate = NULL) 
  {
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    stopifnot( is.character(spec$factor_names_constraints) )
    factor_names <- spec$factor_names_constraints
    factor_list <- .self$factorScores( rebdate = rebdate )
    factor_mat <- do.call(rbind, lapply(factor_list, as.numeric))
    factor_scores_bm <- data$factor_scores_bm[rebdate, factor_names]
    
    Amat <- matrix(factor_mat[factor_names, ], nrow = length(factor_names))
    Amat <- rbind(Amat, Amat)
    rownames(Amat) <- c(paste0(factor_names, "_ub"), paste0(factor_names, "_lb"))
    colnames(Amat) <- selection
    # factor_scores_bm[factor_scores_bm < -1.25] <- -1.25
    # factor_scores_bm[factor_scores_bm > 1.25] <- 1.25
    rhs <- c(as.numeric(factor_scores_bm) + 0.25,
             as.numeric(factor_scores_bm) - 0.25)
    sense <- c(rep("<=", length(factor_names)), rep(">=", length(factor_names)))
    ans <- linearConstraint(rhs = rhs,
                            sense = sense,
                            Amat = Amat)
    return(ans)
  }
  BacktestPanelFactor$methods( constraints.controlbmex25 = BacktestPanelFactor.constraints.controlbmex25 )
  
 
  
  # --------------------------------------------------------------------------
  BacktestPanelFactor.constraints.controlbmex10 <- function(rebdate = NULL) 
  {
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    stopifnot( is.character(spec$factor_names_constraints) )
    factor_names <- spec$factor_names_constraints
    factor_list <- .self$factorScores( rebdate = rebdate )
    factor_mat <- do.call(rbind, lapply(factor_list, as.numeric))
    factor_scores_bm <- data$factor_scores_bm[rebdate, factor_names]
    
    Amat <- matrix(factor_mat[factor_names, ], nrow = length(factor_names))
    Amat <- rbind(Amat, Amat)
    rownames(Amat) <- c(paste0(factor_names, "_ub"), paste0(factor_names, "_lb"))
    colnames(Amat) <- selection
    # factor_scores_bm[factor_scores_bm < -1.25] <- -1.25
    # factor_scores_bm[factor_scores_bm > 1.25] <- 1.25
    rhs <- c(as.numeric(factor_scores_bm) + 0.1,
             as.numeric(factor_scores_bm) - 0.1)
    sense <- c(rep("<=", length(factor_names)), rep(">=", length(factor_names)))
    ans <- linearConstraint(rhs = rhs,
                            sense = sense,
                            Amat = Amat)
    return(ans)
  }
  BacktestPanelFactor$methods( constraints.controlbmex10 = BacktestPanelFactor.constraints.controlbmex10 )
  
  
  
  # --------------------------------------------------------------------------
  BacktestPanelFactor.constraints.controlbmex10ub <- function(rebdate = NULL) 
  {
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    stopifnot( is.character(spec$factor_names_constraints) )
    factor_names <- spec$factor_names_constraints
    factor_list <- .self$factorScores( rebdate = rebdate )
    factor_mat <- do.call(rbind, lapply(factor_list, as.numeric))
    factor_scores_bm <- data$factor_scores_bm[rebdate, factor_names]
    
    Amat <- matrix(factor_mat[factor_names, ], nrow = length(factor_names))
    rownames(Amat) <- paste0(factor_names, "_ub")
    colnames(Amat) <- selection
    rhs <- as.numeric(factor_scores_bm) + 0.1
    sense <- rep("<=", length(factor_names))
    ans <- linearConstraint(rhs = rhs,
                            sense = sense,
                            Amat = Amat)
    return(ans)
  }
  BacktestPanelFactor$methods( constraints.controlbmex10ub = BacktestPanelFactor.constraints.controlbmex10ub )
  
  
  
  # --------------------------------------------------------------------------
  BacktestPanelFactor.constraints.controlbmex10lb <- function(rebdate = NULL) 
  {
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    stopifnot( is.character(spec$factor_names_constraints) )
    factor_names <- spec$factor_names_constraints
    factor_list <- .self$factorScores( rebdate = rebdate )
    factor_mat <- do.call(rbind, lapply(factor_list, as.numeric))
    factor_scores_bm <- data$factor_scores_bm[rebdate, factor_names]
    
    Amat <- matrix(factor_mat[factor_names, ], nrow = length(factor_names))
    rownames(Amat) <- paste0(factor_names, "_lb")
    colnames(Amat) <- selection
    rhs <- as.numeric(factor_scores_bm) - 0.1
    sense <- rep(">=", length(factor_names))
    ans <- linearConstraint(rhs = rhs,
                            sense = sense,
                            Amat = Amat)
    return(ans)
  }
  BacktestPanelFactor$methods( constraints.controlbmex10lb = BacktestPanelFactor.constraints.controlbmex10lb )
  
  
 
  
  # --------------------------------------------------------------------------
  BacktestPanelFactor.constraints.controlbmex00 <- function(rebdate = NULL) 
  {
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    stopifnot( is.character(spec$factor_names_constraints) )
    factor_names <- spec$factor_names_constraints
    factor_list <- .self$factorScores( rebdate = rebdate )
    factor_mat <- do.call(rbind, lapply(factor_list, as.numeric))
    factor_scores_bm <- data$factor_scores_bm[rebdate, factor_names]
    
    Amat <- matrix(factor_mat[factor_names, ], nrow = length(factor_names))
    rownames(Amat) <- paste0(factor_names, "_eb")
    colnames(Amat) <- selection
    # factor_scores_bm[factor_scores_bm < -1.25] <- -1.25
    # factor_scores_bm[factor_scores_bm > 1.25] <- 1.25
    rhs <- as.numeric(factor_scores_bm)
    sense <- rep("=", length(factor_names))
    ans <- linearConstraint(rhs = rhs,
                            sense = sense,
                            Amat = Amat)
    return(ans)
  }
  BacktestPanelFactor$methods( constraints.controlbmex00 = BacktestPanelFactor.constraints.controlbmex00 )
  
  
  # --------------------------------------------------------------------------
  BacktestPanelFactor.constraints.controlbmex00ub <- function(rebdate = NULL) 
  {
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    stopifnot( is.character(spec$factor_names_constraints) )
    factor_names <- spec$factor_names_constraints
    factor_list <- .self$factorScores( rebdate = rebdate )
    factor_mat <- do.call(rbind, lapply(factor_list, as.numeric))
    factor_scores_bm <- data$factor_scores_bm[rebdate, factor_names]
    
    Amat <- matrix(factor_mat[factor_names, ], nrow = length(factor_names))
    rownames(Amat) <- paste0(factor_names, "_ub")
    colnames(Amat) <- selection
    # factor_scores_bm[factor_scores_bm < -1.25] <- -1.25
    # factor_scores_bm[factor_scores_bm > 1.25] <- 1.25
    rhs <- as.numeric(factor_scores_bm)
    sense <- rep("<=", length(factor_names))
    ans <- linearConstraint(rhs = rhs,
                            sense = sense,
                            Amat = Amat)
    return(ans)
  }
  BacktestPanelFactor$methods( constraints.controlbmex00ub = BacktestPanelFactor.constraints.controlbmex00ub )
  
 
  # --------------------------------------------------------------------------
  BacktestPanelFactor.constraints.controlbmex00lb <- function(rebdate = NULL) 
  {
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    stopifnot( is.character(spec$factor_names_constraints) )
    factor_names <- spec$factor_names_constraints
    factor_list <- .self$factorScores( rebdate = rebdate )
    factor_mat <- do.call(rbind, lapply(factor_list, as.numeric))
    factor_scores_bm <- data$factor_scores_bm[rebdate, factor_names]
    
    Amat <- matrix(factor_mat[factor_names, ], nrow = length(factor_names))
    rownames(Amat) <- paste0(factor_names, "_lb")
    colnames(Amat) <- selection
    # factor_scores_bm[factor_scores_bm < -1.25] <- -1.25
    # factor_scores_bm[factor_scores_bm > 1.25] <- 1.25
    rhs <- as.numeric(factor_scores_bm)
    sense <- rep(">=", length(factor_names))
    ans <- linearConstraint(rhs = rhs,
                            sense = sense,
                            Amat = Amat)
    return(ans)
  }
  BacktestPanelFactor$methods( constraints.controlbmex00lb = BacktestPanelFactor.constraints.controlbmex00lb )
  

  
  # --------------------------------------------------------------------------
  BacktestPanelFactor.factorScores <- function(rebdate = NULL)
  {
    # Extract factor scores and standardize them
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    FUN <- function( name, rebdate, selection ) { data$lFactor[[ name ]][ rebdate, selection ] }
    factor_list <- lapply(spec$factor_names, FUN = FUN, rebdate = rebdate, selection = selection )
    # don't standardize in the optimization
    # factor_list <- lapply(factor_list, .self$standardizeFactorScores)
    names(factor_list) <- spec$factor_names
    factor_mat <- do.call(cbind, lapply(factor_list, as.numeric))
    rownames(factor_mat) <- selection
    # Append factor scores to output
    .self$appendOutput(factor_list)
    return(factor_list)
  }
  BacktestPanelFactor$methods( factorScores = BacktestPanelFactor.factorScores )
  

  
  # --------------------------------------------------------------------------
  BacktestPanelFactor.centroidBuckets <- function(buckets, nsim = 1000)
  {
    if ( !is.list(buckets) ) {
      stop("buckets must be a list")
    }
    nassets <- length(unlist(buckets))
    nbuckets <- length(buckets)
    c.hat <- matrix(0, nsim, nbuckets)
    for(i in 1:nsim){
      c.hat[i, ] <- sort(rnorm(nbuckets), decreasing = FALSE)
    }
    mu <- apply(c.hat, 2, mean)
    out <- vector("numeric", nassets)
    for ( j in 1:nbuckets ) {
      out[buckets[[j]]] <- mu[j]
    }
    return(out)
  }
  BacktestPanelFactor$methods( centroidBuckets = BacktestPanelFactor.centroidBuckets )
  
  
  
  # --------------------------------------------------------------------------
  BacktestPanelFactor.standardizeFactorScores <- function( scores, seed = 2019 )
  {
    # Z-score
    scores <- (scores - mean(scores)) / sd(scores)
    # Winsorization
    scores[scores<=-3] <- -3
    scores[scores>=3] <- 3
    # Centroid
    tmp <- lapply( X = sort(unique(scores)), 
                   FUN = function ( X ) { which(scores == X) } )
    if ( !is.null(seed) ) {
      set.seed(seed)
    }
    scores <- setNames( BacktestFactorMinVar.centroidBuckets(buckets = tmp, nsim = 1000), 
                        names(scores) )
    
    return( scores )
  }
  BacktestPanelFactor$methods( standardizeFactorScores = BacktestPanelFactor.standardizeFactorScores )
  
  
  
  # --------------------------------------------------------------------------
  BacktestPanelFactor.quantilePortfolio <- function( rebdate = NULL )
  {
    
    RBE <- rebalenv( rebdate = rebdate )
    selection <- RBE$selection
    wmat_bm <- data$wmat_bm
    q_vec <- spec$q_vec
    weighting_method <- spec$weighting_method
    factor_name <- spec$factor_name
    
    # Extract factor scores
    factor_vec <- data$lFactor[[ factor_name ]][ rebdate, selection ]
    
    # Sort by increasing score
    scores <- sort(factor_vec)
    ordering <- order(scores)
    
    # Identify stocks within quantile
    m <- length(q_vec)
    lIdx <- list()
    lWeights <- list()
    for ( i in 1:m ) {
      k <- ceiling( length(ordering) / m )
      if ( i < m ) {
        idx <- ordering[(i*k-k+1):(i*k)]
      } else {
        idx <- ordering[(i*k-k+1):length(ordering)]
      }
      lIdx[[i]] <- idx
      id <- selection[idx]
      if ( weighting_method == "eqw" ) {
        wghts <- setNames( rep(1/length(id), length(id)), id ) 
      } else {
        wghts_bm <- setNames( as.numeric(wmat_bm[rebdate, id]), id )
        wghts <- wghts_bm / sum(wghts_bm)
      }
      lWeights[[i]] <- wghts
    }
    
    # Make time series and append
    FUN <- function(w) { timeSeries( matrix( w, nrow = 1, dimnames = list(NULL, names(w)) ), rebdate ) }
    lWeights <- lapply( lWeights, FUN = FUN )
    names(lWeights) <- paste0("q", q_vec)
    .self$appendOutput( lWeights )
    
    return( TRUE )
  }
  BacktestPanelFactor$methods( quantilePortfolio = BacktestPanelFactor.quantilePortfolio )
  
  
  # --------------------------------------------------------------------------
  BacktestPanelFactor.decilePortfolio <- function( rebdate = NULL )
  {
    spec$q_vec <<- seq(from = 0.1, to = 1, by = 0.1)
    .self$quantilePortfolio( rebdate = rebdate )
    return( TRUE )
  }
  BacktestPanelFactor$methods( decilePortfolio = BacktestPanelFactor.decilePortfolio )
  
  
  # --------------------------------------------------------------------------
  BacktestPanelFactor.quintilePortfolio <- function( rebdate = NULL )
  {
    spec$q_vec <<- seq(from = 0.2, to = 1, by = 0.2)
    .self$quantilePortfolio( rebdate = rebdate )
    return( TRUE )
  }
  BacktestPanelFactor$methods( quintilePortfolio = BacktestPanelFactor.quintilePortfolio )
  
  
  
  
  
  # --------------------------------------------------------------------------
  BacktestPanelFactor.minTEPortfolio <- function( rebdate = NULL )
  {
    RBE <- rebalenv( rebdate = rebdate )
    .self$benchmarkWeights( rebdate = rebdate )
    RBE$GPS@solver$portfolio <- "mintrackingerror"
    RBE$GPS@solver$obj_quad <- covariance( RBE$GPS )
    RBE$GPO <- GPO::gpo( RBE$GPS )
    w <- getWeights(RBE$GPO)
    wghts <- timeSeries( matrix( w, nrow = 1, dimnames = list(NULL, names(w)) ), rebdate )
    .self$appendOutput( list( weights = wghts ) )
    return( TRUE )
  }
  BacktestPanelFactor$methods( minTEPortfolio = BacktestPanelFactor.minTEPortfolio )
  
  
  
  
  # --------------------------------------------------------------------------
  BacktestPanelFactor.maxExposurePortfolio <- function( rebdate = NULL )
  {
    RBE <- rebalenv( rebdate = rebdate )
    GPS <- RBE$GPS
    selection <- RBE$selection
    factor_name <- spec$factor_name
    direction <- spec$direction
    
    # Extract factor scores
    score <- data$lFactor[[ factor_name ]][ rebdate, selection ]
    
    GPS@solver$obj_lin <- score * direction
    GPS@solver$utility <- list(riskaversion = 0)
    GPS@solver$portfolio <- "meanvariance"
    GPO <- GPO::gpo(GPS)
    wghts <- GPO::getWeights(GPO)
    wghts <- timeSeries( matrix( wghts, nrow = 1, dimnames = list(NULL, names(wghts)) ),
                         rebdate )
    .self$appendOutput( value = list( weights = wghts ) )
    RBE$GPS <- GPS
    RBE$GPO <- GPO
    return( TRUE )
  }
  BacktestPanelFactor$methods( maxExposurePortfolio = BacktestPanelFactor.maxExposurePortfolio )
  
  
  
  # --------------------------------------------------------------------------
  BacktestPanelFactor.constraints.longshort  <- function(rebdate = NULL) 
  {
    ans <- boxConstraint(name = "LongShort")
    return(ans)
  }
  BacktestPanelFactor$methods( constraints.longshort = BacktestPanelFactor.constraints.longshort )
  
  
            