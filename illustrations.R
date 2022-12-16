  
  
  ############################################################################
  ### RANDOM PORTFOLIOS - ILLUSTRATIONS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     09.12.2022
  # First version:    09.12.2022
  # --------------------------------------------------------------------------
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  
  require(rgl)
  require(fBasics)
  require(colorRamps)
  
  
  # --------------------------------------------------------------------------
  # CUSTOM FUNCTIONS AND PARAMETERS
  # --------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  rdirichlet <- function( n, alpha ) 
  {
    # Copy from MCMCpack
    # see ?MCMCpack:rdirichlet
    l <- length(alpha)
    x <- matrix( rgamma(l * n, alpha), ncol = l, byrow = TRUE )
    sm <- x %*% rep(1, l)
    ans <- x / as.vector(sm)
    return( ans )
  }
  
  # --------------------------------------------------------------------------
  ddirichlet <- function ( x, alpha ) 
  {
    # Copy from MCMCpack
    # see ?MCMCpack:ddirichlet
    dirichlet1 <- function(x, alpha) {
      logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
      s <- sum((alpha - 1) * log(x))
      exp(sum(s) - logD)
    }
    if (!is.matrix(x)) 
      if (is.data.frame(x)) 
        x <- as.matrix(x)
    else x <- t(x)
    if (!is.matrix(alpha)) 
      alpha <- matrix(alpha, ncol = length(alpha), nrow = nrow(x), 
                      byrow = TRUE)
    if (any(dim(x) != dim(alpha))) 
      stop("Mismatch between dimensions of x and alpha in ddirichlet().\n")
    pd <- vector(length = nrow(x))
    for ( i in 1:nrow(x) ) {
      pd[i] <- dirichlet1( x[i, ], alpha[i, ] )
    }
    pd[ apply(x, 1, function(z) any(z < 0 | z > 1)) ] <- 0
    pd[ apply(x, 1, function(z) all.equal(sum(z), 1) != TRUE) ] <- 0
    return(pd)
  }
  
  # --------------------------------------------------------------------------
  transformShiftRotate <- function( samples )
  {
    n <- ncol(samples)
    N <- pracma::nullspace( matrix(1, ncol = n) )
    vertices <- (diag(1, n) %*% N)
    colnames(vertices) <- c("x", "y")
    shift <- vertices[1, ]
    vertices_shift <- sweep( vertices, 2, shift, "-" )
    samples_trans <- t( apply( samples, 1, function(x) { x %*% N } ) )
    samples_trans_shift <- sweep( samples_trans, 2, shift, "-" )
    theta <- -RP:::angle( v1 = vertices_shift[2, ], v2 = c(vertices_shift[2, 1], 0) )
    rotation_matrix <- matrix( c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2, byrow = TRUE )
    vertices_rot <- t( rotation_matrix %*% t( vertices_shift ) )
    samples_trans_rot <- t( apply( samples_trans_shift, 1, function(x) { rotation_matrix %*% x } ) )
    
    return( samples_trans_rot )
  }
  
  
  # Parameters
  n_sim <- 10^5
  n <- 3
  # colors_density <- divPalette( n = n_sim, name = "Spectral" )
  colors_density <- matlab.like( n = n_sim )
  
  
  
  par( mfrow = c(2, 3) )
  
  
  # --------------------------------------------------------------------------
  # NAIVE RP: FLAT DIRICHLET (UNIFORM SAMPLING OVER THE SIMPLEX)
  # --------------------------------------------------------------------------
  
  samples_flat <- rdirichlet( n = n_sim, alpha = rep(1, n) )
  samples_flat_2d <- transformShiftRotate( samples_flat )
  
  
  
  # plot3d( x = samples_flat[ ,1], y = samples_flat[ ,2], z = samples_flat[ ,3], col = colors_density[1] )
  plot( x = samples_flat_2d[ ,1], y = samples_flat_2d[ ,2], col = colors_density[1],
        pch = 19, cex = 0.3, ylab = "", xlab = "a)", xaxt = 'n', yaxt = 'n', bty = "n" )
  
  
  # --------------------------------------------------------------------------
  # NAIVE RP: SYMMETRIC DIRICHLET, MORE MASS IN THE CENTER
  # --------------------------------------------------------------------------
  
  lambda <- 4
  alpha <- rep(1, n) * 4
  dens_lambda4 <- apply( samples_flat, 1, function(x) { ddirichlet( x = x, alpha = alpha ) } )
  ordering <- order(dens_lambda4)
  colors <- colors_density
  colors[ordering] <- colors_density
  
  # plot3d( x = samples_flat[ ,1], y = samples_flat[ ,2], z = samples_flat[ ,3], col = colors )
  plot( x = samples_flat_2d[ ,1], y = samples_flat_2d[ ,2], col = colors, pch = 19, 
        cex = 0.3, ylab = "", xlab = "b)", xaxt = 'n', yaxt = 'n', bty = "n" )
  
  
  # # --------------------------------------------------------------------------
  # # NAIVE RP: SYMMETRIC DIRICHLET, MORE MASS AT THE FACES, EDGES AND VERTICES
  # # --------------------------------------------------------------------------
  # 
  # lambda <- 0.8
  # alpha <- rep(1, n) * lambda
  # dens_lambda08 <- apply( samples_flat, 1, function(x) { ddirichlet( x = x, alpha = alpha ) } )
  # ordering <- order(dens_lambda08)
  # colors <- colors_density
  # colors[ordering] <- colors_density
  # 
  # # plot3d( x = samples_flat[ ,1], y = samples_flat[ ,2], z = samples_flat[ ,3], col = colors )
  # plot( x = samples_flat_2d[ ,1], y = samples_flat_2d[ ,2], col = colors, pch = 19, 
  #       cex = 0.3, ylab = "", xlab = "", xaxt = 'n', yaxt = 'n', bty = "n" )
  
  
  # --------------------------------------------------------------------------
  # SIMPLE RP: ASYMMETRIC DIRICHLET (NON-UNIFORM SAMPLING OVER THE SIMPLEX)
  # --------------------------------------------------------------------------
  
  # alpha <- c(0.4, 0.35, 0.25) * n
  alpha <- c(0.5, 0.3, 0.2) * n * 2
  dens_asy <- apply( samples_flat, 1, function(x) { ddirichlet( x = x, alpha = alpha ) } )
  ordering <- order(dens_asy)
  colors <- colors_density
  colors[ordering] <- colors_density
  
  # plot3d( x = samples_flat[ ,1], y = samples_flat[ ,2], z = samples_flat[ ,3], col = colors )
  plot( x = samples_flat_2d[ ,1], y = samples_flat_2d[ ,2], col = colors, pch = 19, 
        cex = 0.3, ylab = "", xlab = "c)", xaxt = 'n', yaxt = 'n', bty = "n" )
  
  
  # --------------------------------------------------------------------------
  # SIMPLY REGULARIZED RP: SHADOW DIRICHLET (NON-UNIFORM SAMPLING OVER THE SIMPLEX)
  # --------------------------------------------------------------------------
  
  # Transformation Matrix
  M <- do.call( cbind, lapply( 1:n, FUN = function(k) { c( rep(0, k-1), rep( 1/(n-k+1), (n-k+1) ) ) } ) )
  samples_sd <- t( apply( samples_flat, 1, function(x) { M %*% x } ) )
  ordering <- order(dens_asy)
  samples_sd_2d <- transformShiftRotate( samples = samples_sd )
  
  # plot3d( x = samples_flat[ ,1], y = samples_flat[ ,2], z = samples_flat[ ,3], 
  #         col = "grey" )
  # points3d( x = samples_sd[ordering, 1], y = samples_sd[ordering, 2], z = samples_sd[ordering, 3], 
  #           col = colors_density, size = 5 )
  plot( x = samples_flat_2d[ ,1], y = samples_flat_2d[ ,2], col = "grey", pch = 19, 
        cex = 0.3, ylab = "", xlab = "d)", xaxt = 'n', yaxt = 'n', bty = "n" )
  points( x = samples_sd_2d[ordering, 1], y = samples_sd_2d[ordering, 2], col = colors_density, pch = 19, 
         cex = 0.3, ylab = "", xlab = "", xaxt = 'n', yaxt = 'n', bty = "n" )
  
  
  
  # --------------------------------------------------------------------------
  # GENERALLY REGULARIZED RP: LINEARLY CONSTRAINED SIMPLEX
  # --------------------------------------------------------------------------
  
  # Using acceptance-rejection
  
  Polytope <- list( A = rbind( c(1, 0, 0), 
                               c(0, 1, 0),
                               c(0, 0, 1),
                               c(1, 1, 0) ),
                    b = c(0.9, 0.5, 0.7, 0.77) )
  acc_rej <- apply( samples_flat, 1, function(x) { all( Polytope$A %*% x <= Polytope$b ) } )
  ordering <- order( dens_asy )
  colors <- colors_density
  colors[ordering] <- colors_density
  colors[acc_rej == FALSE] <- "grey"
  
  # plot3d( x = samples_flat[ ,1], y = samples_flat[ ,2], z = samples_flat[ ,3], col = colors )
  plot( x = samples_flat_2d[ ,1], y = samples_flat_2d[ ,2], col = colors, pch = 19, 
        cex = 0.3, ylab = "", xlab = "d)", xaxt = 'n', yaxt = 'n', bty = "n" )
  points( x = samples_flat_2d[acc_rej, 1], y = samples_flat_2d[acc_rej, 2], col = colors[acc_rej], pch = 19, 
         cex = 0.3, ylab = "", xlab = "d)", xaxt = 'n', yaxt = 'n', bty = "n" )
  
  
  # --------------------------------------------------------------------------
  # GENERALLY REGULARIZED RP: INTERSECTION OF SIMPLEX WITH ELLIPSOID
  # --------------------------------------------------------------------------
  
  covmat <- diag( c(0.5, 1, 0.7) )
  sigma <- apply( samples_flat, 1, function(x) { t(x) %*% covmat %*% x } )
  acc_rej <- abs( sigma - median(sigma) ) < 0.005
  ordering <- order(dens_asy)
  colors <- colors_density
  colors[ordering] <- colors_density
  colors[acc_rej == FALSE] <- "grey"
  
  # plot3d( x = samples_flat[ ,1], y = samples_flat[ ,2], z = samples_flat[ ,3], col = colors )
  plot( x = samples_flat_2d[ ,1], y = samples_flat_2d[ ,2], col = colors, pch = 19, 
        cex = 0.3, ylab = "", xlab = "f)", xaxt = 'n', yaxt = 'n', bty = "n" )
  points( x = samples_flat_2d[acc_rej, 1], y = samples_flat_2d[acc_rej, 2], col = colors[acc_rej], pch = 19, 
          cex = 0.3, ylab = "", xlab = "d)", xaxt = 'n', yaxt = 'n', bty = "n" )
  
  
  
  
  
  
  
  
  
  
  
  ###################################
  
  
  # Debris
  
  
  
  
  
  # n <- 3 # set the dimension
  # sigma <- cov( GPO::Data[ ,1:n] ) # set the covariance that define the ellipsoid x^T \Sigma x = c,  and c the volatility value
  # w_eqw <- rep(1/n, n)
  # cc <- t(w_eqw) %*% sigma %*% w_eqw * 1.01 # set the value of volatility
  #   
  # #constraints of the unit simplex (in general Ax<b could contain any linear inequality that does not make the space infeasible, e.g. regulatory constraints)
  # A <- -diag(n)
  # b <- rep(0, 1)
  # 
  # # \sum x_i = 1
  # Aeq <- matrix(rep(1, n), nrow = 1, ncol = n)
  # beq <- c(1)
  # 
  # N = pracma::nullspace(Aeq) # null space of the equality constraints
  # x0 <- rep(1, n)/n # center of the unit simplex
  # 
  # # apply the transformation
  # b = b - A %*% x0
  # A = A %*% N
  # # now Ax < b is a (n-1)-dimensional polytope
  # 
  # V = diag(n) #vertices of the unit simplex
  # V = V - kronecker(matrix(1, 1, n), matrix(x0, ncol = 1))
  # V = t(N) %*% V # vertices of the simplex in the lower dimension
  # 
  # x0 = rep(0, n-1) # the center of the simplex in the lower dimension
  # 
  # sigma_proj = t(N) %*% sigma %*% N 
  # center = -MASS::ginv(sigma_proj) %*% (t(N) %*% sigma) %*% (rep(1,n)/n)
  # 
  # R = cc + t(center) %*% sigma_proj %*% center - (rep(1,n)/n) %*% sigma %*% (rep(1,n)/n)
  # R = R[1]
  # # the covariance of the projected ellipsoid (x - center)^T \sigma_proj (x - center) = 1
  # sigma_proj = sigma_proj / R
  # 
  # 
  # 
  # 
  # n_sim <- 10^5
  # samples <- rdirichlet( n = n_sim, alpha = rep(1, n) )
  # samples_trans <- t( apply( samples, 1, function(x) { x %*% N } ) )
  # vertices <- (diag(1, n) %*% N)
  # colnames(vertices) <- c("x", "y")
  # shift <- vertices[1, ]
  #   
  #   
  #   
  # alpha <- c(0.5, 0.3, 0.2) * n
  # dens_asy <- apply( samples, 1, function(x) { ddirichlet( x = x, alpha = alpha ) } )
  # ordering <- order(dens_asy)
  # colors <- colors_density
  # colors[ordering] <- colors_density
  # 
  # plot( x = samples_trans[ ,1], y = samples_trans[ ,2], col = colors, pch = 19, cex = 0.3 )
  # points( x = vertices[ ,1], vertices[ ,2], pch = 19, col = 1, cex = 2 )
  # 
  # 
  # 
  # # Shift and rotate
  # 
  # vertices_shift <- sweep( vertices, 2, shift, "-" )
  # samples_trans_shift <- sweep( samples_trans, 2, shift, "-" )
  # 
  # theta <- -RP:::angle( v1 = vertices_shift[2, ], v2 = c(vertices_shift[2, 1], 0) )
  # rotation_matrix <- matrix( c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2, byrow = TRUE )
  # rotation_matrix
  # 
  # vertices_rot <- t( rotation_matrix %*% t( vertices_shift ) )
  # samples_trans_rot <- t( apply( samples_trans_shift, 1, function(x) { rotation_matrix %*% x } ) )
  # 
  # 
  # plot( x = samples_trans_shift[ ,1], y = samples_trans_shift[ ,2], col = colors, pch = 19, cex = 0.3 )
  # points( x = samples_trans_rot[ ,1], samples_trans_rot[ ,2], pch = 19, col = colors, cex = 1 )
  # points( x = vertices_shift[ ,1], vertices_shift[ ,2], pch = 19, col = 1, cex = 2 )
  # points( x = vertices_rot[ ,1], vertices_rot[ ,2], pch = 19, col = 2, cex = 2 )
  
  
  
 
  
  