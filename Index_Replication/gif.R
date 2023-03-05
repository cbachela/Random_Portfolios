  
  
  ############################################################################
  ### ANIMATED CHARTS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     02.03.2023
  # First version:    02.03.2023
  # --------------------------------------------------------------------------
  
  
  
  # https://www.youtube.com/watch?v=ikS1lhAzKWo
  
  
  
  
  # install.packages("gifski")  # for gif output
  # install.packages("av")      # for video output
  # install.packages("gganimate")
  
  require(timetk)
  require(dplyr)
  require(ggplot2)
  require(zoo)
  require(datasets)
  require(gganimate)
  # require(timeSeries)
  require(GPO)
  
  
  # X <- GPO::getMSCIData( universe = "dm" )
  X <- GPO::getMSCIData( universe = "bm" )
  X <- X[isWeekday(time(X)), ]
  X_level <- cumulated(X, "discrete")
  
  plot(X[ ,1])
  
  
  # df <- tk_tbl( X[ ,1] )
  # df <- transform( date = zoo:as.Date(index, frac = 0) )
  
  # dates <- rownames(X)[1:300]
  dates <- rownames( window(X, "2020-01-01", "2020-12-31") )
  # df <- as_tibble( tk_tbl( X[dates, 1] ) )
  df <- as_tibble( tk_tbl( X_level[dates, 1] ) )
  p <- df %>% 
        ggplot( aes( x = index, y = BM) ) + 
        geom_line() + 
        geom_point(size = 2) + 
        labs( x = "Date", y = "Percentage", title = "Returns" ) + 
        transition_reveal( index ) + 
        view_follow( fixed_y = TRUE )
  
  #// fps = frame per second
  animate(p, nframes = nrow(df), fps = 10 ) #, width = 800, height = 500 )  
  
  
  # anime_save( "filename.gif", p )
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Example from https://www.codingfinance.com/post/2020-03-20-create-animations-in-r/
  # --------------------------------------------------------------------------
  
  
  
  # install.packages( c("tidyquant", "tidyverse") )
  library(tidyquant)
  library(tidyverse)
  price_df <- tq_get(c('DIA', 'SPY', 'QQQ', 'IWM'),
                     from = "2016-11-1",
                     get = "stock.prices")
  
  ret_df <- price_df %>%
    group_by(symbol) %>%
    tq_transmute(select = adjusted,
                 mutate_fun = periodReturn,
                 period = "monthly",
                 col_rename = 'ret') %>%
    mutate(ret = if_else(row_number() == 1, 0, ret)) %>%
    mutate(cr = cumprod(1 + ret) - 1)
  
  
  # # Static
  # price_df %>%
  #   group_by(symbol) %>%
  #   tq_transmute(select = adjusted,
  #                mutate_fun = periodReturn,
  #                period = "monthly",
  #                col_rename = 'ret') %>%
  #   mutate(ret = if_else(row_number() == 1, 0, ret)) %>%
  #   mutate(cr = cumprod(1 + ret) - 1) %>%
  #   ggplot(aes(x = date, y = cr, group = symbol, color = symbol)) +
  #   geom_line() +
  #   geom_point() +
  #   geom_point(size = 2) + 
  #   scale_y_continuous(breaks = seq(-0.35,1, 0.1),
  #                      labels = scales::percent) +
  #   coord_cartesian(clip = 'off') + 
  #   labs(title = 'Major Index Returns since Trump\'s Elections in 2016', y = 'Returns (%)') + 
  #   theme_minimal() 
  
    
   
    
    p <- ret_df %>%
      ggplot(aes(x = date, y = cr, group = symbol)) +
      geom_line() +
      geom_segment(aes(xend = ymd("2020-03-25"), yend = cr), linetype = 2, colour = 'grey') + 
      geom_point(size = 2) + 
      geom_text(aes(x = ymd("2020-03-26"), label = symbol), hjust = 0) +
      scale_y_continuous(breaks = seq(-0.35,1, 0.1),
                         labels = scales::percent) +
      transition_reveal(date) +
      coord_cartesian(clip = 'off') + 
      labs(title = 'Major Index Returns since Trump\'s Elections in 2016', x = "Date", y = 'Returns (%)') + 
      theme_minimal() 
    
    animate(p, nframe = 200, end_pause = 20)
      
    
    
  
  
  
    
    
    
    # --------------------------------------------------------------------------
    # 2D density
    # --------------------------------------------------------------------------
    
    require(timetk)
    require(dplyr)
    require(ggplot2)
    require(zoo)
    require(datasets)
    require(gganimate)
    require(plotly)
    require(MASS)
    require(GPO)
    
    
    X <- GPO::getMSCIData( universe = "dm" )
    X_m <- aggMonthly(X)
    dates_m <- rownames(X_m)[-c(1:12)]
    
    # Rolling mean-variance
    mu <- applyRoll( Data = X, 
                     Width = 252,
                     charvec = dates_m,
                     FUN = function(X) { meanGeo(X, scalefactor = 252)} )
    sigma <- applyRoll( Data = X, 
                        Width = 252,
                        charvec = dates_m,
                        FUN = function(X) { apply(X, 2, sd ) * sqrt(252) } )                  
    
    
 
    data <- data.frame( x = as.numeric(sigma[1, ]),
                        y = as.numeric(mu[1, ]) )
    
    data %>% 
      ggplot( aes( x = x, y = y ) ) + 
      geom_density_2d( )
    
    
    data %>% 
      ggplot( aes( x = x, y = y ) ) + 
      geom_density_2d_filled( contour_var = "density" ) + 
      geom_density_2d( color = "orange" )
    
    
    
    
    
    L <- list()
    for ( i in 1:nrow(mu) ) {
      today <- zoo::as.Date(rownames(mu[i, ]))
      L[[i]] <- data.frame( x = as.numeric(sigma[today, ]),
                            y = as.numeric(mu[today, ]), 
                            group = rep(today, ncol(mu)) )
    }
    data <- do.call( rbind, L )
    str(data)
    
    
    
    data %>% 
      ggplot( aes( x = x, y = y ) ) + 
      geom_density_2d( )
    
    
    
    
    df <- as_tibble( tk_tbl( data ) )
    p <- df %>% 
      ggplot( aes( x = x, y = y) ) + 
      geom_density_2d( ) + 
      transition_reveal( group )
    
    #// fps = frame per second
    animate(p, nframes = 100, fps = 10 ) #, width = 800, height = 500 ) 
  

    
    # Check out
    # https://stackoverflow.com/questions/74824236/geom-density-2d-filled-and-gganimate-cumulative-2d-density-estimate-animation-o
    
    
    
    
    
    
    
    