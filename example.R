     
     # initialize and call neighbor identification algorithm
     pkgs <- c('raster','sp','sf','magrittr','geosphere')
     lapply(pkgs, require, character.only = T)
     source('f_sfdalgo.R')
     
     # use continental us counties as toy example
     spatial_df <- st_as_sf(getData('GADM', country='USA', level=2))
     spatial_df <- spatial_df[-grep('^US\\.(AK|HI)', spatial_df$HASC_2),]
     spatial_df <- spatial_df[-which(spatial_df$TYPE_2=='Water body')  ,]     
     
     # construct example variables
     spatial_df$income <- sample(10000:120000, nrow(spatial_df))
     spatial_df$temp   <- rnorm(nrow(spatial_df), 18, 5)
     spatial_df$prec   <- runif(nrow(spatial_df), 10, 40)
     spatial_df$lf     <- abs(round(rnorm(nrow(spatial_df), 5e5, 1e6)))
     
     # Spatial First Difference Transformation
     transformed <- SFD_vars2(spatial_df,  tol = 0,
                              dependent_var    = 'income',
                              independent_vars = c('temp','prec','lf'),
                              rotation = seq(0, 180, 30),
                              needs_balanced = FALSE, balanced_df = NA, 
                              plot = TRUE, "sfdusa.png")
     
     # Inefficiency measure: as sfd, like fe, removes singletons,
     # the more singletons per sequence the more inefficient
     eff.measure <- transformed %>% group_by(angle, group) %>%
          summarise(nobs = n()) %>%
          ungroup() %>% group_by(angle) %>%
          summarise(e = length(which(nobs==1))/sum(nobs))
     