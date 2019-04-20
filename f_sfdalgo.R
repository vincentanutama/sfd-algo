


     # Original algorithm: Hannah Druckenmiller, March 11 2019
     # Some updates by Vincent Tanutama, April 20 2019
     # To do:
     # 1. Update algorithm when needs_balanced = T
     
     SFD_vars2 <- function (spatial_df, tol = 0, dependent_var, independent_vars, 
               rotation = 0, needs_balanced = FALSE, balanced_df = NA, 
               plot = TRUE, pngfilename) 
     {
          print(noquote("Transforming variables..."))
          library(magrittr, quietly = "true")
          library(dplyr, quietly = "true")
          library(sp, quietly = "true")
          library(sf, quietly = "true")
          library(geosphere, quietly = "true")
          
          if (!any(class(spatial_df) %in% c("SpatialPolygonsDataFrame", "sf"))) 
               stop("spatial_df must be a SpatialPolygonsDataFrame")
          if (!class(tol) == "numeric") 
               stop("tol must be a number between zero and one")
          if (!class(dependent_var) == "character") 
               stop("dependent_var must be a character")
          if (!class(independent_vars) == "character") 
               stop("independent_vars must be a character or a list of characters (e.g. independent_vars=c(\"age\", \"gender\"))")
          if (!class(needs_balanced) == "logical") 
               stop("needs_balanced must be a logical")
          if (needs_balanced == TRUE & !class(balanced_df) == "data.frame") 
               stop("balanced_df must be a data.frame containing the obs_var")
          if (!class(plot) == "logical") 
               stop("plot must be a logical")
          if (plot & grep('png|PNG', pngfilename)!=1) 
               stop("plot filename needs to be supplied in png (e.g 'DIRECTORY/sfd.png')")
          if (!class(rotation) == "numeric") 
               stop("rotation must be a numeric vector")
          
          # Rotation function
          rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
          
          # Results holder
          robSDF <- data.frame()
          SP <- spatial_df # SP is now the original spatial_df
          SP$rowno <- 1:nrow(SP)
          
     # Here I defined rotation to take a numeric vector 
     # in degrees, e.g. seq(15, 180, 15)
     for (dg in rotation) {
   
          # Convert sp to sf object for more stablity
          spatial_df <- st_as_sf(SP) 
          st_crs(spatial_df) <- 4326
          # Identify row index of polygon i's neighbors
          print(noquote(paste0('Identifiying neighbors for rotation ', dg,' degree...')))
          neigh <- st_touches(spatial_df, spatial_df, sparse = T)
          # Produce centroid of polygon i for bearing calculation
          cent  <- lapply(st_centroid(st_geometry(spatial_df)), as.vector)
          
          
          # Rotate shapefile
          if (dg!=0) {
               print(noquote(paste0('Rotating polygons by ', dg,' degree...')))
               spdf.geom <- st_geometry(spatial_df) 
               # need to rotate object that contains no crs information because
               # rotation can produce points outside the defined extent of the globe, 
               # i.e. bigger than [-180,180]x[-90,90]
               spdf2 <- spdf.geom * rot(dg*pi/180)
               spatial_df <- st_sf(cbind(as.data.frame(spatial_df),spdf2))
          }     
          
          # Compare centroid-based relative position of each polygon and its neighbors
          pone <- Map(rep, as.list(1:nrow(spatial_df)), lapply(neigh,length))
          p1 <- lapply(pone, as.list)
          p2 <- lapply(neigh, as.list)
          p1 <- unlist(unlist(rapply(p1, function(x) cent[x], how = 'replace'), recursive = F), recursive = F)
          p2 <- unlist(unlist(rapply(p2, function(x) cent[x], how = 'replace'), recursive = F), recursive = F)
          bearing12 <- unlist(Map(bearing, p2, p1)) + 180 - 22.5
          e <- which(sapply(pone, length)==0)
          if  (length(e)!=0) { 
               pone[e]  <- as.list(e)
               for (k in 1:length(e)) {
                    if (e[k]==1) {
                         bearing12 <- c(NA, bearing12)
                    } else {
                         bearing12 <- c(bearing12[1:length(unlist(pone[1:(e[k]-1)]))], NA,
                                        bearing12[length(unlist(pone[1:(e[k])])):length(bearing12)])
                    }
               }
          } 
          
               # Start sampling neighbors starting at leftmost polygon,
               # going in northeast direction until last polygon no longer has neighbors to sample
               # Algorithm is greedy in sampling multiple neighbors at one time
               print(noquote('Sampling counties'))
               left    <- st_bbox(spatial_df)[1]
               bb      <- lapply(st_geometry(spatial_df), st_bbox)
               which(sapply(sapply(bb, function(x) which(x == left)), length) != 0) -> polyleft
               pb1 <- txtProgressBar(min = 1, max = round(nrow(spatial_df)*(1-tol)), initial = 1, style = 3)
               
               # put bearing vector back to a list object per polygon i
               # so each element is the bearings of neighbors of polygon i
               bears <- split(bearing12, unlist(pone)) 
               universe <- 1:nrow(spatial_df)
               # panel is the list whose element is a vector of ordered sampled neighbors
               panel <- list()
               # chopped_map is the shapefile after neighbors have been sampled
               chopped_map <- spatial_df
               chopped_map$no <- 1:nrow(chopped_map)

               
               # Can set tolerance to non-zero, but that will leave some counties unsampled :(
               while (nrow(chopped_map)>(nrow(spatial_df)*tol)) {
  
                    setTxtProgressBar(pb1, length(unlist(panel)))
                    
                    p <- polyleft
                    sample <- polyleft
                    m <- 1
                    
                    while(length(polyleft)!=0) {
                         
                         # find p's next neighbor (.r = row number, .d = bearing)
                         nx.r <- neigh[[p]] 
                         nx.d <- bears[[p]]  
                         
                         # sort bearings of neighbors in the clockwise order starting from north
                         sortma  <- sort(nx.d)
                         # when there are multiple neighbors ...>
                         if (length(sortma)!=1) {
                         #>... grab all neighbors -- this is what I meant by ``greedy'' --
                         #     in the counter-clockwise direction
                         #     so that at the end of sampling, we're back to the upperleftmost polygon
                         #     but want to stop greedy neighbor extraction if sequence of neighbor isn't contiguous
                         #     so need to check for 2nd degree neighbors too
                              m <- rev(match(sortma, nx.d))
                              nm <- lapply(m, function(x) neigh[[nx.r[x]]])     # 2nd degree neighbors of each neighbor of i
                              checknm <- Map(match, as.list(lead(nx.r[m])), nm) # check for 2nd degree neighbors in the set of 2nd degree neighbors
                              checknm <- unlist(checknm)[-length(checknm)]      # NA means in the set of neighbors of i = {A, B, C, D}, B isn't neighbor of A
                              wnm <- which(is.na(checknm))
                              if (length(wnm)!=0) {
                                  
                                   if (max(wnm) %in% c(1,length(checknm))) {
                                        m <- match(sortma, nx.d)[[1]]
                                   } else {
                                        m0 <- max(wnm)
                                        m  <- m[(m0+1):length(checknm)]
                                   }
                                   
                              } 
                              
                         #>... otherwise just grab that neighbor
                         } else {
                              m <- which(nx.d == sortma)
                         }
                         
                         # remove sampled polygons from neigh, bears, universe 
                         # add sampled polygons to object sample
                         if(length(m)!=0) {
                              
                              nxneigh <- nx.r[m]
                              wn <- which(neigh[[p]]%in%nxneigh)
                              universe[c(p, nxneigh)] <- NA
                              ns <- lapply(neigh, function(x) which(x %in% c(p,nxneigh)))
                              wns <- which(sapply(ns, length)!=0)
                              neigh[wns] <- Map(function(x, y) replace(x, y, NA),
                                                neigh[wns], ns[wns])
                              bears[wns]  <- Map(function(x, y) replace(x, y, NA),
                                                 bears[wns], ns[wns])
                              sample <- c(sample, nxneigh)
                              
                         } else {
                              
                              # when there is no more polygons to be sampled, break the while loop
                              # and end this current sequence of neighbor sampling ...(*)
                              polyleft <- integer(0)
                              nxneigh  <- integer(0)
                         }
                         
                         # define polygon for next iteration
                         # nxneigh is all neighbors in current iteration
                         # p is element of nxneigh for next iteration
                         if (length(nxneigh)>1) {
                              p <- nxneigh[length(nxneigh)]
                         } else if (length(nxneigh)==1) {
                              p <- nxneigh
                         } else {
                              p <- integer(0)
                         }
                         
                    }
                    
                    # (*)... redefine next sequence of neighbors
                    chopped_map <- chopped_map[-which(chopped_map$no %in% sample),]
                    left    <- st_bbox(chopped_map)[1]
                    which(sapply(sapply(bb, function(x) which(x == left)), length) != 0) -> polyleft
                    if (length(polyleft)>1) polyleft <- polyleft[which(!polyleft%in%unlist(panel))][[1]]
                    panel <- append(panel, list(spatial_df$rowno[sample]))
                    
               }
               
               # Convert list of sequences of neighbors into dataframe
               # that contains sequence number and ordered polygons within the sequence
               panelgroup <- Map(rep, as.list(1:length(panel)), lapply(panel, length))
               paneldf <- data.frame(group = unlist(panelgroup),
                                     rowno = unlist(panel),
                                     stringsAsFactors = F)
               paneldf$angle <- dg
               # since ordering is paramount, need to have index of ordered rows
               paneldf$index <- 1:nrow(paneldf) 
               robSDF <- rbind(robSDF, paneldf)
               
     }
      
          # Produce transformed variables
               dvars <- setNames(as.data.frame(as.data.frame(SP)[robSDF$rowno,dependent_var]),
                                 dependent_var)
          if (length(independent_vars)==1) {
               ivars <- setNames(as.data.frame(as.data.frame(SP)[robSDF$rowno,independent_vars]),
                                 independent_vars)
          } else {
               ivars <- as.data.frame(SP)[robSDF$rowno,independent_vars]
          }
          robSDF <- cbind(robSDF, dvars, ivars)
          robSDF <- robSDF[,c('angle','group','rowno','index',dependent_var,independent_vars)]
          rowno  <- robSDF[,"rowno"]
          robSDF <- robSDF[,-grep('rowno',colnames(robSDF))]
          minuslag <- function(x,n) x - dplyr::lag(x, n)
          A <- robSDF %>%
               group_by(angle, group) %>%
               arrange(angle, group, index) %>%
               mutate_all(minuslag, 1)
          B <- A %>%
               mutate_all(minuslag, 1)
          A <- setNames(as.data.frame(A[,c(dependent_var,independent_vars)]),
                        paste0('sfd.', c(dependent_var,independent_vars)))
          B <- setNames(as.data.frame(B[,c(dependent_var,independent_vars)]),
                        paste0('sdd.', c(dependent_var,independent_vars)))
          robSDF$rowno <- rowno
          robSDF <- cbind(robSDF, A, B)
          rownames(robSDF) <- 1:nrow(robSDF)
              
          # Produce plot
          if (plot) {
               
               print(noquote(paste0('Producing plot ', pngfilename)))
               pb2 <- txtProgressBar(min = 1, max = length(rotation), initial = 1, style = 3)

               # Determine size of plot
               fact <- function(x) {
                    x <- as.integer(x)
                    div <- seq_len(abs(x))
                    factors <- div[x %% div == 0L]
                    return(factors)
               }
               
               is.prime <- function(x) {
                    x <- as.integer(x)
                    div <- seq_len(abs(x))
                    factors <- div[x %% div == 0L]
                    return(all(factors %in% c(1,x)))
               }
               
               s <- N <- length(rotation)
               boxsize <- st_bbox(SP)
               box.ratio <- unname(
                            (boxsize[3]-boxsize[1])/
                            (boxsize[4]-boxsize[2]))
               while(is.prime(s)&s>=3) {
                    s <- s+1
               }
                    fa <- fact(s)
                    ro <- fa[floor((length(fa)+1)/2)]
                    matsize <- sort(c(ro, s/ro))
                    
                    # plot
                    colors <- paste0('grey',seq(10,70,20))
                    png(pngfilename, width = 1980, 
                        height = 1980/(box.ratio*matsize[2]/matsize[1]))
                    par(mar=rep(0,4),mfrow=matsize)
                    for (i in 1:N) {
                         
                         setTxtProgressBar(pb2, i)
                         
                         plot(st_geometry(SP),col=NA,border=NA,
                              main=paste0(rotation[i],' degree rotation'))
                         paneldf <- robSDF[robSDF$angle==rotation[i],]
                         paneldf <- paneldf[order(paneldf$index),]
                         for (j in 1:max(paneldf$group)) {
                              plot(st_geometry(SP[paneldf$rowno[paneldf$group==j],]),
                                   col=colors[j%%length(colors)+1], border=NA, add=T)
                         }
                    }
                    dev.off()
          }
          
          return(robSDF)
     
     }