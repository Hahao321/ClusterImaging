change_resolution <- function(image_df, x_size)
{
  ## change_resolution(image_df, x_size) subsamples an image to produce
  ## a lower resolution image. Any non-coordinate columns in the data
  ## frame are summarized with their most common value in the larger
  ## grid cell.
  ##
  ## Input:
  ## - image_df: A data frame in wide format. The x-coordinate column MUST
  ##             be named 'x' and the y-coordinate column MUST be named 'y'.
  ##             Further columns have no naming restrictions.
  ## - x_size:   The number of cells in the x-direction. The number of cells
  ##             in the vertical direction will be computed to maintain the 
  ##             perspective. There is no guarantee that the exact number
  ##             of cells in the x-direction is x_size
  ##
  ## Output:
  ## - A data frame with the same column names as image_df, but with fewer 
  ##   entries that corresponds to the reduced resolution image.
  ##
  ## Example:
  ##   library(imager)
  ##   library(dplyr)
  ##   fpath <- system.file('extdata/Leonardo_Birds.jpg',package='imager') 
  ##   im <- load.image(fpath)
  ##   im_dat<- as.data.frame(im,wide = "c") %>% rename(R = c.1, G = c.2, B = c.3) %>%
  ##            select(x,y,R,G,B)
  ##   agg_image <- change_resolution(im_dat, 50)
  
  if(!require(sp)) {
    stop("The sp packages must be installed. Run install.packages(\"sp\") and then try again.")
  }
  if(!require(dplyr)) {
    stop("The dplyr packages must be installed. Run install.packages(\"dplyr\") and then try again.")
  }
  
  sp_dat <- image_df 
  gridded(sp_dat) = ~x+y
  
  persp = (gridparameters(sp_dat)$cells.dim[2]/gridparameters(sp_dat)$cells.dim[1])
  y_size = floor(x_size*persp)
  orig_x_size = gridparameters(sp_dat)$cells.dim[1]
  orig_y_size = gridparameters(sp_dat)$cells.dim[2]
  
  x_res = ceiling(orig_x_size/x_size)
  y_res = ceiling(orig_y_size/y_size)
  
  gt = GridTopology(c(0.5,0.5), c(x_res, y_res),
                    c(floor(orig_x_size/x_res), floor(orig_y_size/y_res)))
  SG = SpatialGrid(gt)
  agg = aggregate(sp_dat, SG, function(x) names(which.max(table(x)))[1] )
  agg@grid@cellsize <- c(1,1)
  df <- agg %>% as.data.frame %>% rename(x = s1, y = s2)  %>% select(colnames(image_df))
  
  return(df)
  
}



process_image <- function(image_file_name, k_list){
  ## process_image(image_file_name, k_list) takes a image and a few values to create
  ## a k means cluster. Includes information on the clustering data, k mean centers
  ## as well as RGB/ DMC values.
  ##
  ## Input:
  ## - image_file_name: any PNG/JPG image
  ## - k_list: a list of integers to act as number of clusters.
  ##
  ## Output:
  ## - a list that includes k clustering data for each of the integers in k_list,
  ## dataframes with RGB/DMC values, and a dataframe in wide format.
  ##
  ## Example:
  ## im <- imager::load.image("C:/Users/haoyi/Desktop/MS.png")
  ## k_list = c(3, 5, 7)
  ##cluster_info <- process_image(im, k_list)
  tidy_dat <- as.data.frame(image_file_name, wide = "c") %>% rename(R = c.1, G = c.2, B = c.3)
  clustdata <- tibble(k_list) %>%
    mutate(
      kclust = map(k_list, ~kmeans(select(tidy_dat,c(-x,-y)), centers = .x, nstart = 6)),
      glanced = map(kclust, glance),
    )
  centres <- list()
  b = 1
  for(i in clustdata$kclust)
  {
    centres[[b]] <- tidy(i)
    centres[[b]] <- centres[[b]] %>% mutate(col = rgb(R,G,B))
    b = b+1
  }
  j = 1
  while(j < length(centres) + 1)
  {
    temp <- list()
    z <- 1
    for(n in centres[[j]]$col)
    {
      temp[z] <- dmc(n)$name
      z = z+1
    }
    dmc_n <- unlist(temp)
    centres[[j]]$name <- c(dmc_n)
    j = j + 1
  }
  
  data <- list(clustdata, centres, tidy_dat) 
}


scree_plot <- function(cluster_info){
  ## scree_plot(cluster_info) takes cluster_info, which is the output of process_image
  ## and uses the data to create a scree_plot which plots the values in k_list on x-axis
  ## and the corresponding tot.withinss(Total within-cluster sum of squares) value
  ##
  ## Input:
  ## - cluster_info, the list of dataframes obtained from process_image
  ##
  ## Output:
  ## - plots a scree_plot based on K-list and tot.withinss
  ##
  ## Example:
  ## im <- imager::load.image("C:/Users/haoyi/Desktop/MS.png")
  ## k_list = c(3, 5, 7)
  ## cluster_info <- process_image(im, k_list)
  ## scree_plot(cluster_info)
  within_list <- list()
  x = 1
  for(i in cluster_info[[1]]$glanced){
    within_list[x] <- i$tot.withinss
    x = x + 1
  }
  withinss <- unlist(within_list)
  
  ggplot(cluster_info[[1]], aes(k_list, withinss)) +
    ggtitle( "Scree-plot for the given values in k_list")+
    ylab("tot.withinss") + xlab("k")+
    geom_line() +
    geom_point() 
}


colour_strips <- function(cluster_info){
  ## colour_strips(cluster_info) takes cluster_info, which is the output of process_image
  ## and uses the RGB data within to create color strips
  ##
  ## Input:
  ## - cluster_info, the list of dataframes obtained from process_image
  ##
  ## Output:
  ## - creates a color strip for each value different clustering in k_list based on 
  ## RGB data
  ##
  ## Example: 
  ## im <- imager::load.image("C:/Users/haoyi/Desktop/MS.png")
  ## k_list = c(3, 5, 7)
  ## cluster_info <- process_image(im, k_list)
  ## colour_strips(cluster_info)
  for (i in cluster_info[[2]])
    show_col(i$col)
}


make_pattern <- function(cluster_info, k, x_size, black_white, background_colour = NULL){
  ## make_pattern(cluster_info, k, x_size, black_white, background_colour = NULL) 
  ## creates a cross stitch based on cluster_info, and a chosen clustering k(from k_list)
  ## changes the resolution based on x_size, which is the number of stitches in the 
  ## x-direction. The stitched pattern can be in color or black and white, based on the
  ## boolean black_white. Finally, the back_ground color can be removed from the stitch
  ## by inputting the hex string, by default this is NULL and the background will not be
  ## removed
  ##
  ## Input:
  ## - cluster_info: the output of process_image, holds the data on the RGB values of the
  ## original image as well as the clustered data
  ## - k: the chosen k value from k_list.
  ## - x_size: the number of desired stitches in the x-direction
  ## - black_white: a boolean that determines if the stitch will include color or not.
  ## - background_colour: removes the color from the stitch based on an inputted RGB
  ##   hex code, by default this is NULL.
  ##
  ## Output:
  ## - A cross stitch with a legend that includes a legend based on what color to stitch.
  ## The image will be based on the parameters above
  ##
  ## Example:
  ## im <- imager::load.image("C:/Users/haoyi/Desktop/MS.png")
  ## k_list = c(3, 5, 7)
  ## cluster_info <- process_image(im, k_list)
  ## make_pattern(cluster_info, 7, 40, TRUE)
  g = 1
  while(g < length(cluster_info[[1]]) + 1)
  {
    if(k == cluster_info[[1]]$k_list[[g]])
      img <- augment(cluster_info[[1]]$kclust[[g]], cluster_info[[3]]) %>% rename(cluster= .cluster)
    sub <- cluster_info[[2]][[g]]
    g=g+1
  }
  new_res <- change_resolution(img, x_size)
  
  if(is.null(background_colour) == FALSE){
    sub$cluster[sub$col == background_colour] <- NA
    sub$col[sub$col == background_colour] <- NA
    print(sub)
    print(cluster_info[[2]])
  }
  
  if(black_white == FALSE){
    ggplot(new_res, aes(x=x, y = y), color = factor(cluster), shape = factor(cluster)) +
      geom_point(aes(col = factor(cluster), shape = factor(cluster))) +
      scale_colour_manual(name = "Color to Stitch",
                          values = sub %>% select(cluster, col) %>% deframe,
                          label =  sub %>% select(cluster, name) %>% deframe)+
      scale_shape_manual(name = "Color to Stitch",
                         values = sub %>% select(cluster) %>% deframe,
                         label = sub %>% select(cluster, name) %>%deframe)+
      scale_y_reverse() + theme_void() +
      background_grid()
  }
  
  
  else{
    ggplot(new_res, aes(x=x, y = y)) +
      geom_point(aes(shape = factor(cluster))) +
      labs(shape = "Color to Stitch") +
      scale_colour_manual(name = "Color to Stitch",
                          values = sub %>% select(cluster, col) %>% deframe,
                          label =  sub %>% select(cluster, name) %>% deframe)+
      scale_shape_manual(name = "Color to Stitch",
                         values = sub %>% select(cluster) %>% deframe,
                         label =  sub %>% select(cluster, name) %>%deframe)+
      scale_y_reverse() + theme_void() +
      background_grid()
  }
}