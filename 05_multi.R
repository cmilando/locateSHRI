library(tidyverse)
library(sf)
library(pbapply)
load("all_objects.RData")

RADIUS <- 20 * 1609.34
potential_sites <- st_coordinates(filtered_points)
population_dots <- st_coordinates(st_transform(population_points_sf, crs = 3857))

head(potential_sites)
head(population_dots)

site_pairs <- expand_grid(site_id = 1:nrow(potential_sites),
                          pop_id = 1:nrow(population_dots))

site_dist <- pbsapply(1:nrow(site_pairs), function(i) {
  site_i = site_pairs$site_id[i]
  pop_i = site_pairs$pop_id[i]
  sqrt((potential_sites[site_i,1] - population_dots[pop_i, 1])^2 +
         (potential_sites[site_i,2] - population_dots[pop_i, 2])^2)
})

# so first you can subset to only ones within the radius
site_pairs$dist <- site_dist
dim(site_pairs)

site_pairs_sub <- site_pairs %>% filter(site_dist < RADIUS)
dim(site_pairs_sub)

# yeah this is now only 93,920 sites-pop sites
length(unique(site_pairs_sub$site_id))
nrow(potential_sites)

# all these variables should be range01 for additivity
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

site_pairs_sub$population <- 0.5

metric_df <- data.frame(pop_id = unique(site_pairs_sub$pop_id))

set.seed(123)
metric_df$metric1 <- range01(rnorm(1:nrow(metric_df)))
site_pairs_sub <- left_join(site_pairs_sub, metric_df)

## this is the function to run over
beta_weights_all <- as.matrix(expand_grid(seq(0.1, 5, 1), seq(0.1, 5, 1)))
nrow(beta_weights_all)

#

get_est <- function(beta_i) {
  
  # get new betas
  beta_weights <- as.vector(beta_weights_all[beta_i,])
  
  score = as.matrix(site_pairs_sub[, c('population', 'metric1')]) %*% beta_weights
  
  site_pairs_sub$score = as.vector(score)
  
  head(site_pairs_sub)
  
  # make sure its arranged
  stopifnot(site_pairs_sub$site_id == sort(site_pairs_sub$site_id))
  
  rlepX <- rle(site_pairs_sub$site_id)
  
  # Compute the indices of the start and end of each run
  end <- cumsum(rlepX$lengths)
  start <- c(1, lag(end)[-1] + 1)
  
  xx <- data.frame(site_id = unique(site_pairs_sub$site_id), start, end)
  head(xx)
  
  #### WHAT ARE YOU PASSING IN
  head(xx)
  
  head(site_pairs_sub)
  
  S <- rep(0, length(unique(site_pairs_sub$site_id)))
  S[c(1:5)] <- 1
  
  penalty <- 1
  
  ## make aa
  aa <- as.matrix(site_pairs_sub, 
                  nrow = nrow(site_pairs_sub),
                  ncol = 6)
  head(aa)
  dim(aa)
  
  ## make bb
  bb = as.matrix(xx, 
                 nrow = nrow(xx),
                 ncol = ncol(xx))
  head(bb)
  dim(bb)
  
  ##
  dyn.load("simann.so")
  oo <- .Fortran("simann",
                 S = as.integer(S),
                 magic_n = as.integer(length(which(S == 1))),
                 np = as.integer(nrow(site_pairs_sub)),
                 nsites = as.integer(length(unique(site_pairs_sub$site_id))),
                 site_pairs_sub = aa,
                 row_lookup = bb,
                 penalty = 1,
                 SCORE = 0.,
                 cooling_rate = 0.95,
                 verbose = as.integer(0)) 
  
  return(oo$S)
  
  
}

library(future)
library(future.apply)
library(pbapply)
plan(multisession)

out_mat_l <- pblapply(1:nrow(beta_weights_all), get_est, cl = 'future')
out_mat <- do.call(rbind, out_mat_l)

dim(out_mat)

out_mat_sum <- colSums(out_mat) / nrow(out_mat)

##
load("all_objects.RData")
filtered_points 

filtered_points_df <- st_as_sf(filtered_points)

filtered_points_df$percent <- ifelse(out_mat_sum == 0., NA, out_mat_sum)

filtered_points_df <- filtered_points_df %>% filter(!is.na(percent))

filtered_points_metric <- st_as_sf(filtered_points_df)

ggplot() +
  geom_sf(data = state_map_sf, fill = "white", color = "black") +
  #geom_sf(data = population_points_sf, color = "forestgreen", size = 0.5) +
  #geom_sf(data = highway_lines_sf, color = "orange", size = 0.1, alpha = 0.5) +
  geom_sf(data = filtered_points, size = 2, shape = 1, color = 'grey',
          alpha = 0.5) +
  geom_sf(data = filtered_points_metric, 
          mapping = aes(color = percent), size = 4, shape = 16) +
  theme_minimal() +
  scale_color_viridis_c(option = 'A', direction = -1) +
  labs(title = "Population Points and Major Highways in Selected Northeastern States")

