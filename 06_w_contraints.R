source("03_prep_simann.R")



# cc
cc <- vector("integer", length(S))
b1 <- 1347
b2 <- 2342
cc[1:b1] <- 1
cc[(b1 + 1):b2] <- 2
cc[(b2 + 1):2547] <- 3
cc <- as.integer(cc)

ggplot() +
  geom_sf(data = state_map_sf, fill = "white", color = "black") +
  #geom_sf(data = population_points_sf, color = "forestgreen", size = 0.5) +
  #geom_sf(data = highway_lines_sf, color = "grey", size = 0.1, alpha = 0.5) +
  geom_sf(data = filtered_points[which(cc == 1), ], color = "red", size = 2) +
  geom_sf(data = filtered_points[which(cc == 2), ], color = "blue", size = 2) +
  geom_sf(data = filtered_points[which(cc == 3), ], color = "green", size = 2) +
  #geom_sf(data = final_buffers, color = "red", size = 5, fill = NA) +
  #geom_sf(data = population_points_sf[pop_pts, ], color = "orange", size = 1) +
  #geom_sf(data = buffers[673,], color = 'red', fill = NA) +
  theme_minimal() +
  labs(title = "Population Points and Major Highways in Selected Northeastern States")

# nzones
nzones = as.integer(3)

# you need to initialize one in each starting zone
S <- vector("integer", length(S))
S[1:length(S)] <- as.integer(0)
S[sample(1:b1, 2)] <- 1
S[sample((b1+1):b2, 2)] <- 1
S[sample((b2+1):2547, 1)] <- 1
which(S == 1)

save(list = ls(), file = "test_objects.Rdata")


#

load("test_objects.Rdata")

# obviously takes slower
system("R CMD SHLIB simann_constrain.f90")
#dyn.unload("simann_constrain.so")
dyn.load("simann_constrain.so")

oo <- .Fortran("simann_constrain",
               S = as.integer(S),
               magic_n = as.integer(length(which(S == 1))),
               np = as.integer(nrow(site_pairs_sub)),
               nsites = as.integer(length(unique(site_pairs_sub$site_id))),
               site_pairs_sub = aa,
               row_lookup = bb,
               penalty = 1,
               SCORE = 0.,
               cooling_rate = 0.95,
               cc = cc,
               nzones = nzones,
               verbose = as.integer(1)
               ) 

# confirming math
which(S == 1)
which(oo$S == 1)
library(tidyverse)
library(sf)
get_score(oo$S)
oo$SCORE

pop_pts <- site_pairs_sub %>% filter(site_id %in% which(oo$S == 1)) %>% pull(pop_id)

final_buffers <- st_buffer(filtered_points[which(oo$S == 1), ], RADIUS <- 20 * 1609.34)
# plot(final_buffers)

ggplot() +
  geom_sf(data = state_map_sf, fill = "white", color = "black") +
  geom_sf(data = population_points_sf, color = "forestgreen", size = 0.5) +
  geom_sf(data = highway_lines_sf, color = "grey", size = 0.1, alpha = 0.5) +
  geom_sf(data = filtered_points[which(oo$S == 1), ], color = "red", size = 2) +
  geom_sf(data = final_buffers, color = "red", size = 5, fill = NA) +
  geom_sf(data = population_points_sf[pop_pts, ], color = "orange", size = 1) +
  #geom_sf(data = buffers[673,], color = 'red', fill = NA) +
  theme_minimal() +
  labs(title = "Population Points and Major Highways in Selected Northeastern States")
