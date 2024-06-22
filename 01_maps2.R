library(tidyverse)
library(maps)
library(sf)
library(tigris)          # fips_codes
library(mapproj)
library(lwgeom)
library(pbapply)

# Define the populations
populations <- data.frame(
  state_code = c("25", "44", "09"), # MA, RI, CT
  population = c(15 * 10^6, 1 * 10^6, 3 * 10^6)
)

# Load state boundaries and filter for the NE states
NE <- c("25", "44", "09")
state_map <- map_data("state")
state_fips_codes <- unique(fips_codes[, c('state_name', 'state_code')]) %>%
  rename(region = state_name) %>%
  mutate(region = tolower(region),
         state_code = as.character(state_code)) %>%
  filter(state_code %in% NE)
state_map <- inner_join(state_map, state_fips_codes, by = "region")
state_map <- state_map %>%
  filter(!(state_code == '25' & group != 21))

# Create a proper polygon structure from the state_map data
state_map <- state_map %>%
  group_by(region, state_code, group) %>%
  summarize(geometry = st_sfc(st_polygon(list(cbind(long, lat))))) %>%
  ungroup()

# Convert to SF object
state_map_sf <- st_as_sf(state_map, crs = 4326)

# Function to generate regular grid of points within a state polygon
generate_points <- function(state_polygon, num_points) {
  bbox <- st_bbox(state_polygon)
  grid_size <- ceiling(sqrt(num_points))
  x_range <- seq(bbox[1], bbox[3], length.out = grid_size)
  y_range <- seq(bbox[2], bbox[4], length.out = grid_size)
  grid <- expand.grid(x = x_range, y = y_range)
  points <- st_as_sf(grid, coords = c("x", "y"), crs = 4326)
  points <- points[st_within(points, state_polygon, sparse = FALSE)[,1], ]
  return(points)
}

# Generate points for each state
population_points <- lapply(1:nrow(state_map_sf), function(i) {
  this_state_code <- state_map_sf$state_code[i]
  state_polygon <- state_map_sf[i, ]
  pop <- populations %>% filter(state_code == this_state_code) %>% pull(population)
  num_points <- pop / 10000
  num_points
  points <- generate_points(state_polygon, num_points)
  points$state_code <- this_state_code
  return(points)
})

# Combine all points into a single SF object
population_points_sf <- do.call(rbind, population_points)

# Plot the state boundaries, population points, and filtered highways
ggplot() +
  geom_sf(data = state_map_sf, fill = "white", color = "black") +
  geom_sf(data = population_points_sf, color = "forestgreen", size = 0.5) +
  theme_minimal() +
  labs(title = "Population Points in Selected Northeastern States")


# Download major highways data using TIGRIS
highways <- primary_secondary_roads(state = "MA") %>% 
  bind_rows(primary_secondary_roads(state = "CT")) %>%
  bind_rows(primary_secondary_roads(state = "RI"))

# Convert highways to points
highway_lines_sf <- st_cast(highways, "LINESTRING")

ggplot() +
  geom_sf(data = state_map_sf, fill = "white", color = "black") +
  geom_sf(data = population_points_sf, color = "forestgreen", size = 0.5) +
  geom_sf(data = highway_lines_sf, color = "blue", size = 0.1, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Population Points in Selected Northeastern States")

# ##############################################
#st_length(highways_transformed)
#max(st_length(highways_transformed))
highways_transformed <- st_transform(highway_lines_sf, crs = 3857)

MILES <- 10

xby = MILES * 1609.34 / max(st_length(highways_transformed))
#xby

highways_sampled <- st_line_sample(highways_transformed, 
                                   sample = units::set_units(seq(0, 1, as.numeric(xby))))

highways_sampled_cast <- st_cast(highways_sampled, "POINT")
length(highways_sampled_cast)

highways_sampled_cast_df <- st_coordinates(highways_sampled_cast)
head(highways_sampled_cast_df)

# Remove points within a 1-mile buffer of each other
buffer_distance <- units::set_units(1, "mile")

# Create buffers around each point
buffers <- st_buffer(highways_sampled_cast, buffer_distance)

ggplot() +
  geom_sf(data = state_map_sf, fill = "white", color = "black") +
  geom_sf(data = population_points_sf, color = "forestgreen", size = 0.5) +
  geom_sf(data = buffers, color = "blue", size = 0.1, alpha = 0.5,
          fill = NA) +
  theme_minimal() +
  labs(title = "Population Points in Selected Northeastern States")


# Function to select the maximum non-overlapping subset of buffers
# this takes ~ 45 minutes for 63,144 nodes

# wait you don't need to do buffers, just do distance

# or you could do st_union

# ********
# -- YOU CAN MAKE THIS FASTER, THIS DOUBLE COUNTS
# ********
# select_non_overlapping_buffers <- function(buffers) {
#   selected_indices <- c()
#   for (i in seq_len(length(buffers))) {
#     if(i && 1000 == 0) timestamp(suffix = paste(" >", i))
#     if (length(selected_indices) == 0 || 
#         all(!st_intersects(buffers[i], buffers[selected_indices], sparse = FALSE))) {
#       selected_indices <- c(selected_indices, i)
#     }
#   }
#   return(selected_indices)
# }
# 
# select_non_overlapping_buffers2 <- function(buffers) {
#   selected_indices <- c()
#   for (i in seq_len(length(buffers))) {
#     if(i && 1000 == 0) timestamp(suffix = paste(" >", i))
#     if (length(selected_indices) == 0 || 
#         !st_intersects(buffers[i], super_buffer, sparse = FALSE)) {
#       selected_indices <- c(selected_indices, i)
#       if(i == 1) {
#         super_buffer = buffers[i]
#       } else {
#         super_buffer <- st_union(super_buffer, buffers[i])
#       }
#     }
#   }
#   return(selected_indices)
# }
# 
# library(pbapply)
# 
# select_non_overlapping_buffers3 <- function(pts) {
#   selected_indices <- c()
#   print(nrow(pts))
#   for (i in seq_len(nrow(pts))) {
#     print(i)
#     if(i && 1000 == 0) timestamp(suffix = paste(" >", i))
#     length0 <- length(selected_indices) == 0
#     if(length(selected_indices) == 0) {
#       within_radius = F
#     } else {
#       # buffers[i], buffers[selected_indices]
#       within_radius <- sapply(selected_indices, function(s_idx) {
#         sqrt((highways_sampled_cast_df[i,1] - highways_sampled_cast_df[s_idx, 1])^2 +
#                (highways_sampled_cast_df[i,2] - highways_sampled_cast_df[s_idx, 2])^2) 
#       })
#       print(within_radius)
#     }
#     if (length0 || all(! within_radius <= 1609.34)) selected_indices <- c(selected_indices, i)
#   }
#   return(selected_indices)
# }

# Select non-overlapping buffers
# doesn't take too long, maybe an hour
# selected_indices <- select_non_overlapping_buffers(buffers)
# selected_indices <- select_non_overlapping_buffers2(buffers)
# selected_indices <- select_non_overlapping_buffers3(highways_sampled_cast_df[1:2,])
# selected_indices

# system("R CMD SHLIB get_buffers.f90")
# dyn.unload("get_buffers.so")
dyn.load("get_buffers.so")
oo <- .Fortran("select_non_overlapping_buffers3",
         pts = as.matrix(highways_sampled_cast_df),
         n = as.integer(nrow(highways_sampled_cast_df)),
         selected_indices = rep(as.integer(0), nrow(highways_sampled_cast_df)),
         buffer_dist = 1609.34 * 2) # 1 mile * 2 for overlapping radii

# head(oo$selected_indices)
# 
# ggplot() + 
#   geom_sf(data = buffers[1], fill = NA) +
#   geom_sf(data = buffers[2:18], fill = NA, color = 'red') +
#   geom_sf(data = buffers[19], fill = NA, color = 'blue')
# 
# plot(buffers[1])
# plot(buffers[2], col = 'red')

# which(oo$selected_indices != 0) 

# Filter points based on the selected indices
filtered_buffers <- buffers[selected_indices, ]
filtered_points <- highways_sampled_cast[selected_indices, ]

# Plot the state boundaries, population points, and filtered highways
ggplot() +
  geom_sf(data = state_map_sf, fill = "white", color = "black") +
  #geom_sf(data = population_points_sf, color = "forestgreen", size = 0.5) +
  geom_sf(data = highway_lines_sf, color = "orange", size = 0.1, alpha = 0.5) +
  geom_sf(data = filtered_points, color = "red", size = 0.5) +
  theme_minimal() +
  labs(title = "Population Points and Major Highways in Selected Northeastern States")

save(list = ls(), file = 'all_objects.RData')

# load("all_objects.RData")
# baseIDX <- selected_indices
# head(baseIDX)
# head(oo$selected_indices)
