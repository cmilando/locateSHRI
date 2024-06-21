library(tidyverse)
load("all_objects.RData")

# population

ggplot() +
  geom_sf(data = state_map_sf, fill = "white", color = "black") +
  geom_sf(data = population_points_sf, color = "forestgreen", size = 0.5) +
  #geom_sf(data = highway_lines_sf, color = "orange", size = 0.1, alpha = 0.5) +
  #geom_sf(data = filtered_points, color = "red", size = 0.5) +
  theme_minimal() +
  labs(title = "Population Points and Major Highways in Selected Northeastern States")

# Points and highways

ggplot() +
  geom_sf(data = state_map_sf, fill = "white", color = "black") +
  geom_sf(data = population_points_sf, color = "forestgreen", size = 0.5) +
  geom_sf(data = highway_lines_sf, color = "orange", size = 0.1, alpha = 0.5) +
  #geom_sf(data = filtered_points, color = "red", size = 0.5) +
  theme_minimal() +
  labs(title = "Population Points and Major Highways in Selected Northeastern States")


# Potential sites

ggplot() +
  geom_sf(data = state_map_sf, fill = "white", color = "black") +
  #geom_sf(data = population_points_sf, color = "forestgreen", size = 0.5) +
  geom_sf(data = highway_lines_sf, color = "orange", size = 0.1, alpha = 0.5) +
  geom_sf(data = filtered_points, color = "red", size = 0.5) +
  theme_minimal() +
  labs(title = "Population Points and Major Highways in Selected Northeastern States")
