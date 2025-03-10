library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)  # For basemaps
library(rosm)       # For OpenStreetMap tiles
library(ggmap)

# Create the map
register_google(key = "AIzaSyChTsYLBX8qJeNRdxhwgSruZfu4eXFKrxI")


# Read data and convert WKT to sf
sites <- read.csv("samplelocations.csv")
sites_sf <- st_as_sf(sites, wkt = "WKT", crs = 4326)  # Convert WKT to sf object

# Register Google API Key (if not already registered)
# register_google(key = "YOUR_GOOGLE_MAPS_API_KEY")

# Ensure sites_sf is an sf object
sites_sf <- st_as_sf(sites_sf, coords = c("Longitude", "Latitude"), crs = 4326)


# Define a custom color palette with HEX codes
custom_colors <- c(
  "NOCA" = "#6247aa",
  "PYRR" = "#2d6a4f"
)

# Define a buffer to zoom out slightly
buffer_factor <- 0.01  # 1% increase in range

# Compute bounding box with buffer
min_long <- min(sites_sf$Longitude, na.rm = TRUE)
max_long <- max(sites_sf$Longitude, na.rm = TRUE)
min_lat <- min(sites_sf$Latitude, na.rm = TRUE)
max_lat <- max(sites_sf$Latitude, na.rm = TRUE)

long_range <- max_long - min_long
lat_range <- max_lat - min_lat

bbox <- c(
  left = min_long - buffer_factor * long_range,
  bottom = min_lat - buffer_factor * lat_range,
  right = max_long + buffer_factor * long_range,
  top = max_lat + buffer_factor * lat_range
)

# Get Google Terrain Map
terrain_map <- get_googlemap(
  center = c(lon = mean(c(bbox["left"], bbox["right"])), lat = mean(c(bbox["bottom"], bbox["top"]))),
  zoom = 7,  # Adjust manually if needed
  maptype = "hybrid",
  scale = 2
)



# Plot with custom colors
ggmap(terrain_map) +
  geom_sf(data = sites_sf, aes(color = Species), inherit.aes = FALSE, size = 3, alpha = 0.8) +
  scale_color_manual(values = custom_colors) +  # Use pre-defined colors
  theme_minimal() +
  labs(color = "Species & Population")  # Legend label


ggsave(
  "allplot.pdf",
  plot = last_plot(),
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)


# plot just urban 

# Ensure you load dplyr and sf packages
library(dplyr)

urban_sites_sf <- sites_sf %>% filter(Population == "Urban")

# Define a buffer to zoom out slightly
buffer_factor <- 0.01  # 1% increase in range

# Compute bounding box with buffer
min_long <- min(urban_sites_sf$Longitude, na.rm = TRUE)
max_long <- max(urban_sites_sf$Longitude, na.rm = TRUE)
min_lat <- min(urban_sites_sf$Latitude, na.rm = TRUE)
max_lat <- max(urban_sites_sf$Latitude, na.rm = TRUE)

long_range <- max_long - min_long
lat_range <- max_lat - min_lat

bbox <- c(
  left = min_long - buffer_factor * long_range,
  bottom = min_lat - buffer_factor * lat_range,
  right = max_long + buffer_factor * long_range,
  top = max_lat + buffer_factor * lat_range
)

# Get Google Terrain Map
terrain_map <- get_googlemap(
  center = c(lon = mean(c(bbox["left"], bbox["right"])), lat = mean(c(bbox["bottom"], bbox["top"]))),
  zoom = 10,  # Adjust manually if needed
  maptype = "hybrid",
  scale = 2
)

# Plot with Google Terrain Map and jittered points
# Subset sites_sf to include only Urban population
urban_sites_sf <- sites_sf %>% filter(Population == "Urban")

# Plot with Google Terrain Map and jittered points
ggmap(terrain_map) +
  geom_jitter(data = urban_sites_sf, aes(x = Longitude, y = Latitude, color = Species), 
              inherit.aes = FALSE, size = 3, alpha = 0.8, width = 0.01, height = 0.01) +  # Adjust jitter width and height
  scale_color_manual(values = custom_colors) +  # Use pre-defined colors
  theme_minimal() +
  labs(color = "Species & Population (Urban Only)")  # Updated legend label


ggsave(
  "urbanplot.pdf",
  plot = last_plot(),
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)
