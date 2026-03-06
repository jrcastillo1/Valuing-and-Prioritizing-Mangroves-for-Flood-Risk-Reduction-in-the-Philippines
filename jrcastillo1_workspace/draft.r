# ---- 0. Setup ----
pkgs <- c(
  "sf","terra","dplyr","purrr","stringr","exactextractr",
  "jsonlite","httr","readr","ggplot2","tmap","units",
  "WDI"  # optional: pull GDP from World Bank API
)

to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(to_install)) install.packages(to_install)

invisible(lapply(pkgs, library, character.only = TRUE))

options(timeout = 600)        # allow big downloads
terraOptions(progress = 1)

dir.create("data", showWarnings = FALSE)
dir.create("data/admin", showWarnings = FALSE, recursive = TRUE)
dir.create("data/flood", showWarnings = FALSE, recursive = TRUE)
dir.create("data/mangroves", showWarnings = FALSE, recursive = TRUE)
dir.create("data/pop", showWarnings = FALSE, recursive = TRUE)
dir.create("data/ntl", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)


# ---- 1. Study area: Philippines + coastal belt ----
# Admin boundaries from GADM via geodata package (downloaded automatically)
# If geodata isn't installed, install it:
if(!"geodata" %in% installed.packages()[,"Package"]) install.packages("geodata")
library(geodata)

# Barangay
phl_adm3 <- geodata::gadm(country = "PHL", level = 3, path = "data/admin") |>
  st_as_sf() |>
  st_make_valid()

ggplot(phl_adm3) + geom_sf(fill = "honeydew", color = "darkgreen") + labs(title = "National Boundary Union")

# Municipality/City
phl_adm2 <- geodata::gadm(country = "PHL", level = 2, path = "data/admin") |>
  st_as_sf() |>
  st_make_valid()

ggplot(phl_adm2) + geom_sf(fill = "honeydew", color = "darkgreen") + labs(title = "National Boundary Union")

# Region
phl_adm1 <- geodata::gadm(country = "PHL", level = 1, path = "data/admin") |>
  st_as_sf() |>
  st_make_valid()

ggplot(phl_adm1) + geom_sf(fill = "honeydew", color = "darkgreen") + labs(title = "National Boundary Union")

# Country
phl_adm0 <- geodata::gadm(country = "PHL", level = 0, path = "data/admin") |>
  st_as_sf() |>
  st_make_valid()

ggplot(phl_adm0) + geom_sf(fill = "honeydew", color = "darkgreen") + labs(title = "National Boundary Union")


phl_poly <- st_union(phl_adm0) |> st_as_sf() |> st_make_valid()

ggplot(phl_poly) + geom_sf()

# Coastal belt: buffer the national boundary line and intersect with land polygon
coastline <- st_cast(st_boundary(phl_poly), "MULTILINESTRING")

ggplot(coastline) + geom_sf()

coastal_km <- 20
coastal_zone <- st_intersection(
  phl_poly,
  st_buffer(st_transform(coastline, 3857), coastal_km * 1000)
  |> st_transform(st_crs(phl_poly))
) |> st_make_valid()

# Use a metric CRS for distance computations & raster alignment
crs_m <- "EPSG:3857"
phl_adm2_m   <- st_transform(phl_adm2, crs_m)
coastal_zone_m <- st_transform(coastal_zone, crs_m)


# NOTE: GO TO BARANGAY LEVEL, BUT FILTER INTERECTION OF COASTAL BOUNDARY WITH REGIONS THEN FILTER TO MUNICITIES THAT INTERSECT COASTAL BOUNDARY THEN FILTER TO BARANGAYS THAT INTERSECT COASTAL BOUNDARY PLUS INTERSECT ONE DEGREE NEIGHBORING BARANGAYS





# ---- 2. Mangroves (GMW v3.0, year 2020) ----
gmw_zip <- "data/mangroves/gmw_v3_2020_gtiff.zip"
gmw_url <- "https://zenodo.org/records/6894273/files/gmw_v3_2020_gtiff.zip?download=1"

if(!file.exists(gmw_zip)) download.file(gmw_url, gmw_zip, mode = "wb")

gmw_dir <- "data/mangroves/gmw2020"
if(!dir.exists(gmw_dir)) dir.create(gmw_dir, recursive = TRUE)

if(length(list.files(gmw_dir, pattern="\\.tif$", recursive=TRUE)) == 0) {
  unzip(gmw_zip, exdir = gmw_dir)
}

gmw_tifs <- list.files(gmw_dir, pattern="\\.tif$", full.names=TRUE, recursive=TRUE)
stopifnot(length(gmw_tifs) >= 1)

# If multiple tiles exist, build a VRT and read as a single raster
gmw_vrt <- "data/mangroves/gmw2020.vrt"
if(length(gmw_tifs) > 1) {
  terra::vrt(gmw_tifs, filename = gmw_vrt, overwrite = TRUE)
  gmw <- rast(gmw_vrt)
} else {
  gmw <- rast(gmw_tifs[1])
}


# Crop to coastal zone bbox (in WGS84), then project to metric CRS
coastal_bbox_wgs <- st_bbox(st_transform(coastal_zone_m, 4326))
gmw_crop <- crop(gmw, ext(coastal_bbox_wgs$xmin, coastal_bbox_wgs$xmax,
                          coastal_bbox_wgs$ymin, coastal_bbox_wgs$ymax))

gmw_m <- project(gmw_crop, crs_m, method="near")
gmw_m <- mask(gmw_m, vect(coastal_zone_m))

# Ensure binary 0/1
gmw_m <- classify(gmw_m, rcl = matrix(c(-Inf,0.5,0, 0.5,Inf,1), ncol=3, byrow=TRUE))
names(gmw_m) <- "mangrove"



###### VIZ
# gmw_df <- as.data.frame(gmw_m, xy = TRUE, na.rm = TRUE)

# ggplot() +
#     # Draw the coastal zone boundary for context
#     geom_sf(data = coastal_zone_m, fill = "gray95", color = "gray80") +
#     # Draw the mangrove pixels
#     # Since it's binary, we use geom_tile or geom_raster
#     geom_raster(data = gmw_plot_df, aes(x = x, y = y, fill = factor(mangrove))) +
#     # Style the colors: 0 (No mangrove) and 1 (Mangrove)
#     scale_fill_manual(
#     values = c("1" = "darkgreen", "0" = "transparent"),
#     labels = c("1" = "Mangrove Presence"),
#     name = "Legend",
#     na.translate = FALSE
#     ) +
#     # Use coord_sf to ensure the metric CRS (3857) is handled correctly
#     coord_sf(crs = st_crs(3857)) +
#     theme_minimal() +
#     theme(aspect.ratio = 1) +
#     labs(
#     title = "Mangrove Distribution (GMW 2020)",
#     subtitle = "Clipped to 20km Coastal Zone (Philippines)",
#     x = "Easting (m)",
#     y = "Northing (m)"
#     )


# 1. Aggregate the raster: factor = 10 means 10x10 pixels become 1 larger pixel
# fun = "max" ensures if there was ANY mangrove in that 10x10 area, the new pixel is 1.
gmw_low_res <- aggregate(gmw_m, fact = 10, fun = "max")
# 2. Now convert the much smaller raster to a data frame
gmw_df <- as.data.frame(gmw_low_res, xy = TRUE, na.rm = TRUE)
# 3. Check the size - it should be much smaller now
nrow(gmw_df)

ggplot() +
    geom_sf(data = phl_adm2_m, fill = "gray95", color = "gray80", size = 0.1) +
    geom_point(data = gmw_df, aes(x = x, y = y),  color = "darkgreen", shape = 15, size = 1.5) +    coord_sf(crs = st_crs(3857)) +
    theme_minimal() +
    labs(title = "Mangrove Presence (Downsampled 10x)",
    subtitle = "Showing the physical distribution along the coast")


library(exactextractr)
library(viridis)

# 1. Calculate the sum of mangrove pixels for each municipality
# Since pixels are 1 (mangrove) or 0 (not), the sum is the count of mangrove pixels.
phl_adm2_m$mangrove_pixels <- exact_extract(gmw_m, phl_adm2_m, 'sum')

# 2. Convert pixel count to Area (Square Kilometers)
# We calculate: (Pixel Count * Pixel Area) / 1,000,000
# Get the resolution (width/height) of a pixel in meters
res_val <- res(gmw_m)[1] 
phl_adm2_m$mangrove_km2 <- (phl_adm2_m$mangrove_pixels * (res_val^2)) / 1000000

# 3. Create the Choropleth Map
ggplot(phl_adm2_m) +
  geom_sf(aes(fill = mangrove_km2), color = "gray90", linewidth = 0.05) +
  # Custom gradient from white to dark green
  scale_fill_gradient(
    low = "white", 
    high = "darkgreen", 
    name = expression(km^2),
    na.value = "white" # Keeps non-mangrove areas clean
  ) +
  theme_minimal() +
  theme(
    aspect.ratio = 1, 
    legend.position = "right",
    # Adding a light blue background for the 'sea' makes the map pop
    panel.background = element_rect(fill = "aliceblue", color = NA) 
  ) +
  labs(
    title = "Mangrove Extent by Municipality",
    subtitle = "Total area of mangroves within the 20km coastal belt",
    caption = "Data: Global Mangrove Watch 2020"
  )
