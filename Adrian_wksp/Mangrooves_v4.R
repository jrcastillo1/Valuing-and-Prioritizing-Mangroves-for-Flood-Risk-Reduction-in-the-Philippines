###############################################################################
# Mangroves, Coastal Flood Risk, and Spatial Economic Exposure — Philippines
# Full R script for RStudio
#
# FIXES + LIGHTWEIGHT/CACHE IMPROVEMENTS (minimal but effective):
#  1) Flood (Deltares): Planetary Computer STAC + signing (works) + local caching of RP GeoTIFFs
#  2) Mangroves: avoid 25m classify; downscale to flood template (1km) + cache fraction + distance rasters
#  3) Night lights: fix OpenLandMap STAC reading (relative href + proper JSON Accept header)
#     + cache the processed (cropped+projected+resampled) nightlights raster
#  4) Population: cache the processed (cropped+projected+resampled) population raster
###############################################################################

# ---- Fix PROJ db path (macOS) ----
if (Sys.getenv("PROJ_LIB") == "") {
  hits <- unlist(lapply(.libPaths(), function(p) {
    list.files(p, pattern = "proj\\.db$", recursive = TRUE, full.names = TRUE)
  }))
  if (length(hits) > 0) {
    Sys.setenv(PROJ_LIB = dirname(hits[1]))
    message("Set PROJ_LIB to: ", Sys.getenv("PROJ_LIB"))
  }
}

# ---- GDAL settings (remote reads) ----
Sys.setenv(
  GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR",
  CPL_VSIL_CURL_ALLOWED_EXTENSIONS = ".nc,.tif,.tiff,.json",
  CPL_VSIL_CURL_USE_HEAD = "NO",
  GDAL_HTTP_MULTIRANGE = "YES",
  GDAL_HTTP_MERGE_CONSECUTIVE_RANGES = "YES"
)

# ---- Terra temp folder (optional, helps when disk is tight) ----
dir.create("~/terra_tmp", showWarnings = FALSE)
terraOptions(tempdir = "~/terra_tmp", progress = 1)

# ---- 0) Packages & folders ---------------------------------------------------
pkgs <- c(
  "sf","terra","dplyr","purrr","stringr","exactextractr",
  "jsonlite","httr","readr","ggplot2","tmap","units",
  "WDI","geodata","tidyr","tibble","rlang","curl"
)

to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

`%||%` <- function(x, y) if(is.null(x) || length(x) == 0 || all(is.na(x))) y else x

options(timeout = 600)

dir.create("data", showWarnings = FALSE)
dir.create("data/admin", showWarnings = FALSE, recursive = TRUE)
dir.create("data/flood", showWarnings = FALSE, recursive = TRUE)
dir.create("data/flood/cache", showWarnings = FALSE, recursive = TRUE)
dir.create("data/mangroves", showWarnings = FALSE, recursive = TRUE)
dir.create("data/pop", showWarnings = FALSE, recursive = TRUE)
dir.create("data/ntl", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

# ---- Helper: robust JSON fetch (Accept: application/json) --------------------
safe_get_json <- function(url) {
  tryCatch({
    resp <- httr::GET(
      url,
      httr::add_headers(Accept = "application/json"),
      httr::timeout(120)
    )
    httr::stop_for_status(resp)
    txt <- httr::content(resp, "text", encoding = "UTF-8")
    jsonlite::fromJSON(txt)
  }, error = function(e) NULL)
}

# ---- Helper: make absolute href from a base URL ------------------------------
abs_url <- function(href, base_url) {
  if (is.null(href) || is.na(href) || href == "") return(href)
  if (grepl("^https?://", href)) return(href)
  base_dir <- dirname(base_url)
  href <- sub("^\\./", "", href)
  paste0(base_dir, "/", href)
}

# ---- 1) Study area: Philippines admin + coastal belt -------------------------
crs_m <- "EPSG:3857"

phl_adm2 <- geodata::gadm(country = "PHL", level = 2, path = "data/admin") |>
  st_as_sf() |> st_make_valid()

phl_adm1 <- geodata::gadm(country = "PHL", level = 1, path = "data/admin") |>
  st_as_sf() |> st_make_valid()

phl_poly <- st_union(phl_adm1) |> st_as_sf() |> st_make_valid()

coastline <- st_cast(st_boundary(phl_poly), "MULTILINESTRING")

coastal_km <- 20
coastal_zone <- st_intersection(
  phl_poly,
  st_buffer(st_transform(coastline, 3857), coastal_km * 1000) |> st_transform(st_crs(phl_poly))
) |> st_make_valid()

phl_adm2_m     <- st_transform(phl_adm2, crs_m)
phl_adm1_m     <- st_transform(phl_adm1, crs_m)
coastal_zone_m <- st_transform(coastal_zone, crs_m)

coastal_bbox_wgs <- st_bbox(st_transform(coastal_zone_m, 4326))

# ---- 2) Mangroves: Global Mangrove Watch 2020 (GMW v3.0) ---------------------
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

gmw_vrt <- "data/mangroves/gmw2020.vrt"
if(length(gmw_tifs) > 1) {
  terra::vrt(gmw_tifs, filename = gmw_vrt, overwrite = TRUE)
  gmw <- rast(gmw_vrt)
} else {
  gmw <- rast(gmw_tifs[1])
}

gmw_crop <- crop(
  gmw,
  ext(coastal_bbox_wgs$xmin, coastal_bbox_wgs$xmax, coastal_bbox_wgs$ymin, coastal_bbox_wgs$ymax)
)

gmw_m <- project(gmw_crop, crs_m, method="near")
gmw_m <- mask(gmw_m, vect(coastal_zone_m))
names(gmw_m) <- "mangrove_raw"   # keep raw (no classify at 25m)

# ---- 3) Flood hazard: Deltares (Planetary Computer STAC + sign) + caching ----
pc_stac_search <- "https://planetarycomputer.microsoft.com/api/stac/v1/search"
pc_sign_api    <- "https://planetarycomputer.microsoft.com/api/sas/v1/sign"

sign_pc_href <- function(href) {
  resp <- httr::GET(pc_sign_api, query = list(href = href), httr::timeout(120))
  httr::stop_for_status(resp)
  jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"))$href
}

get_deltares_data_href <- function(dem_name = "NASADEM",
                                   resolution = "1km",
                                   sea_level_year = 2018,
                                   return_period = 2) {
  body <- list(
    collections = list("deltares-floods"),
    limit = 1,
    query = list(
      "deltares:dem_name"       = list(eq = dem_name),
      "deltares:resolution"     = list(eq = resolution),
      "deltares:sea_level_year" = list(eq = sea_level_year),
      "deltares:return_period"  = list(eq = return_period)
    )
  )
  
  resp <- httr::POST(pc_stac_search, body = body, encode = "json", httr::timeout(120))
  httr::stop_for_status(resp)
  
  j <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"), simplifyVector = FALSE)
  
  if (is.null(j$features) || length(j$features) == 0) {
    stop("No Deltares STAC item found for: ", dem_name, ", ", resolution,
         ", year=", sea_level_year, ", rp=", return_period)
  }
  
  feat <- j$features[[1]]
  href <- feat$assets$data$href
  if (is.null(href) || is.na(href) || href == "") stop("Deltares STAC item has no data asset href.")
  href
}

read_deltares_inun <- function(dem_name = "NASADEM",
                               resolution = "1km",
                               sea_level_year = 2018,
                               return_period = 2,
                               var = "inun",
                               cache_dir = "data/flood/cache") {
  
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  
  href        <- get_deltares_data_href(dem_name, resolution, sea_level_year, return_period)
  signed_href <- sign_pc_href(href)
  
  # 1) try streaming
  src1 <- sprintf('NETCDF:"/vsicurl_streaming/%s":%s', signed_href, var)
  r <- try(terra::rast(src1), silent = TRUE)
  if (!inherits(r, "try-error")) return(r)
  
  # 2) fallback: download once
  base <- basename(sub("\\?.*$", "", signed_href))
  dest <- file.path(cache_dir, base)
  
  if (!file.exists(dest)) {
    message("Downloading (one-time): ", base)
    curl::curl_download(signed_href, destfile = dest, quiet = FALSE, mode = "wb")
  }
  
  src2 <- sprintf('NETCDF:"%s":%s', dest, var)
  terra::rast(src2)
}

# Parameters
return_periods <- c(2, 5, 10, 25, 50, 100, 250)
dem_name <- "NASADEM"
res_km   <- "1km"
slr_year <- 2018

flood_list <- list()

for (rp in return_periods) {
  out_tif <- file.path("data/flood", paste0("phl_coastal_inun_rp", rp, "_", res_km, ".tif"))
  
  if (file.exists(out_tif)) {
    message("Loading cached flood RP=", rp, " ...")
    r_m <- rast(out_tif)
  } else {
    message("Fetching flood RP=", rp, " ...")
    r <- read_deltares_inun(
      dem_name = dem_name,
      resolution = res_km,
      sea_level_year = slr_year,
      return_period = rp
    )
    
    r <- crop(r, ext(coastal_bbox_wgs$xmin, coastal_bbox_wgs$xmax,
                     coastal_bbox_wgs$ymin, coastal_bbox_wgs$ymax))
    r_m <- project(r, crs_m, method = "bilinear")
    r_m <- mask(r_m, vect(coastal_zone_m))
    names(r_m) <- paste0("inun_rp", rp)
    
    writeRaster(r_m, out_tif, overwrite = TRUE)
  }
  
  names(r_m) <- paste0("inun_rp", rp)
  flood_list[[as.character(rp)]] <- r_m
}

template <- flood_list[["10"]]

# ---- Downscale mangroves to flood grid (1km) + cache -------------------------
gmw_frac_path <- file.path("data/mangroves", paste0("gmw2020_frac_", res_km, "_phl.tif"))
dist_path     <- file.path("data/mangroves", paste0("dist_to_mangrove_", res_km, "_phl.tif"))

if (file.exists(gmw_frac_path)) {
  message("Loading cached mangrove fraction: ", gmw_frac_path)
  gmw_frac_1km <- rast(gmw_frac_path)
  # ensure aligned (safety)
  if (!compareGeom(gmw_frac_1km, template, stopOnError = FALSE)) {
    gmw_frac_1km <- resample(gmw_frac_1km, template, method="bilinear")
  }
} else {
  message("Creating mangrove fraction at template resolution...")
  gmw_frac_1km <- resample(gmw_m, template, method = "bilinear")
  gmw_frac_1km <- clamp(gmw_frac_1km, 0, 1)
  writeRaster(gmw_frac_1km, gmw_frac_path, overwrite = TRUE)
}
names(gmw_frac_1km) <- "mangrove_frac"
mangrove_presence_1km <- gmw_frac_1km > 0.05

if (file.exists(dist_path)) {
  message("Loading cached distance-to-mangrove: ", dist_path)
  dist_to_mangrove <- rast(dist_path)
  if (!compareGeom(dist_to_mangrove, template, stopOnError = FALSE)) {
    dist_to_mangrove <- resample(dist_to_mangrove, template, method="bilinear")
  }
} else {
  message("Computing distance-to-mangrove (template grid)...")
  dist_to_mangrove <- distance(mangrove_presence_1km)
  writeRaster(dist_to_mangrove, dist_path, overwrite = TRUE)
}
names(dist_to_mangrove) <- "dist_m"

# ---- 4) Population: WorldPop Philippines 2020 + cache processed raster -------
pop_raw_tif  <- "data/pop/phl_ppp_2020.tif"
pop_url      <- "https://data.worldpop.org/GIS/Population/Global_2000_2020/2020/PHL/phl_ppp_2020.tif"
pop_proc_tif <- file.path("data/pop", paste0("phl_ppp_2020_", res_km, "_phl_3857.tif"))

if(!file.exists(pop_raw_tif)) download.file(pop_url, pop_raw_tif, mode = "wb")

if (file.exists(pop_proc_tif)) {
  message("Loading cached processed population: ", pop_proc_tif)
  pop_m <- rast(pop_proc_tif)
} else {
  message("Processing population to template grid...")
  pop <- rast(pop_raw_tif)
  pop <- crop(pop, ext(coastal_bbox_wgs$xmin, coastal_bbox_wgs$xmax,
                       coastal_bbox_wgs$ymin, coastal_bbox_wgs$ymax))
  pop_m <- project(pop, crs_m, method="bilinear")
  pop_m <- mask(pop_m, vect(coastal_zone_m))
  pop_m <- resample(pop_m, template, method="bilinear")
  writeRaster(pop_m, pop_proc_tif, overwrite = TRUE)
}
names(pop_m) <- "pop"

# ---- 5) Night lights: OpenLandMap STAC (fix relative href + JSON header) -----
# NOTE: This keeps your original idea but fixes the common failure:
# - collection/item href can be relative
# - GET must request JSON (Accept: application/json)
catalog_url   <- "http://s3.eu-central-1.wasabisys.com/stac/openlandmap/catalog.json"
collection_id <- "nightlights.average_viirs.v21"
year_pattern  <- "20180101_20181231"
ntl_year      <- 2018

ntl_raw_path  <- file.path("data/ntl", paste0("viirs_annual_", ntl_year, ".tif"))
ntl_proc_path <- file.path("data/ntl", paste0("viirs_annual_", ntl_year, "_", res_km, "_phl_3857.tif"))

if (file.exists(ntl_proc_path)) {
  message("Loading cached processed nightlights: ", ntl_proc_path)
  ntl_m <- rast(ntl_proc_path)
} else {
  # If raw already exists, we can skip STAC and just process it.
  if (!file.exists(ntl_raw_path)) {
    message("Finding OpenLandMap nightlights asset via STAC...")
    cat_json <- safe_get_json(catalog_url)
    if (is.null(cat_json)) stop("Could not read OpenLandMap STAC catalog JSON.")
    
    links <- cat_json$links
    # links often arrive as data.frame
    rels <- links$rel
    hrefs <- links$href
    
    coll_idx <- which(rels == "child" & grepl(collection_id, hrefs))
    if (length(coll_idx) == 0) stop("Could not find nightlights collection link in catalog.")
    coll_url <- abs_url(hrefs[coll_idx[1]], catalog_url)
    
    coll <- safe_get_json(coll_url)
    if (is.null(coll)) stop("Could not read the nightlights collection JSON from STAC (collection URL).")
    
    coll_links <- coll$links
    rels2 <- coll_links$rel
    hrefs2 <- coll_links$href
    
    item_idx <- which(rels2 %in% c("item","child") & grepl(year_pattern, hrefs2))
    if (length(item_idx) == 0) {
      # fallback: take first item link (not ideal, but prevents hard failure)
      item_idx <- which(rels2 %in% c("item","child"))
    }
    if (length(item_idx) == 0) stop("No item links found inside the nightlights collection.")
    
    item_url <- abs_url(hrefs2[item_idx[1]], coll_url)
    item <- safe_get_json(item_url)
    if (is.null(item)) stop("Could not read nightlights item JSON from STAC (item URL).")
    
    # Extract a GeoTIFF asset href
    assets <- item$assets
    asset_hrefs <- purrr::map_chr(assets, ~ .x$href %||% NA_character_)
    tif_href <- asset_hrefs[grepl("\\.tif(f)?$", asset_hrefs, ignore.case = TRUE)][1]
    if (is.na(tif_href) || tif_href == "") {
      # fallback: first asset
      tif_href <- asset_hrefs[1]
    }
    tif_href <- abs_url(tif_href, item_url)
    
    message("Downloading nightlights (one-time): ", basename(tif_href))
    curl::curl_download(tif_href, destfile = ntl_raw_path, quiet = FALSE, mode = "wb")
  } else {
    message("Raw nightlights already exists, processing to template grid...")
  }
  
  # Process raw to template CRS/grid
  ntl <- rast(ntl_raw_path)
  ntl <- crop(ntl, ext(coastal_bbox_wgs$xmin, coastal_bbox_wgs$xmax,
                       coastal_bbox_wgs$ymin, coastal_bbox_wgs$ymax))
  ntl_m <- project(ntl, crs_m, method="bilinear")
  ntl_m <- mask(ntl_m, vect(coastal_zone_m))
  ntl_m <- resample(ntl_m, template, method="bilinear")
  writeRaster(ntl_m, ntl_proc_path, overwrite = TRUE)
}
names(ntl_m) <- "ntl"

# ---- 6) Counterfactual: mangrove attenuation model ---------------------------
# ---- Robust attenuation + protection (handles minor geometry mismatches) ----
attenuation_fraction <- function(dist_m, max_red = 0.35, L = 2000, max_range = 10000) {
  # uses terra::ifel instead of r[dist>...] <- 0 (more robust)
  r <- max_red * exp(-dist_m / L)
  r <- terra::ifel(dist_m > max_range, 0, r)
  terra::clamp(r, 0, 0.95)
}

apply_mangrove_protection <- function(depth_raster, dist_raster,
                                      max_red=0.35, L=2000, max_range=10000) {
  # align distance raster to depth raster if needed
  if (!terra::compareGeom(depth_raster, dist_raster, stopOnError = FALSE)) {
    dist_raster <- terra::resample(dist_raster, depth_raster, method = "bilinear")
  }
  red <- attenuation_fraction(dist_raster, max_red, L, max_range)
  depth_raster * (1 - red)
}

scenario <- list(no_mangroves = flood_list, with_mangroves = list())
for(rp in names(flood_list)) {
  scenario$with_mangroves[[rp]] <- apply_mangrove_protection(
    flood_list[[rp]], dist_to_mangrove, max_red=0.35, L=2000, max_range=10000
  )
  names(scenario$with_mangroves[[rp]]) <- names(flood_list[[rp]])
}

# ---- 7) Exposure: people and night-lights exposed (RP10 by default) ----------
flood_threshold_m <- 0.10

exposure_by_rp <- function(depth_r, pop_r, ntl_r, polygons_sf) {
  flooded <- depth_r > flood_threshold_m
  pop_exp <- exactextractr::exact_extract(pop_r * flooded, polygons_sf, "sum")
  ntl_exp <- exactextractr::exact_extract(ntl_r * flooded, polygons_sf, "sum")
  tibble::tibble(pop_exposed = pop_exp, ntl_exposed = ntl_exp)
}

rp_focus <- "10"

exp_nomang <- exposure_by_rp(scenario$no_mangroves[[rp_focus]], pop_m, ntl_m, phl_adm2_m)
exp_with   <- exposure_by_rp(scenario$with_mangroves[[rp_focus]], pop_m, ntl_m, phl_adm2_m)

exposure_tbl <- phl_adm2_m |>
  st_drop_geometry() |>
  dplyr::transmute(GID_2, NAME_1, NAME_2) |>
  dplyr::bind_cols(
    exp_nomang |> dplyr::rename_with(~paste0(.x,"_nomang")),
    exp_with   |> dplyr::rename_with(~paste0(.x,"_with"))
  ) |>
  dplyr::mutate(
    pop_avoided = pop_exposed_nomang - pop_exposed_with,
    ntl_avoided = ntl_exposed_nomang - ntl_exposed_with
  )

# ---- 8) Monetary valuation: GDP proxy, depth-damage, EAD, avoided EAD --------
gdp_year <- 2018
gdp_usd <- tryCatch({
  WDI::WDI(country="PHL", indicator="NY.GDP.MKTP.CD", start=gdp_year, end=gdp_year) |>
    dplyr::pull(NY.GDP.MKTP.CD) |> as.numeric()
}, error = function(e) NA_real_)
if(is.na(gdp_usd)) gdp_usd <- 376e9

ntl_phl_total <- exactextractr::exact_extract(
  ntl_m,
  (phl_adm1_m |> st_union() |> st_as_sf()),
  "sum"
)[1]

gdp_per_ntl <- gdp_usd / max(ntl_phl_total, 1e-9)
gdp_grid <- ntl_m * gdp_per_ntl
names(gdp_grid) <- "gdp_proxy_usd"

K_Y <- 2.5
asset_value <- gdp_grid * K_Y
names(asset_value) <- "asset_usd"

damage_frac <- function(depth_m) {
  approx(
    x = c(0, 0.1, 0.5, 1, 2, 4, 6),
    y = c(0, 0.01, 0.10, 0.25, 0.45, 0.75, 1.00),
    xout = depth_m,
    rule = 2
  )$y
}

damage_raster <- function(depth_r, asset_r) {
  frac <- app(depth_r, fun = damage_frac)
  frac * asset_r
}

compute_damage_long <- function(depth_maps, polygons_sf, asset_r) {
  n_poly <- nrow(polygons_sf)
  purrr::map_dfr(names(depth_maps), function(rp_chr) {
    rp <- as.numeric(rp_chr)
    dmg_map <- damage_raster(depth_maps[[rp_chr]], asset_r)
    dmg_sum <- exactextractr::exact_extract(dmg_map, polygons_sf, "sum")
    tibble::tibble(poly_id = seq_len(n_poly), rp = rp, damage_usd = as.numeric(dmg_sum))
  })
}

risk_integral_ead <- function(rp, damage) {
  p <- 1 / rp
  ord <- order(p, decreasing = TRUE)
  p <- p[ord]; damage <- damage[ord]
  p <- c(1, p, 0)
  damage <- c(0, damage, tail(damage, 1))
  sum((damage[-length(damage)] + damage[-1]) / 2 * (p[-length(p)] - p[-1]))
}

damage_to_ead <- function(dmg_long) {
  dmg_long |>
    dplyr::group_by(poly_id) |>
    dplyr::summarise(EAD_usd = risk_integral_ead(rp, damage_usd), .groups = "drop")
}

dmg_nomang_long <- compute_damage_long(scenario$no_mangroves,   phl_adm2_m, asset_value)
dmg_with_long   <- compute_damage_long(scenario$with_mangroves, phl_adm2_m, asset_value)

ead_nomang <- damage_to_ead(dmg_nomang_long) |> dplyr::rename(EAD_nomang_usd = EAD_usd)
ead_with   <- damage_to_ead(dmg_with_long)   |> dplyr::rename(EAD_with_usd   = EAD_usd)

valuation_tbl <- phl_adm2_m |>
  st_drop_geometry() |>
  dplyr::transmute(poly_id = dplyr::row_number(), GID_2, NAME_1, NAME_2) |>
  dplyr::left_join(ead_nomang, by = "poly_id") |>
  dplyr::left_join(ead_with,   by = "poly_id") |>
  dplyr::mutate(EAD_avoided_usd = EAD_nomang_usd - EAD_with_usd)

# ---- 9) Mangrove hectares + CBA (NPV, BCR) ----------------------------------
cell_area_ha <- cellSize(template, unit="ha")
mangrove_area_ha_r <- cell_area_ha * gmw_frac_1km
mangrove_ha_by_muni <- exactextractr::exact_extract(mangrove_area_ha_r, phl_adm2_m, "sum")

T_years <- 30
disc <- 0.05
benefit_growth <- 0.00

npv_annuity <- function(annual_benefit, r, T, g = 0) {
  # annual_benefit can be a vector; compute PV factor once, then multiply
  pv_factor <- sum(((1 + g)^(0:(T - 1))) / ((1 + r)^(1:T)))
  annual_benefit * pv_factor
}

cba_tbl <- valuation_tbl |>
  dplyr::mutate(
    mangrove_ha = as.numeric(mangrove_ha_by_muni),
    benefit_per_ha_year = EAD_avoided_usd / pmax(mangrove_ha, 1e-9),
    NPV_benefits_usd = npv_annuity(EAD_avoided_usd, disc, T_years, benefit_growth)
  )

cost_low  <- 1500
cost_mid  <- 4000
cost_high <- 10000

cba_tbl <- cba_tbl |>
  dplyr::mutate(
    NPV_cost_low  = mangrove_ha * cost_low,
    NPV_cost_mid  = mangrove_ha * cost_mid,
    NPV_cost_high = mangrove_ha * cost_high,
    BCR_low  = NPV_benefits_usd / pmax(NPV_cost_low,  1e-9),
    BCR_mid  = NPV_benefits_usd / pmax(NPV_cost_mid,  1e-9),
    BCR_high = NPV_benefits_usd / pmax(NPV_cost_high, 1e-9)
  )

final_tbl <- cba_tbl |>
  dplyr::left_join(exposure_tbl, by = c("GID_2","NAME_1","NAME_2"))

readr::write_csv(final_tbl, "outputs/municipality_mangrove_flood_cba.csv")

# ---- 10) Distributional analysis --------------------------------------------
ntl_by_muni <- exactextractr::exact_extract(ntl_m, phl_adm2_m, "sum")
pop_by_muni <- exactextractr::exact_extract(pop_m, phl_adm2_m, "sum")

dist_tbl <- final_tbl |>
  dplyr::mutate(
    ntl_sum = as.numeric(ntl_by_muni),
    pop_sum = as.numeric(pop_by_muni),
    ntl_pc  = ntl_sum / pmax(pop_sum, 1e-9),
    benefit_pc = EAD_avoided_usd / pmax(pop_sum, 1e-9)
  ) |>
  dplyr::arrange(ntl_pc) |>
  dplyr::mutate(
    pop_share = pop_sum / sum(pop_sum, na.rm=TRUE),
    cum_pop   = cumsum(pop_share),
    benefit_share = EAD_avoided_usd / sum(EAD_avoided_usd, na.rm=TRUE),
    cum_benefit   = cumsum(replace_na(benefit_share, 0))
  )

bottom40_benefit_share <- dist_tbl |>
  dplyr::filter(cum_pop <= 0.40) |>
  dplyr::summarise(share = sum(benefit_share, na.rm=TRUE)) |>
  dplyr::pull(share)

message(sprintf("Share of avoided EAD accruing to bottom 40%% (ntl/capita proxy): %.3f", bottom40_benefit_share))

# ---- 11) Maps ----------------------------------------------------------------
tmap_mode("plot")

phl_map <- phl_adm2_m |>
  dplyr::left_join(final_tbl, by = c("GID_2","NAME_1","NAME_2"))

p1 <- tm_shape(gmw_frac_1km) + tm_raster(title=paste0("Mangrove fraction (", res_km, ")")) +
  tm_shape(coastal_zone_m) + tm_borders(lwd=1) +
  tm_layout(main.title="Mangrove cover (fraction) in coastal belt")

p2 <- tm_shape(scenario$no_mangroves[[rp_focus]]) + tm_raster(title=paste0("Flood depth (m), RP", rp_focus)) +
  tm_shape(coastal_zone_m) + tm_borders(lwd=1) +
  tm_layout(main.title=paste0("Coastal inundation depth (RP", rp_focus, ", baseline)"))

p3 <- tm_shape(phl_map) + tm_polygons("EAD_avoided_usd", style="quantile", n=5,
                                      title="Avoided EAD (USD/yr)") +
  tm_layout(main.title="Annual expected damages avoided by mangrove protection")

# Use base png device (no Cairo)
# Guaranteed fallback (base graphics device)
grDevices::png("outputs/map_mangroves.png", width=1800, height=1200, res=300)
print(p1)
dev.off()

grDevices::png("outputs/map_flood_rp10.png", width=1800, height=1200, res=300)
print(p2)
dev.off()

grDevices::png("outputs/map_avoided_ead.png", width=1800, height=1200, res=300)
print(p3)
dev.off()

# ---- 12) Sensitivity analysis ------------------------------------------------
grid <- tidyr::expand_grid(
  max_red = c(0.20, 0.35, 0.50),
  L       = c(1000, 2000, 4000)
)

run_one <- function(max_red, L) {
  with_maps <- lapply(flood_list, \(r) apply_mangrove_protection(r, dist_to_mangrove,
                                                                 max_red=max_red, L=L, max_range=10000))
  dmg_n <- compute_damage_long(flood_list, phl_adm2_m, asset_value)
  dmg_w <- compute_damage_long(with_maps,  phl_adm2_m, asset_value)
  
  ead_n <- damage_to_ead(dmg_n)
  ead_w <- damage_to_ead(dmg_w)
  
  sum(ead_n$EAD_usd - ead_w$EAD_usd, na.rm=TRUE)
}

sens <- grid |>
  dplyr::mutate(avoided_EAD_total = purrr::pmap_dbl(list(max_red, L), run_one))

readr::write_csv(sens, "outputs/sensitivity_avoided_ead_total.csv")

message("Done. Key file: outputs/municipality_mangrove_flood_cba.csv")
###############################################################################