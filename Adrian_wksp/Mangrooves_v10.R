###############################################################################
# Mangrove Protection, Coastal Flood Risk, and Economic Exposure — Philippines
###############################################################################

# =========================
# 0) SETTINGS, WORKING DIRECTORY
# =========================
RUN_SENS_ALPHA_L <- TRUE   # Table 3 (alpha-L grid)
RUN_SENS_DISC    <- TRUE   # discount-rate grid for NPV/BCR
RP_FOCUS         <- 10      # for exposure + flood map
COASTAL_KM       <- 20      # coastal belt width in km (PDF: 20 km)

if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  } else if (!is.null(sys.frames()[[1]]$ofile)) {
    script_dir <- dirname(normalizePath(sys.frames()[[1]]$ofile))
  } else {
    script_dir <- getwd()
  }
}

setwd(script_dir)
cat("Working directory is:", getwd(), "\n")


# =========================
# 0.1) PACKAGES & SYSTEM
# =========================
pkgs <- c(
  "sf","terra","dplyr","purrr","stringr","exactextractr",
  "jsonlite","httr","readr","tmap","WDI","geodata","tidyr","tibble","rlang","curl", "ggplot2"
)
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

`%||%` <- function(x, y) if(is.null(x) || length(x)==0 || all(is.na(x))) y else x
options(timeout = 600)

# Fix PROJ on macOS (if needed)
if (Sys.getenv("PROJ_LIB") == "") {
  hits <- unlist(lapply(.libPaths(), function(p) {
    list.files(p, pattern = "proj\\.db$", recursive = TRUE, full.names = TRUE)
  }))
  if (length(hits) > 0) Sys.setenv(PROJ_LIB = dirname(hits[1]))
}

# GDAL settings for remote reads
Sys.setenv(
  GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR",
  CPL_VSIL_CURL_ALLOWED_EXTENSIONS = ".nc,.tif,.tiff,.json",
  CPL_VSIL_CURL_USE_HEAD = "NO",
  GDAL_HTTP_MULTIRANGE = "YES",
  GDAL_HTTP_MERGE_CONSECUTIVE_RANGES = "YES"
)

# Terra temp (stable runs)
dir.create("~/terra_tmp", showWarnings = FALSE)
terra::terraOptions(tempdir="~/terra_tmp", progress=1, todisk=TRUE, memfrac=0.6)

# Folders
dir.create("data", showWarnings = FALSE)
dir.create("data/admin", showWarnings = FALSE, recursive = TRUE)
dir.create("data/flood", showWarnings = FALSE, recursive = TRUE)
dir.create("data/flood/cache", showWarnings = FALSE, recursive = TRUE)
dir.create("data/mangroves", showWarnings = FALSE, recursive = TRUE)
dir.create("data/pop", showWarnings = FALSE, recursive = TRUE)
dir.create("data/ntl", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

# =========================
# 0.2) PAPER PARAMETERS (PDF)
# =========================
# Section 4.2–4.6 parameters (central case + sensitivity grids) :contentReference[oaicite:1]{index=1}
tau     <- 0.05     # Eq (1)
h0      <- 0.10     # Eq (6)-(11)
alpha0  <- 0.25     # Eq (3) central
L0      <- 1000     # Eq (3) central (m)
dmax    <- 10000    # Eq (3) cutoff (m)
K_Y     <- 2.0      # Eq (13)
T_years <- 30       # Eq (19)
disc    <- 0.05     # Eq (19)
g       <- 0.00     # Eq (19)

# Table 3 grid (PDF)
alpha_grid <- c(0.15, 0.25, 0.40)
L_grid     <- c(500, 1000, 2000)
disc_grid  <- c(0.03, 0.05, 0.08)

# Flood hazard parameters (PDF: NASADEM 1km, 2018 baseline, RP set)
return_periods <- c(2, 5, 10, 25, 50, 100, 250)
dem_name <- "NASADEM"
res_km   <- "1km"
slr_year <- 2018
crs_m    <- "EPSG:3857"

# =========================
# 0.3) HELPER FUNCTIONS
# =========================

# ---- robust JSON fetch (STAC)
safe_get_json <- function(url) {
  tryCatch({
    resp <- httr::GET(url, httr::add_headers(Accept="application/json"), httr::timeout(120))
    httr::stop_for_status(resp)
    jsonlite::fromJSON(httr::content(resp, "text", encoding="UTF-8"))
  }, error = function(e) NULL)
}

# ---- convert relative STAC href to absolute
abs_url <- function(href, base_url) {
  if (is.null(href) || is.na(href) || href == "") return(href)
  if (grepl("^https?://", href)) return(href)
  base_dir <- dirname(base_url)
  href <- sub("^\\./", "", href)
  paste0(base_dir, "/", href)
}

# ---- cache helper for rasters
read_or_build_raster <- function(path, build_fun) {
  if (file.exists(path)) return(terra::rast(path))
  r <- build_fun()
  terra::writeRaster(r, path, overwrite=TRUE)
  terra::rast(path)
}

# =========================
# 1) DATA: STUDY AREA (PDF 3.1)
# =========================
# Build coastal belt and admin boundaries (municipality level)

phl_adm2 <- geodata::gadm(country="PHL", level=2, path="data/admin") |> st_as_sf() |> st_make_valid()
phl_adm1 <- geodata::gadm(country="PHL", level=1, path="data/admin") |> st_as_sf() |> st_make_valid()
phl_poly <- st_union(phl_adm1) |> st_as_sf() |> st_make_valid()

coastline <- st_cast(st_boundary(phl_poly), "MULTILINESTRING")

coastal_zone <- st_intersection(
  phl_poly,
  st_buffer(st_transform(coastline, 3857), COASTAL_KM * 1000) |> st_transform(st_crs(phl_poly))
) |> st_make_valid()

phl_adm2_m     <- st_transform(phl_adm2, crs_m)
phl_adm1_m     <- st_transform(phl_adm1, crs_m)
coastal_zone_m <- st_transform(coastal_zone, crs_m)
coastal_bbox_wgs <- st_bbox(st_transform(coastal_zone_m, 4326))

# =========================
# 2) DATA: MANGROVES (PDF 3.2)
# =========================
gmw_zip <- "data/mangroves/gmw_v3_2020_gtiff.zip"
gmw_url <- "https://zenodo.org/records/6894273/files/gmw_v3_2020_gtiff.zip?download=1"
if(!file.exists(gmw_zip)) download.file(gmw_url, gmw_zip, mode="wb")

gmw_dir <- "data/mangroves/gmw2020"
dir.create(gmw_dir, showWarnings=FALSE, recursive=TRUE)
if(length(list.files(gmw_dir, pattern="\\.tif$", recursive=TRUE)) == 0) unzip(gmw_zip, exdir=gmw_dir)

gmw_tifs <- list.files(gmw_dir, pattern="\\.tif$", full.names=TRUE, recursive=TRUE)
stopifnot(length(gmw_tifs) >= 1)

gmw_vrt <- "data/mangroves/gmw2020.vrt"
if(length(gmw_tifs) > 1) {
  terra::vrt(gmw_tifs, filename = gmw_vrt, overwrite = TRUE)
  gmw <- terra::rast(gmw_vrt)
} else {
  gmw <- terra::rast(gmw_tifs[1])
}

gmw_m <- gmw |>
  terra::crop(terra::ext(coastal_bbox_wgs$xmin, coastal_bbox_wgs$xmax,
                         coastal_bbox_wgs$ymin, coastal_bbox_wgs$ymax)) |>
  terra::project(crs_m, method="near") |>
  terra::mask(terra::vect(coastal_zone_m))
names(gmw_m) <- "mangrove_raw"

# =========================
# 3) DATA: FLOOD HAZARD (PDF 3.3)
# =========================
pc_stac_search <- "https://planetarycomputer.microsoft.com/api/stac/v1/search"
pc_sign_api    <- "https://planetarycomputer.microsoft.com/api/sas/v1/sign"

sign_pc_href <- function(href) {
  resp <- httr::GET(pc_sign_api, query=list(href=href), httr::timeout(120))
  httr::stop_for_status(resp)
  jsonlite::fromJSON(httr::content(resp, "text", encoding="UTF-8"))$href
}

get_deltares_data_href <- function(dem_name, resolution, sea_level_year, return_period) {
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
  resp <- httr::POST(pc_stac_search, body=body, encode="json", httr::timeout(120))
  httr::stop_for_status(resp)
  j <- jsonlite::fromJSON(httr::content(resp, "text", encoding="UTF-8"), simplifyVector=FALSE)
  feat <- j$features[[1]]
  feat$assets$data$href
}

read_deltares_inun <- function(dem_name, resolution, sea_level_year, return_period,
                               var="inun", cache_dir="data/flood/cache") {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  href <- get_deltares_data_href(dem_name, resolution, sea_level_year, return_period)
  signed <- sign_pc_href(href)
  
  # Try streaming first
  src <- sprintf('NETCDF:"/vsicurl_streaming/%s":%s', signed, var)
  r <- try(terra::rast(src), silent=TRUE)
  if (!inherits(r, "try-error")) return(r)
  
  # Fallback: download once
  base <- basename(sub("\\?.*$", "", signed))
  dest <- file.path(cache_dir, base)
  if(!file.exists(dest)) curl::curl_download(signed, destfile=dest, quiet=FALSE, mode="wb")
  terra::rast(sprintf('NETCDF:"%s":%s', dest, var))
}

# Build flood_list with caching (one tif per RP)
flood_list <- list()
for (rp in return_periods) {
  out_tif <- file.path("data/flood", paste0("phl_coastal_inun_rp", rp, "_", res_km, ".tif"))
  
  r_m <- if (file.exists(out_tif)) {
    message("Loading cached flood RP=", rp)
    terra::rast(out_tif)
  } else {
    message("Fetching flood RP=", rp)
    r <- read_deltares_inun(dem_name, res_km, slr_year, rp)
    
    r_m <- r |>
      terra::crop(terra::ext(coastal_bbox_wgs$xmin, coastal_bbox_wgs$xmax,
                             coastal_bbox_wgs$ymin, coastal_bbox_wgs$ymax)) |>
      terra::project(crs_m, method="bilinear") |>
      terra::mask(terra::vect(coastal_zone_m))
    
    names(r_m) <- paste0("inun_rp", rp)
    terra::writeRaster(r_m, out_tif, overwrite=TRUE)
    terra::rast(out_tif)
  }
  
  names(r_m) <- paste0("inun_rp", rp)
  flood_list[[as.character(rp)]] <- r_m
}

template <- flood_list[[as.character(RP_FOCUS)]]

# =========================
# 4) DATA: POPULATION (PDF 3.4) + CACHING
# =========================
pop_raw <- "data/pop/phl_ppp_2020.tif"
pop_url <- "https://data.worldpop.org/GIS/Population/Global_2000_2020/2020/PHL/phl_ppp_2020.tif"
if(!file.exists(pop_raw)) download.file(pop_url, pop_raw, mode="wb")

pop_proc <- file.path("data/pop", paste0("phl_ppp_2020_", res_km, "_phl_3857.tif"))
pop_m <- read_or_build_raster(pop_proc, function() {
  terra::rast(pop_raw) |>
    terra::crop(terra::ext(coastal_bbox_wgs$xmin, coastal_bbox_wgs$xmax,
                           coastal_bbox_wgs$ymin, coastal_bbox_wgs$ymax)) |>
    terra::project(crs_m, method="bilinear") |>
    terra::mask(terra::vect(coastal_zone_m)) |>
    terra::resample(template, method="bilinear")
})
names(pop_m) <- "pop"

# =========================
# 5) DATA: NIGHT LIGHTS (PDF 3.5) + CACHING
# =========================
ntl_year <- 2018
ntl_raw  <- file.path("data/ntl", paste0("viirs_annual_", ntl_year, ".tif"))
ntl_proc <- file.path("data/ntl", paste0("viirs_annual_", ntl_year, "_", res_km, "_phl_3857.tif"))

if (!file.exists(ntl_proc)) {
  if (!file.exists(ntl_raw)) {
    catalog_url   <- "http://s3.eu-central-1.wasabisys.com/stac/openlandmap/catalog.json"
    collection_id <- "nightlights.average_viirs.v21"
    year_pattern  <- "20180101_20181231"
    
    message("Finding OpenLandMap nightlights asset via STAC...")
    cat_json <- safe_get_json(catalog_url); stopifnot(!is.null(cat_json))
    
    links <- cat_json$links
    coll_idx <- which(links$rel == "child" & grepl(collection_id, links$href))
    stopifnot(length(coll_idx) > 0)
    coll_url <- abs_url(links$href[coll_idx[1]], catalog_url)
    
    coll <- safe_get_json(coll_url); stopifnot(!is.null(coll))
    coll_links <- coll$links
    item_idx <- which(coll_links$rel %in% c("item","child") & grepl(year_pattern, coll_links$href))
    if (length(item_idx) == 0) item_idx <- which(coll_links$rel %in% c("item","child"))
    stopifnot(length(item_idx) > 0)
    
    item_url <- abs_url(coll_links$href[item_idx[1]], coll_url)
    item <- safe_get_json(item_url); stopifnot(!is.null(item))
    
    assets <- item$assets
    asset_hrefs <- purrr::map_chr(assets, ~ .x$href %||% NA_character_)
    tif_href <- asset_hrefs[grepl("\\.tif(f)?$", asset_hrefs, ignore.case = TRUE)][1]
    if (is.na(tif_href) || tif_href == "") tif_href <- asset_hrefs[1]
    tif_href <- abs_url(tif_href, item_url)
    
    message("Downloading nightlights (one-time): ", basename(tif_href))
    curl::curl_download(tif_href, destfile=ntl_raw, quiet=FALSE, mode="wb")
  }
  
  # Process raw -> coastal belt -> template grid
  ntl_m <- terra::rast(ntl_raw) |>
    terra::crop(terra::ext(coastal_bbox_wgs$xmin, coastal_bbox_wgs$xmax,
                           coastal_bbox_wgs$ymin, coastal_bbox_wgs$ymax)) |>
    terra::project(crs_m, method="bilinear") |>
    terra::mask(terra::vect(coastal_zone_m)) |>
    terra::resample(template, method="bilinear")
  
  terra::writeRaster(ntl_m, ntl_proc, overwrite=TRUE)
}

ntl_m <- terra::rast(ntl_proc)
names(ntl_m) <- "ntl"

# =========================
# 6) METHODOLOGY: MANGROVE COUNTERFACTUAL (PDF 4.2, Eq 1–5)
# =========================

# Eq (1) MangrovePresence(x) = 1(MangroveFraction(x) > tau)
gmw_frac_path <- file.path("data/mangroves", paste0("gmw2020_frac_", res_km, "_phl.tif"))
gmw_frac_1km <- read_or_build_raster(gmw_frac_path, function() {
  terra::resample(gmw_m, template, method="bilinear") |> terra::clamp(0, 1)
})
names(gmw_frac_1km) <- "mangrove_frac"
mangrove_presence <- gmw_frac_1km > tau

# Eq (2) d(x) = distance to nearest mangrove cell (needs NA-based source!)
dist_path <- file.path("data/mangroves", paste0("dist_to_mangrove_", res_km, "_phl.tif"))
dist_to_mangrove <- read_or_build_raster(dist_path, function() {
  src <- terra::ifel(mangrove_presence, 1, NA)   # mangrove cells are non-NA
  terra::distance(src)
})
# Sanity check: if cached file was created incorrectly earlier, fix automatically
maxd <- terra::global(dist_to_mangrove, "max", na.rm=TRUE)[1,1]
if (is.na(maxd) || maxd < 1) {
  message("dist_to_mangrove cache looks invalid (max<1). Recomputing...")
  src <- terra::ifel(mangrove_presence, 1, NA)
  dist_to_mangrove <- terra::distance(src)
  terra::writeRaster(dist_to_mangrove, dist_path, overwrite=TRUE)
  dist_to_mangrove <- terra::rast(dist_path)
}
names(dist_to_mangrove) <- "dist_m"

# Eq (3) attenuation kernel r(x)
attenuation_r <- function(dist_m, alpha=alpha0, L=L0, dmax_=dmax) {
  r <- alpha * exp(-dist_m / L)
  r <- terra::ifel(dist_m > dmax_, 0, r)
  terra::clamp(r, 0, 0.95)
}

# Eq (4) Depth_with = Depth_no * (1 - r(x))
apply_protection <- function(depth_r, dist_r, alpha=alpha0, L=L0, dmax_=dmax) {
  if (!terra::compareGeom(depth_r, dist_r, stopOnError = FALSE)) {
    dist_r <- terra::resample(dist_r, depth_r, method="bilinear")
  }
  red <- attenuation_r(dist_r, alpha, L, dmax_)
  depth_r * (1 - red)
}

# Build scenario maps (no vs with) for all return periods
scenario <- list(no_mangroves = flood_list, with_mangroves = list())
for (rp in names(flood_list)) {
  scenario$with_mangroves[[rp]] <- apply_protection(flood_list[[rp]], dist_to_mangrove,
                                                    alpha=alpha0, L=L0, dmax_=dmax)
  names(scenario$with_mangroves[[rp]]) <- names(flood_list[[rp]])
}

# =========================
# 7) METHODOLOGY: EXPOSURE (PDF 4.3, Eq 6–11)
# =========================
exposure_by_rp <- function(depth_r, pop_r, ntl_r, polys_sf, h=h0) {
  flooded <- depth_r > h
  tibble::tibble(
    pop_exposed = exactextractr::exact_extract(pop_r * flooded, polys_sf, "sum"),
    ntl_exposed = exactextractr::exact_extract(ntl_r * flooded, polys_sf, "sum")
  )
}

rp_focus_chr <- as.character(RP_FOCUS)
exp_nomang <- exposure_by_rp(scenario$no_mangroves[[rp_focus_chr]], pop_m, ntl_m, phl_adm2_m, h=h0)
exp_with   <- exposure_by_rp(scenario$with_mangroves[[rp_focus_chr]], pop_m, ntl_m, phl_adm2_m, h=h0)

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

# =========================
# 8) METHODOLOGY: MONETARY VALUATION (PDF 4.4, Eq 12–17 + Table 2)
# =========================

# Eq (12) GDP(x) = NTL(x) * GDP_PHL / sum_{x in C} NTL(x)
gdp_year <- 2018
gdp_usd <- tryCatch({
  WDI::WDI(country="PHL", indicator="NY.GDP.MKTP.CD", start=gdp_year, end=gdp_year) |>
    dplyr::pull(NY.GDP.MKTP.CD) |> as.numeric()
}, error=function(e) NA_real_)
if (is.na(gdp_usd)) gdp_usd <- 376e9   # PDF value

ntl_sum_coast <- terra::global(ntl_m, "sum", na.rm=TRUE)[1,1]
gdp_per_ntl <- gdp_usd / max(ntl_sum_coast, 1e-9)
gdp_grid <- ntl_m * gdp_per_ntl
names(gdp_grid) <- "gdp_proxy_usd"

# Eq (13) Assets(x) = GDP(x) * (K/Y)
asset_value <- gdp_grid * K_Y
names(asset_value) <- "asset_usd"

# Table 2: Huizinga et al. (2017) Asia residential/mixed depth-damage control points
damage_frac <- function(depth_m) {
  approx(
    x = c(0.0, 0.1, 0.5, 1.0, 2.0, 4.0, 6.0),
    y = c(0.00, 0.07, 0.22, 0.38, 0.57, 0.80, 1.00),
    xout = depth_m,
    rule = 2
  )$y
}

# Eq (14) Damage_RP(x) = f(depth) * Assets(x)
damage_raster <- function(depth_r, asset_r) {
  frac <- terra::app(depth_r, fun = damage_frac)
  frac * asset_r
}

# Eq (15) D_{i,RP} = sum_{x in i} Damage_RP(x)
compute_damage_long <- function(depth_maps, polygons_sf, asset_r) {
  n_poly <- nrow(polygons_sf)
  purrr::map_dfr(names(depth_maps), function(rp_chr) {
    rp <- as.numeric(rp_chr)
    dmg_map <- damage_raster(depth_maps[[rp_chr]], asset_r)
    dmg_sum <- exactextractr::exact_extract(dmg_map, polygons_sf, "sum")
    rm(dmg_map); gc()
    tibble::tibble(poly_id = seq_len(n_poly), rp = rp, damage_usd = as.numeric(dmg_sum))
  })
}

# Eq (16) EAD integral: trapezoidal approximation at the 7 return periods (no extra p=1 anchor)
risk_integral_ead <- function(rp, damage) {
  p <- 1 / rp
  ord <- order(p, decreasing=TRUE)
  p <- p[ord]; damage <- damage[ord]
  
  # Tail extension to p=0 using max-RP damage (common in flood-risk practice)
  p <- c(p, 0)
  damage <- c(damage, tail(damage, 1))
  
  sum((damage[-length(damage)] + damage[-1]) / 2 * (p[-length(p)] - p[-1]))
}

damage_to_ead <- function(dmg_long) {
  dmg_long |>
    dplyr::group_by(poly_id) |>
    dplyr::summarise(EAD_usd = risk_integral_ead(rp, damage_usd), .groups="drop")
}

dmg_nomang_long <- compute_damage_long(scenario$no_mangroves,   phl_adm2_m, asset_value)
dmg_with_long   <- compute_damage_long(scenario$with_mangroves, phl_adm2_m, asset_value)

ead_nomang <- damage_to_ead(dmg_nomang_long) |> dplyr::rename(EAD_nomang_usd = EAD_usd)
ead_with   <- damage_to_ead(dmg_with_long)   |> dplyr::rename(EAD_with_usd   = EAD_usd)

# Eq (17) ΔEAD_i = EAD_no - EAD_with
valuation_tbl <- phl_adm2_m |>
  st_drop_geometry() |>
  dplyr::transmute(poly_id = dplyr::row_number(), GID_2, NAME_1, NAME_2) |>
  dplyr::left_join(ead_nomang, by="poly_id") |>
  dplyr::left_join(ead_with,   by="poly_id") |>
  dplyr::mutate(EAD_avoided_usd = EAD_nomang_usd - EAD_with_usd)

# =========================
# 9) CBA (PDF 4.5, Eq 18–20)
# =========================

# Eq (18) MangroveHa_i = sum MangroveFraction(x) * CellAreaHa(x)
cell_area_ha <- terra::cellSize(template, unit="ha")
mangrove_area_ha_r <- cell_area_ha * gmw_frac_1km
mangrove_ha_by_muni <- exactextractr::exact_extract(mangrove_area_ha_r, phl_adm2_m, "sum")

# Eq (19) NPV_i = sum_t ΔEAD_i*(1+g)^(t-1)/(1+r)^t
npv_annuity <- function(annual_benefit, r, T, g=0) {
  pv_factor <- sum(((1 + g)^(0:(T - 1))) / ((1 + r)^(1:T)))
  annual_benefit * pv_factor
}

cost_low  <- 1500
cost_mid  <- 4000
cost_high <- 10000

cba_tbl <- valuation_tbl |>
  dplyr::mutate(
    mangrove_ha = as.numeric(mangrove_ha_by_muni),
    NPV_benefits_usd = npv_annuity(EAD_avoided_usd, disc, T_years, g),
    # Eq (20) BCR = NPV / (CostPerHa * MangroveHa)
    NPV_cost_low  = mangrove_ha * cost_low,
    NPV_cost_mid  = mangrove_ha * cost_mid,
    NPV_cost_high = mangrove_ha * cost_high,
    BCR_low  = NPV_benefits_usd / pmax(NPV_cost_low,  1e-9),
    BCR_mid  = NPV_benefits_usd / pmax(NPV_cost_mid,  1e-9),
    BCR_high = NPV_benefits_usd / pmax(NPV_cost_high, 1e-9)
  )

final_tbl <- cba_tbl |>
  dplyr::left_join(exposure_tbl, by=c("GID_2","NAME_1","NAME_2"))

readr::write_csv(final_tbl, "outputs/municipality_mangrove_flood_cba.csv")

# --- Top 20 municipalities by BCR ---------------------------------------

top20 <- final_tbl |>
  dplyr::arrange(desc(BCR_low)) |>
  dplyr::slice(1:20)

readr::write_csv(top20, "outputs/top20_bcr_municipalities.csv")


library(dplyr)
library(readr)

top20_clean <- readr::read_csv(
  "outputs/top20_bcr_municipalities.csv",
  show_col_types = FALSE
) |>
  # remove only these ID columns
  dplyr::select(-any_of(c("poly_id", "GID_2"))) |>
  # remove only columns that are entirely 0 or NA
  dplyr::select(
    where(~ !(is.numeric(.x) && all(is.na(.x) | .x == 0)))
  )

readr::write_csv(
  top20_clean,
  "outputs/top20_bcr_municipalities_clean.csv"
)

# =========================
# 10) DISTRIBUTIONAL ANALYSIS (PDF 4.6)
# =========================
ntl_by_muni <- exactextractr::exact_extract(ntl_m, phl_adm2_m, "sum")
pop_by_muni <- exactextractr::exact_extract(pop_m, phl_adm2_m, "sum")

dist_tbl <- final_tbl |>
  dplyr::mutate(
    ntl_sum = as.numeric(ntl_by_muni),
    pop_sum = as.numeric(pop_by_muni),
    ntl_pc  = ntl_sum / pmax(pop_sum, 1e-9),
    benefit_share = EAD_avoided_usd / sum(EAD_avoided_usd, na.rm=TRUE),
    pop_share = pop_sum / sum(pop_sum, na.rm=TRUE)
  ) |>
  dplyr::arrange(ntl_pc) |>
  dplyr::mutate(cum_pop = cumsum(pop_share))

bottom40_benefit_share <- dist_tbl |>
  dplyr::filter(cum_pop <= 0.40) |>
  dplyr::summarise(share = sum(benefit_share, na.rm=TRUE)) |>
  dplyr::pull(share)

message(sprintf("Share of avoided EAD accruing to bottom 40%% (ntl/capita proxy): %.3f",
                bottom40_benefit_share))

# =========================
# 11) MAPS AND FIGURES
#   Simplified: fewer helper functions, more self-contained sections
# =========================

if (!requireNamespace("png", quietly = TRUE)) install.packages("png")

tmap::tmap_mode("plot")

# -------------------------
# 11.1 COMMON HELPERS
# -------------------------
save_tmap_png <- function(tm_obj, filename, width = 2200, height = 1400, res = 300) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  grDevices::png(filename, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(tm_obj)
}

thicken_plot_raster <- function(r, fact = 2, win = 3) {
  r_agg <- terra::aggregate(r, fact = fact, fun = "max", na.rm = TRUE)
  r_thick <- terra::focal(
    r_agg,
    w = matrix(1, nrow = win, ncol = win),
    fun = "max",
    na.rm = TRUE,
    fillvalue = NA
  )
  terra::ifel(r_thick <= 0, NA, r_thick)
}

make_exposure_raster <- function(depth_r, var_r, h = h0) {
  flooded <- depth_r > h
  var_r * flooded
}

# -------------------------
# 11.2 COMMON LAYERS / LAYOUT
# -------------------------
phl_outline <- phl_adm1_m |>
  sf::st_union() |>
  sf::st_as_sf() |>
  sf::st_simplify(dTolerance = 2000)

coastal_outline <- coastal_zone_m |>
  sf::st_union() |>
  sf::st_as_sf() |>
  sf::st_simplify(dTolerance = 1000)

layout_common <- tm_layout(
  frame = FALSE,
  legend.outside = TRUE,
  legend.outside.position = "right",
  legend.bg.color = "white",
  legend.bg.alpha = 0.90,
  legend.title.size = 1.00,
  legend.text.size = 0.80,
  inner.margins = c(0.02, 0.02, 0.02, 0.02)
)

comp_pos  <- c(0.88, 0.22)
scale_pos <- c(0.88, 0.10)

border_main_lwd <- 0.9
border_main_col <- "grey35"

# -------------------------
# 11.3 MAP 1: MANGROVES
# -------------------------
mangrove_plot <- terra::ifel(mangrove_presence, 1, NA)
names(mangrove_plot) <- "mangrove_presence"
mangrove_plot_thick <- thicken_plot_raster(mangrove_plot, fact = 2, win = 3)

p1 <- tm_shape(phl_outline) +
  tm_borders(lwd = 0.3, col = "white") +
  tm_shape(mangrove_plot_thick) +
  tm_raster(
    col.scale = tm_scale_categorical(values = "#39FF14", labels = "Present"),
    col.legend = tm_legend(title = "Mangroves"),
    col_alpha = 0.90
  ) +
  layout_common +
  tm_title("Mangrove presence", color = "white", size = .8) +
  tm_layout(
    legend.show = FALSE,
    bg.color = "black",
    outer.bg.color = "black",
    inner.margins = c(0, 0, 0, 0)
  )

{
  tmp_map <- tempfile(fileext = ".png")
  
  save_tmap_png(
    p1,
    tmp_map,
    width = 2200,
    height = 1600,
    res = 400
  )
  
  img <- png::readPNG(tmp_map)
  
  # ---- crop internal black margins --------------------------------------
  rgb_max <- apply(img[,,1:3, drop = FALSE], c(1, 2), max)
  keep_rows <- which(apply(rgb_max > 0.01, 1, any))
  keep_cols <- which(apply(rgb_max > 0.01, 2, any))
  
  if (length(keep_rows) > 0 && length(keep_cols) > 0) {
    pad <- 8
    r1 <- max(min(keep_rows) - pad, 1)
    r2 <- min(max(keep_rows) + pad, dim(img)[1])
    c1 <- max(min(keep_cols) - pad, 1)
    c2 <- min(max(keep_cols) + pad, dim(img)[2])
    img <- img[r1:r2, c1:c2, , drop = FALSE]
  }
  
  # ---- compute narrower final width, no stretching ----------------------
  final_height <- 1600
  map_height_frac <- 0.98
  map_drawn_px <- final_height * map_height_frac
  aspect_ratio <- ncol(img) / nrow(img)
  
  final_width <- ceiling(map_drawn_px * aspect_ratio + 80)
  final_width <- max(final_width, 1300)
  
  grDevices::png(
    "outputs/map_mangroves.png",
    width = final_width,
    height = final_height,
    res = 300,
    bg = "black"
  )
  
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "black", col = NA))
  
  grid::grid.raster(
    img,
    x = 0.5, y = 0.5,
    height = 0.98,
    just = "center",
    interpolate = FALSE
  )
  
  grDevices::dev.off()
  unlink(tmp_map)
  }

# -------------------------
# 11.7 MAP 5: EAD COMPARE
# -------------------------
phl_map_ead <- phl_adm2_m |>
  dplyr::left_join(final_tbl, by = c("GID_2", "NAME_1", "NAME_2"))

# friendlier facet names
phl_map_ead$No_mangroves <- phl_map_ead$EAD_nomang_usd / 1e4
phl_map_ead$With_mangroves <- phl_map_ead$EAD_with_usd / 1e4
phl_map_ead$Avoided <- pmax(phl_map_ead$No_mangroves - phl_map_ead$With_mangroves, 0)

ead_vals <- c(
  phl_map_ead$No_mangroves,
  phl_map_ead$With_mangroves,
  phl_map_ead$Avoided
)
ead_vals <- ead_vals[is.finite(ead_vals) & ead_vals >= 0]
ead_breaks <- c(0, 1, 10, 100, 1000, 5000, 20000, 50000, Inf)

# common colors for all three facets
n_ead <- length(ead_breaks) - 1
ead_cols <- grDevices::colorRampPalette(
  c("#d9ecff", "#9fc5f8", "#5b9bd5", "#0b78c8", "#084c8d")
)(n_ead)

ead_labels <- sapply(seq_len(length(ead_breaks) - 1), function(i) {
  paste0(
    format(round(ead_breaks[i], 1), big.mark = ",", scientific = FALSE, trim = TRUE),
    " - ",
    format(round(ead_breaks[i + 1], 1), big.mark = ",", scientific = FALSE, trim = TRUE)
  )
})

p_ead_compare <- tm_shape(phl_outline) +
  tm_borders(lwd = 0.8, col = "white") +
  tm_shape(phl_map_ead) +
  tm_polygons(
    c("No_mangroves", "With_mangroves", "Avoided"),
    style = "fixed",
    breaks = ead_breaks,
    palette = ead_cols,
    border.col = "grey65",
    lwd = 0.15,
    colorNA = "grey20",
    textNA = "Missing",
    title = "EAD (10k USD/yr)"
  ) +
  tm_facets(ncol = 3) +
  tm_title("Expected annual damages: without, with, and avoided", color = "white", size = 1.2) +
  tm_layout(
    frame = FALSE,
    legend.show = FALSE,
    bg.color = "black",
    outer.bg.color = "black",
    inner.margins = c(0, 0, 0, 0)
  )

{
  tmp_map <- tempfile(fileext = ".png")
  
  save_tmap_png(
    p_ead_compare,
    tmp_map,
    width = 2600,
    height = round(1600 * 0.93),
    res = 400
  )
  
  img <- png::readPNG(tmp_map)
  
  # ---- crop internal black border from temporary map --------------------
  rgb_max <- apply(img[,,1:3, drop = FALSE], c(1, 2), max)
  keep_rows <- which(apply(rgb_max > 0.01, 1, any))
  keep_cols <- which(apply(rgb_max > 0.01, 2, any))
  
  if (length(keep_rows) > 0 && length(keep_cols) > 0) {
    pad <- 8
    r1 <- max(min(keep_rows) - pad, 1)
    r2 <- min(max(keep_rows) + pad, dim(img)[1])
    c1 <- max(min(keep_cols) - pad, 1)
    c2 <- min(max(keep_cols) + pad, dim(img)[2])
    img <- img[r1:r2, c1:c2, , drop = FALSE]
  }
  
  # ---- make final canvas narrower, same style as RP10 map --------------
  final_height <- 1600
  top_frac <- 0.93
  map_height_frac <- 0.98
  
  top_panel_px <- final_height * top_frac
  map_drawn_px <- top_panel_px * map_height_frac
  
  aspect_ratio <- ncol(img) / nrow(img)
  
  final_width <- ceiling(map_drawn_px * aspect_ratio + 120)
  final_width <- max(final_width, 2000)
  
  grDevices::png(
    "outputs/map_ead_compare.png",
    width = final_width,
    height = final_height,
    res = 300,
    bg = "black"
  )
  
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "black", col = NA))
  
  lay <- grid::grid.layout(
    nrow = 2,
    ncol = 1,
    heights = grid::unit(c(0.93, 0.07), "npc")
  )
  grid::pushViewport(grid::viewport(layout = lay))
  
  # map
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.raster(
    img,
    x = 0.5, y = 0.5,
    height = 0.98,
    just = "center",
    interpolate = FALSE
  )
  grid::popViewport()
  
  # bottom legend
  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid::grid.text(
    "EAD (10k USD/yr)",
    x = 0.5, y = 0.72,
    just = c("center", "center"),
    gp = grid::gpar(fontsize = 9, col = "white")
  )
  
  legend_labels <- c(
    "0-50k", "50k-100k", "100k-150k", "150k-200k"
  )
  
  x_pos <- seq(0.08, 0.82, length.out = length(legend_labels))
  box_w <- 0.014
  box_h <- 0.20
  
  for (i in seq_along(legend_labels)) {
    grid::grid.roundrect(
      x = x_pos[i], y = 0.28,
      width = box_w, height = box_h,
      r = grid::unit(0.02, "snpc"),
      just = c("left", "center"),
      gp = grid::gpar(fill = ead_cols[i], col = ead_cols[i])
    )
    grid::grid.text(
      legend_labels[i],
      x = x_pos[i] + 0.020, y = 0.28,
      just = c("left", "center"),
      gp = grid::gpar(fontsize = 7, col = "white")
    )
  }
  
  grid::popViewport()
  grid::popViewport()
  grDevices::dev.off()
  unlink(tmp_map)
  }

# -------------------------
# 11.5 MAP 3: MULTI-RP FLOOD DARK
# -------------------------
flood_cols_blue <- c(
  "#9fe8ff", "#6fddff", "#3ed0ff", "#14c2ff",
  "#00adf0", "#0088e0", "#0060c7"
)

selected_rps <- c(2, 10, 50, 250)

rp_stack <- terra::rast(
  lapply(selected_rps, function(rp) {
    r <- scenario$no_mangroves[[as.character(rp)]]
    r <- terra::ifel(r <= 0, NA, r)
    r <- thicken_plot_raster(r, fact = 2, win = 3)
    #r
  })
)

names(rp_stack) <- paste0("RP ", selected_rps)

p2b_dark <- tm_shape(phl_outline) +
  tm_borders(lwd = 0.15, col = "white") +
  tm_shape(rp_stack) +
  tm_raster(
    col.scale = tm_scale_intervals(
      style = "fixed",
      breaks = flood_breaks,
      values = flood_cols_blue
    ),
    col.legend = tm_legend(show = FALSE, title = "Flood depth (m)")
  ) +
  tm_facets(ncol = 2, free.scales = FALSE) +
  tm_title("Flood hazard across return periods (no mangroves)", color = "white") +
  tm_layout(
    frame = FALSE,
    legend.show = FALSE,
    bg.color = "black",
    outer.bg.color = "black",
    inner.margins = c(0, 0, 0, 0)
  )

{
  tmp_map <- tempfile(fileext = ".png")
  
  save_tmap_png(
    p2b_dark,
    tmp_map,
    width = 2400,
    height = 1800,
    res = 400
  )
  
  img <- png::readPNG(tmp_map)
  
  # ---- crop internal black border from temporary map --------------------
  rgb_max <- apply(img[,,1:3, drop = FALSE], c(1, 2), max)
  keep_rows <- which(apply(rgb_max > 0.01, 1, any))
  keep_cols <- which(apply(rgb_max > 0.01, 2, any))
  
  if (length(keep_rows) > 0 && length(keep_cols) > 0) {
    pad <- 8
    r1 <- max(min(keep_rows) - pad, 1)
    r2 <- min(max(keep_rows) + pad, dim(img)[1])
    c1 <- max(min(keep_cols) - pad, 1)
    c2 <- min(max(keep_cols) + pad, dim(img)[2])
    img <- img[r1:r2, c1:c2, , drop = FALSE]
  }
  
  # ---- make final canvas narrower, without stretching map ---------------
  final_height <- 1800
  top_frac <- 0.94
  map_height_frac <- 0.98
  
  top_panel_px <- final_height * top_frac
  map_drawn_px <- top_panel_px * map_height_frac
  
  aspect_ratio <- ncol(img) / nrow(img)
  
  final_width <- ceiling(map_drawn_px * aspect_ratio + 120)
  
  # keep enough width so the legend fits comfortably
  final_width <- max(final_width, 1700)
  
  grDevices::png(
    "outputs/map_flood_multi_rp_dark.png",
    width = final_width,
    height = final_height,
    res = 400,
    bg = "black"
  )
  
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "black", col = NA))
  
  lay <- grid::grid.layout(
    nrow = 2,
    ncol = 1,
    heights = grid::unit(c(0.94, 0.06), "npc")
  )
  grid::pushViewport(grid::viewport(layout = lay))
  
  # map
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.raster(
    img,
    x = 0.5, y = 0.5,
    height = 0.98,
    just = "center",
    interpolate = FALSE
  )
  grid::popViewport()
  
  # legend
  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid::grid.text(
    "Flood depth (m)",
    x = 0.5, y = 0.72,
    just = c("center", "center"),
    gp = grid::gpar(fontsize = 8, col = "white")
  )
  
  legend_labels <- c(
    "0.0 - 0.4", "0.5 - 0.9", "1.0 - 1.9",
    "2.0 - 3.9", "4.0 - 5.9", "6.0 - 9.9", "≥ 10.0"
  )
  
  x_pos <- seq(0.08, 0.88, length.out = length(legend_labels))
  box_w <- 0.014
  box_h <- 0.16
  
  for (i in seq_along(legend_labels)) {
    grid::grid.roundrect(
      x = x_pos[i], y = 0.28,
      width = box_w, height = box_h,
      r = grid::unit(0.02, "snpc"),
      just = c("left", "center"),
      gp = grid::gpar(fill = flood_cols_blue[i], col = flood_cols_blue[i])
    )
    grid::grid.text(
      legend_labels[i],
      x = x_pos[i] + 0.018, y = 0.28,
      just = c("left", "center"),
      gp = grid::gpar(fontsize = 7, col = "white")
    )
  }
  
  grid::popViewport()
  grid::popViewport()
  grDevices::dev.off()
  unlink(tmp_map)
  }
# -------------------------
# 11.6 MAP 4: EXPOSURE DARK
# -------------------------
pop_nomang_map <- make_exposure_raster(scenario$no_mangroves[[rp_focus_chr]], pop_m)
pop_with_map   <- make_exposure_raster(scenario$with_mangroves[[rp_focus_chr]], pop_m)

ntl_nomang_map <- make_exposure_raster(scenario$no_mangroves[[rp_focus_chr]], ntl_m)
ntl_with_map   <- make_exposure_raster(scenario$with_mangroves[[rp_focus_chr]], ntl_m)

pop_diff_map <- pop_nomang_map - pop_with_map
ntl_diff_map <- ntl_nomang_map - ntl_with_map

pop_nomang_map <- terra::ifel(pop_nomang_map <= 0, NA, pop_nomang_map)
pop_with_map   <- terra::ifel(pop_with_map   <= 0, NA, pop_with_map)
pop_diff_map   <- terra::ifel(pop_diff_map   <= 0, NA, pop_diff_map)

ntl_nomang_map <- terra::ifel(ntl_nomang_map <= 0, NA, ntl_nomang_map)
ntl_with_map   <- terra::ifel(ntl_with_map   <= 0, NA, ntl_with_map)
ntl_diff_map   <- terra::ifel(ntl_diff_map   <= 0, NA, ntl_diff_map)

pop_nomang_plot <- thicken_plot_raster(pop_nomang_map, fact = 2, win = 5)
pop_with_plot   <- thicken_plot_raster(pop_with_map,   fact = 2, win = 5)
pop_diff_plot   <- thicken_plot_raster(pop_diff_map,   fact = 2, win = 5)

ntl_nomang_plot <- thicken_plot_raster(ntl_nomang_map, fact = 2, win = 5)
ntl_with_plot   <- thicken_plot_raster(ntl_with_map,   fact = 2, win = 5)
ntl_diff_plot   <- thicken_plot_raster(ntl_diff_map,   fact = 2, win = 5)

pop_stack <- c(pop_nomang_plot, pop_with_plot, pop_diff_plot)
names(pop_stack) <- c("Pop NoMang", "Pop WithMang", "Pop Avoided")

ntl_stack <- c(ntl_nomang_plot, ntl_with_plot, ntl_diff_plot)
names(ntl_stack) <- c("NTL NoMang", "NTL WithMang", "NTL Avoided")

# --- Population breaks and labels computed locally here ------------------
pop_vals <- terra::values(pop_stack, mat = FALSE, na.rm = TRUE)
pop_vals <- pop_vals[is.finite(pop_vals) & pop_vals > 0]

if (length(pop_vals) == 0) {
  pop_breaks <- c(0, 1)
} else {
  pop_breaks <- as.numeric(stats::quantile(
    pop_vals,
    probs = c(0, 0.25, 0.50, 0.75, 0.90, 1.00),
    na.rm = TRUE,
    names = FALSE
  ))
  pop_breaks <- unique(pop_breaks)
  
  if (length(pop_breaks) < 2) {
    pop_breaks <- c(0, max(pop_vals, na.rm = TRUE))
  }
  
  if (length(pop_breaks) < 6) {
    pop_breaks <- pretty(c(0, max(pop_vals, na.rm = TRUE)), n = 5)
    pop_breaks <- unique(pop_breaks)
  }
  
  pop_breaks[1] <- 0
  
  if (tail(pop_breaks, 1) < max(pop_vals, na.rm = TRUE)) {
    pop_breaks <- c(pop_breaks, max(pop_vals, na.rm = TRUE))
  }
  
  pop_breaks <- unique(pop_breaks)
}

pop_labels <- sapply(seq_len(length(pop_breaks) - 1), function(i) {
  paste0(
    format(round(pop_breaks[i], 1), big.mark = ",", scientific = FALSE, trim = TRUE),
    " - ",
    format(round(pop_breaks[i + 1], 1), big.mark = ",", scientific = FALSE, trim = TRUE)
  )
})

# --- NTL breaks and labels computed locally here -------------------------
ntl_vals <- terra::values(ntl_stack, mat = FALSE, na.rm = TRUE)
ntl_vals <- ntl_vals[is.finite(ntl_vals) & ntl_vals > 0]

if (length(ntl_vals) == 0) {
  ntl_breaks <- c(0, 1)
} else {
  ntl_breaks <- as.numeric(stats::quantile(
    ntl_vals,
    probs = c(0, 0.25, 0.50, 0.75, 0.90, 1.00),
    na.rm = TRUE,
    names = FALSE
  ))
  ntl_breaks <- unique(ntl_breaks)
  
  if (length(ntl_breaks) < 2) {
    ntl_breaks <- c(0, max(ntl_vals, na.rm = TRUE))
  }
  
  if (length(ntl_breaks) < 6) {
    ntl_breaks <- pretty(c(0, max(ntl_vals, na.rm = TRUE)), n = 5)
    ntl_breaks <- unique(ntl_breaks)
  }
  
  ntl_breaks[1] <- 0
  
  if (tail(ntl_breaks, 1) < max(ntl_vals, na.rm = TRUE)) {
    ntl_breaks <- c(ntl_breaks, max(ntl_vals, na.rm = TRUE))
  }
  
  ntl_breaks <- unique(ntl_breaks)
}

ntl_labels <- sapply(seq_len(length(ntl_breaks) - 1), function(i) {
  paste0(
    format(round(ntl_breaks[i], 1), big.mark = ",", scientific = FALSE, trim = TRUE),
    " - ",
    format(round(ntl_breaks[i + 1], 1), big.mark = ",", scientific = FALSE, trim = TRUE)
  )
})

exp_cols <- c("#ffb3ff", "#ff80ff", "#ff4dff", "#e600e6", "#730073")

layout_exp_dark <- tm_layout(
  frame = FALSE,
  legend.show = FALSE,
  bg.color = "black",
  outer.bg.color = "black",
  inner.margins = c(0.005, 0.005, 0.005, 0.005)
)

p_exp_pop <- tm_shape(phl_outline) +
  tm_borders(lwd = 0.3, col = "white") +
  tm_shape(pop_stack) +
  tm_raster(
    col.scale = tm_scale_intervals(
      style = "fixed",
      breaks = pop_breaks,
      values = exp_cols,
      labels = pop_labels
    ),
    col.legend = tm_legend(title = "Exposed population")
  ) +
  tm_facets(ncol = 3, free.scales = FALSE) +
  layout_exp_dark

p_exp_ntl <- tm_shape(phl_outline) +
  tm_borders(lwd = 0.3, col = "white") +
  tm_shape(ntl_stack) +
  tm_raster(
    col.scale = tm_scale_intervals(
      style = "fixed",
      breaks = ntl_breaks,
      values = exp_cols,
      labels = ntl_labels
    ),
    col.legend = tm_legend(title = "Exposed NTL")
  ) +
  tm_facets(ncol = 3, free.scales = FALSE) +
  layout_exp_dark

tmp_pop <- tempfile(fileext = ".png")
tmp_ntl <- tempfile(fileext = ".png")

# common layout settings for each exposure row
row_height <- 950
legend_frac <- 0.16
map_frac <- 0.84
map_height_frac <- 0.96

# -------------------------
# population exposure row
# -------------------------
{
  tmp_map <- tempfile(fileext = ".png")
  
  save_tmap_png(
    p_exp_pop,
    tmp_map,
    width = round(3200 * 0.72),
    height = row_height,
    res = 300
  )
  
  img <- png::readPNG(tmp_map)
  
  # ---- crop internal black border from temporary map --------------------
  rgb_max <- apply(img[,,1:3, drop = FALSE], c(1, 2), max)
  keep_rows <- which(apply(rgb_max > 0.01, 1, any))
  keep_cols <- which(apply(rgb_max > 0.01, 2, any))
  
  if (length(keep_rows) > 0 && length(keep_cols) > 0) {
    pad <- 8
    r1 <- max(min(keep_rows) - pad, 1)
    r2 <- min(max(keep_rows) + pad, dim(img)[1])
    c1 <- max(min(keep_cols) - pad, 1)
    c2 <- min(max(keep_cols) + pad, dim(img)[2])
    img <- img[r1:r2, c1:c2, , drop = FALSE]
  }
  
  # ---- compute narrower final row width ---------------------------------
  aspect_ratio <- ncol(img) / nrow(img)
  map_drawn_px <- row_height * map_height_frac
  row_width_pop <- ceiling((map_drawn_px * aspect_ratio + 30) / map_frac)
  row_width_pop <- max(row_width_pop, 2400)
  
  grDevices::png(tmp_pop, width = row_width_pop, height = row_height, res = 300, bg = "black")
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "black", col = NA))
  
  lay <- grid::grid.layout(
    nrow = 1, ncol = 2,
    widths = grid::unit(c(legend_frac, map_frac), "npc")
  )
  grid::pushViewport(grid::viewport(layout = lay))
  
  # legend panel
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.roundrect(
    x = 0.05, y = 0.50,
    width = 0.90, height = 0.58,
    just = c("left", "center"),
    r = grid::unit(0.02, "snpc"),
    gp = grid::gpar(fill = "black", col = "white", lwd = 1)
  )
  grid::grid.text(
    "Exposed population",
    x = 0.12, y = 0.76,
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 10, col = "white")
  )
  
  y_pos <- seq(0.64, 0.26, length.out = length(pop_labels))
  for (i in seq_along(pop_labels)) {
    grid::grid.roundrect(
      x = 0.12, y = y_pos[i],
      width = 0.14, height = 0.07,
      just = c("left", "center"),
      r = grid::unit(0.02, "snpc"),
      gp = grid::gpar(fill = exp_cols[i], col = exp_cols[i])
    )
    grid::grid.text(
      pop_labels[i],
      x = 0.33, y = y_pos[i],
      just = c("left", "center"),
      gp = grid::gpar(fontsize = 8, col = "white")
    )
  }
  grid::popViewport()
  
  # map panel
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::grid.raster(
    img,
    x = 0.49, y = 0.46,
    height = map_height_frac,
    just = "center",
    interpolate = FALSE
  )
  grid::popViewport()
  
  grid::popViewport()
  grDevices::dev.off()
  unlink(tmp_map)
}

# -------------------------
# NTL exposure row
# -------------------------
{
  tmp_map <- tempfile(fileext = ".png")
  
  save_tmap_png(
    p_exp_ntl,
    tmp_map,
    width = round(3200 * 0.72),
    height = row_height,
    res = 300
  )
  
  img <- png::readPNG(tmp_map)
  
  # ---- crop internal black border from temporary map --------------------
  rgb_max <- apply(img[,,1:3, drop = FALSE], c(1, 2), max)
  keep_rows <- which(apply(rgb_max > 0.01, 1, any))
  keep_cols <- which(apply(rgb_max > 0.01, 2, any))
  
  if (length(keep_rows) > 0 && length(keep_cols) > 0) {
    pad <- 8
    r1 <- max(min(keep_rows) - pad, 1)
    r2 <- min(max(keep_rows) + pad, dim(img)[1])
    c1 <- max(min(keep_cols) - pad, 1)
    c2 <- min(max(keep_cols) + pad, dim(img)[2])
    img <- img[r1:r2, c1:c2, , drop = FALSE]
  }
  
  # ---- compute narrower final row width ---------------------------------
  aspect_ratio <- ncol(img) / nrow(img)
  map_drawn_px <- row_height * map_height_frac
  row_width_ntl <- ceiling((map_drawn_px * aspect_ratio + 30) / map_frac)
  row_width_ntl <- max(row_width_ntl, 2400)
  
  grDevices::png(tmp_ntl, width = row_width_ntl, height = row_height, res = 300, bg = "black")
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "black", col = NA))
  
  lay <- grid::grid.layout(
    nrow = 1, ncol = 2,
    widths = grid::unit(c(legend_frac, map_frac), "npc")
  )
  grid::pushViewport(grid::viewport(layout = lay))
  
  # legend panel
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.roundrect(
    x = 0.05, y = 0.50,
    width = 0.90, height = 0.58,
    just = c("left", "center"),
    r = grid::unit(0.02, "snpc"),
    gp = grid::gpar(fill = "black", col = "white", lwd = 1)
  )
  grid::grid.text(
    "Exposed NTL",
    x = 0.12, y = 0.76,
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 10, col = "white")
  )
  
  y_pos <- seq(0.64, 0.26, length.out = length(ntl_labels))
  for (i in seq_along(ntl_labels)) {
    grid::grid.roundrect(
      x = 0.12, y = y_pos[i],
      width = 0.14, height = 0.07,
      just = c("left", "center"),
      r = grid::unit(0.02, "snpc"),
      gp = grid::gpar(fill = exp_cols[i], col = exp_cols[i])
    )
    grid::grid.text(
      ntl_labels[i],
      x = 0.33, y = y_pos[i],
      just = c("left", "center"),
      gp = grid::gpar(fontsize = 8, col = "white")
    )
  }
  grid::popViewport()
  
  # map panel
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::grid.raster(
    img,
    x = 0.49, y = 0.46,
    height = map_height_frac,
    just = "center",
    interpolate = FALSE
  )
  grid::popViewport()
  
  grid::popViewport()
  grDevices::dev.off()
  unlink(tmp_map)
}

# -------------------------
# combine both rows
# -------------------------
{
  img_top <- png::readPNG(tmp_pop)
  img_bottom <- png::readPNG(tmp_ntl)
  
  final_width <- max(dim(img_top)[2], dim(img_bottom)[2])
  
  grDevices::png("outputs/map_exposure_dark.png", width = final_width, height = 1900, res = 300, bg = "black")
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "black", col = NA))
  
  lay <- grid::grid.layout(
    nrow = 2,
    ncol = 1,
    heights = grid::unit(c(0.5, 0.5), "npc")
  )
  grid::pushViewport(grid::viewport(layout = lay))
  
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.raster(
    img_top,
    x = 0.5, y = 0.5,
    width = 1, height = 1,
    just = "center",
    interpolate = FALSE
  )
  grid::popViewport()
  
  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid::grid.raster(
    img_bottom,
    x = 0.5, y = 0.5,
    width = 1, height = 1,
    just = "center",
    interpolate = FALSE
  )
  grid::popViewport()
  
  grid::popViewport()
  grDevices::dev.off()
}

unlink(c(tmp_pop, tmp_ntl))

# -------------------------
# 11.7 MAP 5: EAD COMPARE
#   Separate breaks per panel, readable legends in USD/yr
# -------------------------
phl_map_ead <- phl_adm2_m |>
  dplyr::left_join(final_tbl, by = c("GID_2", "NAME_1", "NAME_2"))

# use RAW USD values (no division by 1e4)
phl_map_ead$No_mangroves   <- phl_map_ead$EAD_nomang_usd
phl_map_ead$With_mangroves <- phl_map_ead$EAD_with_usd
phl_map_ead$Avoided        <- pmax(phl_map_ead$EAD_nomang_usd - phl_map_ead$EAD_with_usd, 0)

# ---- panel-specific breaks ---------------------------------------------
calc_panel_breaks <- function(x) {
  x <- x[is.finite(x) & x >= 0]
  if (length(x) == 0) return(c(0, 1))
  
  qs <- as.numeric(stats::quantile(
    x,
    probs = c(0, 0.50, 0.90, 0.95, 0.99, 1.00),
    na.rm = TRUE,
    names = FALSE
  ))
  qs <- unique(qs)
  
  if (length(qs) < 3) {
    qs <- unique(pretty(c(0, max(x, na.rm = TRUE)), n = 5))
  }
  
  qs[1] <- 0
  
  if (tail(qs, 1) < max(x, na.rm = TRUE)) {
    qs <- c(qs, max(x, na.rm = TRUE))
  }
  
  unique(qs)
}

format_usd_short <- function(v) {
  if (v >= 1e9) {
    paste0(round(v / 1e9, 1), "bn")
  } else if (v >= 1e6) {
    paste0(round(v / 1e6, 1), "m")
  } else if (v >= 1e3) {
    paste0(round(v / 1e3, 0), "k")
  } else {
    as.character(round(v, 0))
  }
}

make_panel_labels <- function(brks) {
  labs <- sapply(seq_len(length(brks) - 1), function(i) {
    lo <- format_usd_short(brks[i])
    hi <- format_usd_short(brks[i + 1])
    paste0(lo, "-", hi)
  })
  labs[length(labs)] <- paste0(format_usd_short(brks[length(brks) - 1]), "+")
  labs
}

brks_nomang <- calc_panel_breaks(phl_map_ead$No_mangroves)
brks_with   <- calc_panel_breaks(phl_map_ead$With_mangroves)
brks_avoid  <- calc_panel_breaks(phl_map_ead$Avoided)

labs_nomang <- make_panel_labels(brks_nomang)
labs_with   <- make_panel_labels(brks_with)
labs_avoid  <- make_panel_labels(brks_avoid)

cols_nomang <- grDevices::colorRampPalette(
  c("#d9ecff", "#9fc5f8", "#5b9bd5", "#0b78c8", "#084c8d")
)(length(brks_nomang) - 1)

cols_with <- grDevices::colorRampPalette(
  c("#d9ecff", "#9fc5f8", "#5b9bd5", "#0b78c8", "#084c8d")
)(length(brks_with) - 1)

cols_avoid <- grDevices::colorRampPalette(
  c("#f3e5ff", "#d8b4fe", "#c084fc", "#a855f7", "#7e22ce")
)(length(brks_avoid) - 1)

# ---- build separate maps -----------------------------------------------
ead_map_layout <- tm_layout(
  frame = FALSE,
  legend.show = FALSE,
  bg.color = "black",
  outer.bg.color = "black",
  inner.margins = c(0, 0, 0, 0)
)

p_ead_nomang <- tm_shape(phl_outline) +
  tm_borders(lwd = 0.8, col = "white") +
  tm_shape(phl_map_ead) +
  tm_polygons(
    "No_mangroves",
    style = "fixed",
    breaks = brks_nomang,
    palette = cols_nomang,
    border.col = "grey65",
    lwd = 0.15,
    colorNA = "grey20",
    textNA = "Missing",
    title = "EAD (USD/yr)"
  ) +
  tm_title("No mangroves", color = "white", size = 1.0) +
  ead_map_layout

p_ead_with <- tm_shape(phl_outline) +
  tm_borders(lwd = 0.8, col = "white") +
  tm_shape(phl_map_ead) +
  tm_polygons(
    "With_mangroves",
    style = "fixed",
    breaks = brks_with,
    palette = cols_with,
    border.col = "grey65",
    lwd = 0.15,
    colorNA = "grey20",
    textNA = "Missing",
    title = "EAD (USD/yr)"
  ) +
  tm_title("With mangroves", color = "white", size = 1.0) +
  ead_map_layout

p_ead_avoid <- tm_shape(phl_outline) +
  tm_borders(lwd = 0.8, col = "white") +
  tm_shape(phl_map_ead) +
  tm_polygons(
    "Avoided",
    style = "fixed",
    breaks = brks_avoid,
    palette = cols_avoid,
    border.col = "grey65",
    lwd = 0.15,
    colorNA = "grey20",
    textNA = "Missing",
    title = "EAD (USD/yr)"
  ) +
  tm_title("Avoided", color = "white", size = 1.0) +
  ead_map_layout

# ---- render temporary panel maps ---------------------------------------
tmp_nomang <- tempfile(fileext = ".png")
tmp_with   <- tempfile(fileext = ".png")
tmp_avoid  <- tempfile(fileext = ".png")

save_tmap_png(p_ead_nomang, tmp_nomang, width = 900, height = 1150, res = 300)
save_tmap_png(p_ead_with,   tmp_with,   width = 900, height = 1150, res = 300)
save_tmap_png(p_ead_avoid,  tmp_avoid,  width = 900, height = 1150, res = 300)

img_nomang <- png::readPNG(tmp_nomang)
img_with   <- png::readPNG(tmp_with)
img_avoid  <- png::readPNG(tmp_avoid)

# ---- final combined figure ---------------------------------------------
grDevices::png(
  "outputs/map_ead_compare.png",
  width = 3000,
  height = 1600,
  res = 300,
  bg = "black"
)

grid::grid.newpage()
grid::grid.rect(gp = grid::gpar(fill = "black", col = NA))

outer_lay <- grid::grid.layout(
  nrow = 2,
  ncol = 1,
  heights = grid::unit(c(0.83, 0.17), "npc")
)
grid::pushViewport(grid::viewport(layout = outer_lay))

# ---- top section: title + 3 maps ---------------------------------------
grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))

grid::grid.text(
  "Expected annual damages",
  x = 0.5, y = 0.98,
  just = c("center", "top"),
  gp = grid::gpar(fontsize = 18, col = "white")
)

map_lay <- grid::grid.layout(
  nrow = 1,
  ncol = 3,
  widths = grid::unit(c(1, 1, 1), "null")
)
grid::pushViewport(grid::viewport(
  x = 0.5, y = 0.44, width = 0.98, height = 0.82,
  just = "center",
  layout = map_lay
))

grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
grid::grid.raster(img_nomang, x = 0.5, y = 0.5, width = 0.96, height = 0.96, just = "center", interpolate = FALSE)
grid::popViewport()

grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
grid::grid.raster(img_with, x = 0.5, y = 0.5, width = 0.96, height = 0.96, just = "center", interpolate = FALSE)
grid::popViewport()

grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
grid::grid.raster(img_avoid, x = 0.5, y = 0.5, width = 0.96, height = 0.96, just = "center", interpolate = FALSE)
grid::popViewport()

grid::popViewport()
grid::popViewport()

# ---- bottom section: 3 separate legends --------------------------------
grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))

legend_lay <- grid::grid.layout(
  nrow = 1,
  ncol = 3,
  widths = grid::unit(c(1, 1, 1), "null")
)
grid::pushViewport(grid::viewport(layout = legend_lay))

draw_panel_legend <- function(col_id, title_txt, legend_labels, legend_cols) {
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = col_id))
  grid::grid.text(
    title_txt,
    x = 0.5, y = 0.92,
    just = c("center", "center"),
    gp = grid::gpar(fontsize = 9, col = "white")
  )
  
  x_pos <- seq(0.06, 0.84, length.out = length(legend_labels))
  box_w <- min(0.018, 0.70 / length(legend_labels))
  box_h <- 0.16
  txt_dx <- box_w + 0.010
  
  for (i in seq_along(legend_labels)) {
    grid::grid.roundrect(
      x = x_pos[i], y = 0.48,
      width = box_w, height = box_h,
      r = grid::unit(0.02, "snpc"),
      just = c("left", "center"),
      gp = grid::gpar(fill = legend_cols[i], col = legend_cols[i])
    )
    grid::grid.text(
      legend_labels[i],
      x = x_pos[i] + txt_dx, y = 0.48,
      just = c("left", "center"),
      gp = grid::gpar(fontsize = 6.5, col = "white")
    )
  }
  grid::popViewport()
}

draw_panel_legend(1, "No mangroves (USD/yr)", labs_nomang, cols_nomang)
draw_panel_legend(2, "With mangroves (USD/yr)", labs_with, cols_with)
draw_panel_legend(3, "Avoided (USD/yr)", labs_avoid, cols_avoid)

grid::popViewport()
grid::popViewport()

grid::popViewport()
grDevices::dev.off()

unlink(c(tmp_nomang, tmp_with, tmp_avoid))

# -------------------------
# 11.8 MAP 6: AVOIDED EAD
# -------------------------
phl_map_avoided <- phl_adm2_m |>
  dplyr::left_join(final_tbl, by = c("GID_2", "NAME_1", "NAME_2"))

# use raw USD
phl_map_avoided$Avoided <- pmax(phl_map_avoided$EAD_avoided_usd, 0)

# ---- local helper functions --------------------------------------------
calc_panel_breaks <- function(x) {
  x <- x[is.finite(x) & x >= 0]
  if (length(x) == 0) return(c(0, 1))
  
  qs <- as.numeric(stats::quantile(
    x,
    probs = c(0, 0.50, 0.90, 0.95, 0.99, 1.00),
    na.rm = TRUE,
    names = FALSE
  ))
  qs <- unique(qs)
  
  if (length(qs) < 3) {
    qs <- unique(pretty(c(0, max(x, na.rm = TRUE)), n = 5))
  }
  
  qs[1] <- 0
  
  if (tail(qs, 1) < max(x, na.rm = TRUE)) {
    qs <- c(qs, max(x, na.rm = TRUE))
  }
  
  unique(qs)
}

format_usd_short <- function(v) {
  if (v >= 1e9) {
    paste0(round(v / 1e9, 1), "bn")
  } else if (v >= 1e6) {
    paste0(round(v / 1e6, 1), "m")
  } else if (v >= 1e3) {
    paste0(round(v / 1e3, 0), "k")
  } else {
    as.character(round(v, 0))
  }
}

make_panel_labels <- function(brks) {
  labs <- sapply(seq_len(length(brks) - 1), function(i) {
    lo <- format_usd_short(brks[i])
    hi <- format_usd_short(brks[i + 1])
    paste0(lo, "-", hi)
  })
  labs[length(labs)] <- paste0(format_usd_short(brks[length(brks) - 1]), "+")
  labs
}

# ---- avoided-specific breaks / labels / colors -------------------------
brks_avoid <- calc_panel_breaks(phl_map_avoided$Avoided)
labs_avoid <- make_panel_labels(brks_avoid)

cols_avoid <- grDevices::colorRampPalette(
  c("#f3e5ff", "#d8b4fe", "#c084fc", "#a855f7", "#7e22ce")
)(length(brks_avoid) - 1)

# ---- build map ----------------------------------------------------------
p3 <- tm_shape(phl_outline) +
  tm_borders(lwd = 0.8, col = "white") +
  tm_shape(phl_map_avoided) +
  tm_polygons(
    "Avoided",
    style = "fixed",
    breaks = brks_avoid,
    palette = cols_avoid,
    border.col = "grey65",
    lwd = 0.15,
    colorNA = "grey20",
    textNA = "Missing",
    title = "Avoided EAD (USD/yr)"
  ) +
  tm_title(
    paste0("Avoided EAD (alpha=", alpha0, ", L=", L0, "m, K/Y=", K_Y, ")"),
    color = "white",
    size = 1.0
  ) +
  tm_layout(
    frame = FALSE,
    legend.show = FALSE,
    bg.color = "black",
    outer.bg.color = "black",
    inner.margins = c(0, 0, 0, 0)
  )

# ---- render temporary map ----------------------------------------------
{
  tmp_map <- tempfile(fileext = ".png")
  
  save_tmap_png(
    p3,
    tmp_map,
    width = 2200,
    height = round(1600 * 0.93),
    res = 400
  )
  
  img <- png::readPNG(tmp_map)
  
  # ---- crop internal black border --------------------------------------
  rgb_max <- apply(img[, , 1:3, drop = FALSE], c(1, 2), max)
  keep_rows <- which(apply(rgb_max > 0.01, 1, any))
  keep_cols <- which(apply(rgb_max > 0.01, 2, any))
  
  if (length(keep_rows) > 0 && length(keep_cols) > 0) {
    pad <- 8
    r1 <- max(min(keep_rows) - pad, 1)
    r2 <- min(max(keep_rows) + pad, dim(img)[1])
    c1 <- max(min(keep_cols) - pad, 1)
    c2 <- min(max(keep_cols) + pad, dim(img)[2])
    img <- img[r1:r2, c1:c2, , drop = FALSE]
  }
  
  # ---- compute compact final canvas ------------------------------------
  final_height <- 1600
  top_frac <- 0.93
  map_height_frac <- 0.98
  
  top_panel_px <- final_height * top_frac
  map_drawn_px <- top_panel_px * map_height_frac
  
  aspect_ratio <- ncol(img) / nrow(img)
  final_width <- ceiling(map_drawn_px * aspect_ratio + 120)
  final_width <- max(final_width, 1500)
  
  grDevices::png(
    "outputs/map_avoided_ead.png",
    width = final_width,
    height = final_height,
    res = 300,
    bg = "black"
  )
  
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "black", col = NA))
  
  lay <- grid::grid.layout(
    nrow = 2,
    ncol = 1,
    heights = grid::unit(c(0.93, 0.07), "npc")
  )
  grid::pushViewport(grid::viewport(layout = lay))
  
  # map
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.raster(
    img,
    x = 0.5, y = 0.5,
    height = 0.98,
    just = "center",
    interpolate = FALSE
  )
  grid::popViewport()
  
  # legend
  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid::grid.text(
    "Avoided EAD (USD/yr)",
    x = 0.5, y = 0.72,
    just = c("center", "center"),
    gp = grid::gpar(fontsize = 9, col = "white")
  )
  
  n_lab <- length(labs_avoid)
  x_pos <- seq(0.06, 0.90, length.out = n_lab)
  box_w <- min(0.018, 0.70 / max(n_lab, 1))
  box_h <- 0.20
  txt_dx <- box_w + 0.010
  
  for (i in seq_along(labs_avoid)) {
    grid::grid.roundrect(
      x = x_pos[i], y = 0.28,
      width = box_w, height = box_h,
      r = grid::unit(0.02, "snpc"),
      just = c("left", "center"),
      gp = grid::gpar(fill = cols_avoid[i], col = cols_avoid[i])
    )
    grid::grid.text(
      labs_avoid[i],
      x = x_pos[i] + txt_dx, y = 0.28,
      just = c("left", "center"),
      gp = grid::gpar(fontsize = 7, col = "white")
    )
  }
  
  grid::popViewport()
  grid::popViewport()
  grDevices::dev.off()
  unlink(tmp_map)
}

# -------------------------
# 11.9 MAP 7: BCR LOW
#   Option A: show BCR only where local mangroves and local cost exist
# -------------------------

crop_black_img <- function(img, threshold = 0.01, pad = 8) {
  rgb_max <- apply(img[, , 1:3, drop = FALSE], c(1, 2), max)
  keep_rows <- which(apply(rgb_max > threshold, 1, any))
  keep_cols <- which(apply(rgb_max > threshold, 2, any))
  
  if (length(keep_rows) == 0 || length(keep_cols) == 0) return(img)
  
  r1 <- max(min(keep_rows) - pad, 1)
  r2 <- min(max(keep_rows) + pad, dim(img)[1])
  c1 <- max(min(keep_cols) - pad, 1)
  c2 <- min(max(keep_cols) + pad, dim(img)[2])
  
  img[r1:r2, c1:c2, , drop = FALSE]
}

phl_map_bcr <- phl_adm2_m |>
  dplyr::left_join(final_tbl, by = c("GID_2", "NAME_1", "NAME_2"))

# Option A:
# only compute BCR where municipality has local mangroves and positive local cost
phl_map_bcr$BCR_low_A <- dplyr::if_else(
  phl_map_bcr$mangrove_ha > 0 & phl_map_bcr$NPV_cost_low > 0,
  phl_map_bcr$NPV_benefits_usd / phl_map_bcr$NPV_cost_low,
  NA_real_
)

# include NA category in legend
bcr_breaks <- c(0, 0.5, 1, 2, 5, 10, Inf)
bcr_labels <- c("0–0.5", "0.5–1", "1–2", "2–5", "5–10", "≥10", "No mangroves")
bcr_cols <- c(
  "#d8ffe1", "#a8ffb8", "#6dff8e",
  "#39ff6a", "#00c853", "#007a33",
  "grey20"
)

p_bcr <- tm_shape(phl_outline) +
  tm_borders(lwd = border_main_lwd, col = "white") +
  tm_shape(phl_map_bcr) +
  tm_polygons(
    "BCR_low_A",
    style = "fixed",
    breaks = bcr_breaks,
    palette = bcr_cols[1:6],   # only class colors go here
    border.col = "grey65",
    lwd = 0.15,
    colorNA = "grey20",
    textNA = "No mangroves",
    title = "BCR"
  ) +
  tm_title("BCR (low cost, local mangroves only)", color = "white", size = 1.0) +
  tm_layout(
    frame = FALSE,
    legend.show = FALSE,
    bg.color = "black",
    outer.bg.color = "black",
    inner.margins = c(0, 0, 0, 0)
  )

{
  tmp_map <- tempfile(fileext = ".png")
  
  save_tmap_png(
    p_bcr,
    tmp_map,
    width = 2200,
    height = round(1600 * 0.93),
    res = 400
  )
  
  img <- png::readPNG(tmp_map)
  img <- crop_black_img(img, threshold = 0.01, pad = 8)
  
  final_height <- 1600
  top_frac <- 0.93
  map_height_frac <- 0.98
  
  top_panel_px <- final_height * top_frac
  map_drawn_px <- top_panel_px * map_height_frac
  aspect_ratio <- ncol(img) / nrow(img)
  
  final_width <- ceiling(map_drawn_px * aspect_ratio + 120)
  final_width <- max(final_width, 1500)
  
  grDevices::png(
    "outputs/map_bcr_low.png",
    width = final_width,
    height = final_height,
    res = 300,
    bg = "black"
  )
  
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "black", col = NA))
  
  lay <- grid::grid.layout(
    nrow = 2,
    ncol = 1,
    heights = grid::unit(c(0.93, 0.07), "npc")
  )
  grid::pushViewport(grid::viewport(layout = lay))
  
  # map
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.raster(
    img,
    x = 0.5, y = 0.5,
    height = 0.98,
    just = "center",
    interpolate = FALSE
  )
  grid::popViewport()
  
  # bottom legend
  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid::grid.text(
    "BCR (low cost)",
    x = 0.5, y = 0.72,
    just = c("center", "center"),
    gp = grid::gpar(fontsize = 10, col = "white")
  )
  
  x_pos <- seq(0.04, 0.84, length.out = length(bcr_labels))
  box_w <- 0.017
  box_h <- 0.22
  
  for (i in seq_along(bcr_labels)) {
    grid::grid.roundrect(
      x = x_pos[i], y = 0.28,
      width = box_w, height = box_h,
      r = grid::unit(0.02, "snpc"),
      just = c("left", "center"),
      gp = grid::gpar(fill = bcr_cols[i], col = bcr_cols[i])
    )
    grid::grid.text(
      bcr_labels[i],
      x = x_pos[i] + 0.022, y = 0.28,
      just = c("left", "center"),
      gp = grid::gpar(fontsize = 7.3, col = "white")
    )
  }
  
  grid::popViewport()
  grid::popViewport()
  grDevices::dev.off()
  unlink(tmp_map)
  }



# -------------------------
# 11.10 SCATTER PLOTS
#   Harmonised font sizes and export sizes
# -------------------------

scatter_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "plain"),
    axis.title = element_text(size = 13),
    axis.text  = element_text(size = 10),
    panel.grid.minor = element_line(linewidth = 0.25),
    panel.grid.major = element_line(linewidth = 0.35)
  )

dist_tbl_ntl_plot <- dist_tbl |>
  dplyr::filter(pop_sum > 0)

p_ntl <- ggplot(dist_tbl_ntl_plot, aes(x = ntl_pc, y = EAD_avoided_usd)) +
  geom_point(alpha = 0.35, size = 1.6) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 0.01, 0.1, 1, 10),
    labels = function(x) scales::number(x, accuracy = 0.01, big.mark = ",")
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 1e5, 1e6, 1e7, 2e7, 4e7),
    labels = function(y) scales::number(y / 1e6, accuracy = 0.1, big.mark = ",")
  ) +
  labs(
    title = "Avoided EAD vs NTL per capita",
    x = "NTL per capita",
    y = "Avoided EAD (million USD)"
  ) +
  scatter_theme

dist_tbl_pop_plot <- dist_tbl |>
  dplyr::filter(pop_sum > 0)

p_pop <- ggplot(dist_tbl_pop_plot, aes(x = pop_sum, y = EAD_avoided_usd)) +
  geom_point(alpha = 0.5, size = 1.6) +
  scale_x_continuous(
    trans = "log1p",
    breaks = c(0, 1e4, 5e4, 1e5, 2e5, 5e5, 1e6),
    labels = function(x) scales::number(x / 1e4, accuracy = 0.1, big.mark = ",")
  ) +
  scale_y_continuous(
    trans = "log1p",
    breaks = c(0, 1e5, 1e6, 1e7, 2e7, 4e7),
    labels = function(y) scales::number(y / 1e6, accuracy = 0.1, big.mark = ",")
  ) +
  labs(
    title = "Avoided EAD vs population",
    x = "Population (10,000 persons)",
    y = "Avoided EAD (million USD)"
  ) +
  scatter_theme

dist_tbl_ntl_log <- dist_tbl |>
  dplyr::filter(pop_sum > 0, ntl_pc > 0, EAD_avoided_usd > 0) |>
  dplyr::mutate(
    log10_ntl_pc = log10(ntl_pc),
    log10_ead_avoided_usd = log10(EAD_avoided_usd)
  )

p_ntl_log <- ggplot(
  dist_tbl_ntl_log,
  aes(x = log10_ntl_pc, y = log10_ead_avoided_usd)
) +
  geom_point(alpha = 0.35, size = 1.6) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  labs(
    title = "Avoided EAD vs NTL per capita (log10 values)",
    x = "log10(NTL per capita)",
    y = "log10(Avoided EAD in USD)"
  ) +
  scatter_theme

dist_tbl_pop_log <- dist_tbl |>
  dplyr::filter(pop_sum > 0, EAD_avoided_usd > 0) |>
  dplyr::mutate(
    log10_pop_sum = log10(pop_sum),
    log10_ead_avoided_usd = log10(EAD_avoided_usd)
  )

p_pop_log <- ggplot(
  dist_tbl_pop_log,
  aes(x = log10_pop_sum, y = log10_ead_avoided_usd)
) +
  geom_point(alpha = 0.5, size = 1.6) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  labs(
    title = "Avoided EAD vs population (log10 values)",
    x = "log10(Population)",
    y = "log10(Avoided EAD in USD)"
  ) +
  scatter_theme

# same export size for all four figures
ggsave("outputs/dist_ntl.png",       p_ntl,     width = 7, height = 5.5, dpi = 300)
ggsave("outputs/dist_pop.png",       p_pop,     width = 7, height = 5.5, dpi = 300)
ggsave("outputs/dist_ntl_log10.png", p_ntl_log, width = 7, height = 5.5, dpi = 300)
ggsave("outputs/dist_pop_log10.png", p_pop_log, width = 7, height = 5.5, dpi = 300)

message("Done. Maps and figures saved in ./outputs/")

# =========================
# 12) SENSITIVITY (PDF Table 3 + discount-rate)
# =========================

# --- Helper: NATIONAL EAD total from depth maps (for Table 3)
# This computes total EAD over the coastal belt raster domain (not by municipality).
ead_total_national <- function(depth_maps, asset_r) {
  rp_vec <- as.numeric(names(depth_maps))
  dmg_tot <- numeric(length(rp_vec))
  for (k in seq_along(rp_vec)) {
    rp <- rp_vec[k]
    dmg_map <- damage_raster(depth_maps[[as.character(rp)]], asset_r)
    dmg_tot[k] <- terra::global(dmg_map, "sum", na.rm=TRUE)[1,1]
    rm(dmg_map); gc()
  }
  risk_integral_ead(rp_vec, dmg_tot)
}

if (RUN_SENS_ALPHA_L) {
  # Baseline national EAD (no mangroves) computed once
  ead_nat_nomang <- ead_total_national(flood_list, asset_value)
  
  # Grid over alpha-L (Table 3)
  grid <- tidyr::expand_grid(alpha = alpha_grid, L = L_grid)
  
  # For each combo, compute national EAD with mangroves, then avoided = baseline - with
  sens <- grid |>
    dplyr::mutate(
      ead_nat_with = purrr::pmap_dbl(list(alpha, L), function(alpha, L) {
        with_maps <- lapply(flood_list, \(r) apply_protection(r, dist_to_mangrove,
                                                              alpha=alpha, L=L, dmax_=dmax))
        names(with_maps) <- names(flood_list)
        ead_total_national(with_maps, asset_value)
      }),
      avoided_EAD_total_usd  = ead_nat_nomang - ead_nat_with,
      avoided_EAD_total_musd = avoided_EAD_total_usd / 1e6
    )
  
  readr::write_csv(sens, "outputs/sensitivity_avoided_ead_total.csv")
  
  # Table-3-ready wide matrix (USD million/year)
  table3 <- sens |>
    dplyr::select(alpha, L, avoided_EAD_total_musd) |>
    tidyr::pivot_wider(names_from=L, values_from=avoided_EAD_total_musd) |>
    dplyr::arrange(alpha)
  
  readr::write_csv(table3, "outputs/table3_avoided_ead_alpha_L_musd.csv")
}

if (RUN_SENS_DISC) {
  # Discount-rate sensitivity affects NPV/BCR (Eq 19–20)
  sens_disc <- purrr::map_dfr(disc_grid, function(rdisc) {
    pv_factor <- sum(((1 + g)^(0:(T_years - 1))) / ((1 + rdisc)^(1:T_years)))
    
    tmp <- valuation_tbl |>
      dplyr::mutate(
        mangrove_ha = as.numeric(mangrove_ha_by_muni),
        NPV_benefits_usd = EAD_avoided_usd * pv_factor,
        NPV_cost_mid = mangrove_ha * cost_mid,
        BCR_mid = NPV_benefits_usd / pmax(NPV_cost_mid, 1e-9)
      )
    
    tibble::tibble(
      discount_rate = rdisc,
      total_NPV_benefits_usd = sum(tmp$NPV_benefits_usd, na.rm=TRUE),
      median_BCR_mid = median(tmp$BCR_mid, na.rm=TRUE)
    )
  })
  
  readr::write_csv(sens_disc, "outputs/sensitivity_discount_rate.csv")
}

message("Done. Main output: outputs/municipality_mangrove_flood_cba.csv")
###############################################################################