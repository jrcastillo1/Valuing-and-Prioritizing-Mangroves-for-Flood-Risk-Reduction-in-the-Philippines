###############################################################################
# Mangrove Protection, Coastal Flood Risk, and Economic Exposure — Philippines
# CLEAN + INTERPRETABLE SCRIPT (matches PDF formulas, tables, parameters)
#
# Outputs:
#  - outputs/municipality_mangrove_flood_cba.csv
#  - outputs/sensitivity_avoided_ead_total.csv
#  - outputs/table3_avoided_ead_alpha_L_musd.csv
#  - outputs/sensitivity_discount_rate.csv
#  - outputs/map_mangroves.png
#  - outputs/map_flood_rp10.png
#  - outputs/map_avoided_ead.png
###############################################################################

# =========================
# 0) SETTINGS / SWITCHES
# =========================
MAKE_MAPS        <- TRUE
RUN_SENS_ALPHA_L <- TRUE   # Table 3 (alpha-L grid)
RUN_SENS_DISC    <- TRUE   # discount-rate grid for NPV/BCR
RP_FOCUS         <- 10      # for exposure + flood map
COASTAL_KM       <- 20      # coastal belt width in km (PDF: 20 km)

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
# 11) MAPS (paper-ready + compass/scale aligned under legend on the right)
#   MODS APPLIED (as suggested):
#     1) Reserve a larger right margin so legend sits in a clean “right column” (0.24 -> 0.30)
#     2) Nudge compass + scale a bit further right (0.86 -> 0.88)
#     3) Lighten borders slightly (less heavy ink)
# =========================
if (MAKE_MAPS) {
  tmap::tmap_mode("plot")
  
  # --- Prepare cleaner layers --------------------------------------------
  phl_outline <- phl_adm1_m |> sf::st_union() |> sf::st_as_sf() |> sf::st_simplify(dTolerance = 2000)
  coastal_outline <- coastal_zone_m |> sf::st_union() |> sf::st_as_sf() |> sf::st_simplify(dTolerance = 1000)
  
  mangrove_plot <- terra::ifel(mangrove_presence, 1, NA)
  names(mangrove_plot) <- "mangrove_presence"
  
  flood_plot <- scenario$no_mangroves[[rp_focus_chr]]
  flood_plot <- terra::ifel(flood_plot <= 0, NA, flood_plot)
  
  phl_map <- phl_adm2_m |>
    dplyr::left_join(final_tbl, by = c("GID_2","NAME_1","NAME_2"))
  phl_map$EAD_avoided_usd <- pmax(phl_map$EAD_avoided_usd, 0)
  
  # --- Layout: move legend to the right OUTSIDE (no overlap) ----------------
  layout_common <- tm_layout(
    frame = FALSE,
    legend.outside = TRUE,                 # <-- changed
    legend.outside.position = "right",     # <-- changed
    # legend.position removed (only relevant when legend.outside = FALSE)
    legend.bg.color = "white",
    legend.bg.alpha = 0.90,
    legend.title.size = 1.00,
    legend.text.size  = 0.80,
    main.title.size = 1.10,
    inner.margins = c(0.02, 0.02, 0.02, 0.02)   # <-- optional but recommended when legend is outside
  )
  
  # Place compass and scale in the reserved right margin, below the legend
  comp_pos  <- c(0.88, 0.22)   # <-- nudged right (was 0.86)
  scale_pos <- c(0.88, 0.10)   # <-- nudged right (was 0.86)
  
  # Border styling (lighter than before)
  border_main_lwd <- 0.9       # was 1.2
  border_main_col <- "grey35"
  border_coast_lwd <- 0.8      # was 1.0
  border_coast_col <- "grey55"
  
  # --- Map 1: Mangroves presence -----------------------------------------
  p1 <- tm_shape(phl_outline) +
    tm_borders(lwd = border_main_lwd, col = border_main_col) +
    tm_shape(coastal_outline) +
    tm_borders(lwd = border_coast_lwd, col = border_coast_col) +
    tm_shape(mangrove_plot) +
    tm_raster(
      title = "Mangroves\n(1km cells)",
      style = "cat",
      labels = "Present",
      alpha = 0.90
    ) +
    tm_compass(type = "8star", position = comp_pos, size = 1.7) +
    tm_scale_bar(position = scale_pos, text.size = 0.65) +
    layout_common +
    tm_layout(main.title = paste0("Mangrove presence in coastal belt (\\tau=", tau, ")"))
  
  # --- Map 2: Flood depth RP10 --------------------------------------------
  flood_breaks <- c(0, 0.5, 1, 2, 4, 6, 10, Inf)
  
  p2 <- tm_shape(phl_outline) +
    tm_borders(lwd = border_main_lwd, col = border_main_col) +
    tm_shape(coastal_outline) +
    tm_borders(lwd = border_coast_lwd, col = border_coast_col) +
    tm_shape(flood_plot) +
    tm_raster(
      title = paste0("Flood depth (m)\nRP", RP_FOCUS),
      style = "fixed",
      breaks = flood_breaks,
      alpha = 0.90
    ) +
    tm_compass(type = "8star", position = comp_pos, size = 1.7) +
    tm_scale_bar(position = scale_pos, text.size = 0.65) +
    layout_common +
    tm_layout(main.title = paste0("Coastal inundation depth (RP", RP_FOCUS, ", baseline)"))

  # --- Map 3: Multi-RP flood hazard -----------------------
  rp_maps <- lapply(c(2,10,50,250), function(rp) {
    r <- scenario$no_mangroves[[as.character(rp)]]
    r <- terra::ifel(r <= 0, NA, r)
    names(r) <- paste0("RP", rp)
    r
  })

  rp_stack <- terra::rast(rp_maps)

  p2b <- tm_shape(phl_outline) +
    tm_borders(lwd = border_main_lwd, col = border_main_col) +
    tm_shape(rp_stack) +
    tm_raster(
      style = "fixed",
      breaks = flood_breaks,
      title = "Flood depth (m)"
    ) +
    tm_facets(ncol = 2) +
    layout_common +
    tm_layout(main.title = "Flood hazard across return periods (no mangroves)")

  save_tmap_png(p2b, "outputs/map_flood_multi_rp.png")

  # --- Map 4: Exposure maps (population + NTL, with/without mangroves) ------------

  make_exposure_raster <- function(depth_r, var_r, h=h0) {
    flooded <- depth_r > h
    var_r * flooded
  }

  pop_nomang_map <- make_exposure_raster(scenario$no_mangroves[[rp_focus_chr]], pop_m)
  pop_with_map   <- make_exposure_raster(scenario$with_mangroves[[rp_focus_chr]], pop_m)

  ntl_nomang_map <- make_exposure_raster(scenario$no_mangroves[[rp_focus_chr]], ntl_m)
  ntl_with_map   <- make_exposure_raster(scenario$with_mangroves[[rp_focus_chr]], ntl_m)

  exp_stack <- c(pop_nomang_map, ntl_nomang_map, pop_with_map, ntl_with_map)
  names(exp_stack) <- c("Pop NoMang", "NTL NoMang", "Pop WithMang", "NTL WithMang")

  p_exp <- tm_shape(exp_stack) +
    tm_raster(style="quantile") +
    tm_facets(ncol=2) +
    layout_common +
    tm_layout(main.title="Exposure to flooding (population and economic activity)")

  save_tmap_png(p_exp, "outputs/map_exposure.png")

  # --- Map 5: Expected annual damages maps ---------------------------------------

  phl_map$EAD_nomang_10k <- phl_map$EAD_nomang_usd / 1e4
  phl_map$EAD_with_10k   <- phl_map$EAD_with_usd / 1e4

  p_ead_compare <- tm_shape(phl_outline) +
    tm_borders(lwd = border_main_lwd, col = border_main_col) +
    tm_shape(phl_map) +
    tm_polygons(c("EAD_nomang_10k","EAD_with_10k"),
                style="quantile",
                title="EAD (10k USD/yr)") +
    tm_facets(ncol=2) +
    layout_common +
    tm_layout(main.title="Expected annual damages: with vs without mangroves")

  save_tmap_png(p_ead_compare, "outputs/map_ead_compare.png")
  
  # --- Map 6: Avoided EAD (clean legend, fewer digits) ---------------------
  phl_map$EAD_avoided_10k <- phl_map$EAD_avoided_usd / 1e4
  ead10k_breaks <- c(0, 0.1, 1, 10, 50, 200, Inf)
  
  p3 <- tm_shape(phl_outline) +
    tm_borders(lwd = border_main_lwd, col = border_main_col) +
    tm_shape(phl_map) +
    tm_polygons(
      "EAD_avoided_10k",
      title = "Avoided EAD\n(10k USD/yr)",
      style = "fixed",
      breaks = ead10k_breaks,
      border.col = "grey75",
      lwd = 0.15,
      colorNA = "grey92",
      textNA = "Missing"
    ) +
    tm_compass(type = "8star", position = comp_pos, size = 1.7) +
    tm_scale_bar(position = scale_pos, text.size = 0.65) +
    layout_common +
    tm_layout(main.title = paste0("Avoided EAD (\\alpha=", alpha0, ", L=", L0, "m, K/Y=", K_Y, ")"))
  
  # --- Map 7: BCR map (low-cost scenario) ----------------------------------------

  phl_map$BCR_low <- final_tbl$BCR_low

  p_bcr <- tm_shape(phl_outline) +
    tm_borders(lwd = border_main_lwd, col = border_main_col) +
    tm_shape(phl_map) +
    tm_polygons("BCR_low",
                style="quantile",
                title="BCR (low cost)") +
    layout_common +
    tm_layout(main.title="Benefit-cost ratio (low restoration cost scenario)")

  save_tmap_png(p_bcr, "outputs/map_bcr_low.png")

  # --- NEW: Distributional scatter plots ---------------------------------------

  p_ntl <- ggplot(dist_tbl, aes(x = ntl_pc, y = EAD_avoided_usd)) +
    geom_point(alpha=0.5) +
    scale_x_continuous(trans="log1p") +
    scale_y_continuous(trans="log1p") +
    labs(title="Avoided EAD vs NTL per capita",
        x="NTL per capita",
        y="Avoided EAD (USD)")

  p_pop <- ggplot(dist_tbl, aes(x = pop_sum, y = EAD_avoided_usd)) +
    geom_point(alpha=0.5) +
    scale_x_continuous(trans="log1p") +
    scale_y_continuous(trans="log1p") +
    labs(title="Avoided EAD vs population",
        x="Population",
        y="Avoided EAD (USD)")

  ggsave("outputs/dist_ntl.png", p_ntl, width=6, height=5)
  ggsave("outputs/dist_pop.png", p_pop, width=6, height=5)

  # --- F) Save (base PNG) ----------------------------------------------------
  save_tmap_png <- function(tm_obj, filename, width = 2200, height = 1400, res = 300) {
    dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
    grDevices::png(filename, width = width, height = height, res = res)
    print(tm_obj)
    grDevices::dev.off()
  }
  
  save_tmap_png(p1, "outputs/map_mangroves.png")
  save_tmap_png(p2, "outputs/map_flood_rp10.png")
  save_tmap_png(p3, "outputs/map_avoided_ead.png")
}
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