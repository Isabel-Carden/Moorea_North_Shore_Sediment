# 04 — Lean paper figures (PI-scoped)
#
# Purpose: three main-text figures only — interpolation map, in-bay vs
#          out-of-bay contrast, and in-bay linear models (no grain size).
#
# Response: spreadsheet Fine Fraction
# In-bay definition: Cook's Bay (C*) + Opunohu Bay (O*); all other sites
#                    are out-of-bay (Settings: in_bay_prefixes).
#
# Inputs:
#   Data/sediment_raw_data_snapshot.csv
#   outputs/02_osm_moorea.rds   (basemap cache from script 02)
#
# Outputs (outputs/):
#   04_analysis_ready.csv
#   04_fig1_kriging_fine_fraction.png
#   04_fig2_in_vs_out_bay.png
#   04_fig2_ttest.csv
#   04_fig3_inbay_lm.png          (outlet-only; main text)
#   04_fig3_lm_results.csv
#   04_run_metadata.txt
#   04_summary.md
#
# Run from project root: Sediment_Project_2026/
# Note: script 03 remains the exploratory archive; this is the paper path.

# ⚙️ Settings ----
# =============================================================================

snapshot_file <- "Data/sediment_raw_data_snapshot.csv"
out_dir <- "outputs"
set.seed(42)

# In-bay = Cook's (C) + Opunohu (O); everything else = out-of-bay
in_bay_prefixes <- c("C", "O")

# Fig 3: "outlet" only for now; set TRUE later if PIs want shore too
fig3_include_shore <- FALSE

# Map / kriging (same island frame as scripts 02–03)
map_xlim <- c(-149.94, -149.74)
map_ylim <- c(-17.615, -17.465)
grid_n <- 80
osm_cache <- file.path(out_dir, "02_osm_moorea.rds")
utm_crs <- 32706 # WGS84 / UTM zone 6S
lagoon_hull_ratio <- 0.22
moorea_label <- data.frame(
  Long = -149.828,
  Lat = -17.533,
  name = "Mo'orea"
)

col_lagoon <- "#9FD4E0"
col_reef <- "#E8D4A8"
col_island <- "#8FBF6A"
col_ocean <- "#5E9FBF"
col_in <- "#2A6F7A"
col_out <- "#C4A574"

# 📁 Paths / packages ----
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(sf)
library(sp)
library(gstat)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "Pacific/Tahiti")

# 1️⃣ Load & analysis-ready table ----
# =============================================================================

raw <- read_csv(snapshot_file, show_col_types = FALSE)

names(raw) <- names(raw) |>
  str_replace("^\\s+", "") |>
  str_replace(
    "Grain Size <=63 Mass \\(g\\)",
    "Grain Size le63 Mass (g)"
  )

clean_num <- function(x) {
  x_chr <- str_trim(as.character(x))
  bad <- x_chr %in% c("", "N/A", "NA", "#DIV/0!", "#VALUE!", "#REF!", "#N/A")
  x_chr[bad] <- NA_character_
  suppressWarnings(as.numeric(x_chr))
}

bay_from_site <- function(site) {
  site <- str_trim(as.character(site))
  case_when(
    str_detect(site, "^C") ~ "C",
    str_detect(site, "^O") ~ "O",
    str_detect(site, "^H") ~ "H",
    str_detect(site, "^P") ~ "P",
    str_detect(site, "^G") ~ "G",
    str_detect(site, "^J") ~ "J",
    str_detect(site, "^F") ~ "F",
    str_detect(site, "^A") ~ "A",
    str_detect(site, "^B") ~ "B",
    str_detect(site, "^S") ~ "S",
    str_detect(site, "^M") ~ "M",
    TRUE ~ "other"
  )
}

dat <- raw |>
  mutate(
    across(
      c(
        Lat, Long,
        `Distance From Shore (m)`,
        `Distance From Nearest major outlet (m)`,
        `Fine Fraction`
      ),
      clean_num
    ),
    Name = str_trim(as.character(Name)),
    Site = str_trim(as.character(Site)),
    bay = bay_from_site(Site),
    dist_outlet_m = `Distance From Nearest major outlet (m)`,
    dist_shore_m = `Distance From Shore (m)`,
    fine_fraction = `Fine Fraction`,
    # WHY: PI definition — embayment vs open lagoon/reef
    bay_group = if_else(
      bay %in% in_bay_prefixes,
      "In bay (Cook's / Opunohu)",
      "Out of bay"
    ),
    bay_group = factor(
      bay_group,
      levels = c("In bay (Cook's / Opunohu)", "Out of bay")
    ),
    log1p_fine = log1p(fine_fraction)
  )

analysis <- dat |>
  filter(
    !is.na(Site), Site != "",
    !is.na(Lat), !is.na(Long),
    !is.na(fine_fraction), is.finite(fine_fraction),
    !is.na(dist_outlet_m),
    !is.na(dist_shore_m)
  ) |>
  select(
    Site, Name, bay, bay_group, Lat, Long,
    dist_outlet_m, dist_shore_m,
    fine_fraction, log1p_fine
  )

write_csv(analysis, file.path(out_dir, "04_analysis_ready.csv"))

# 2️⃣ Fig 1 — kriging map of Fine Fraction ----
# =============================================================================

has_sf <- function(x) {
  !is.null(x) && inherits(x, "sf") && nrow(x) > 0 && !all(st_is_empty(x))
}

if (!file.exists(osm_cache)) {
  stop(
    "Missing OSM cache: ", osm_cache,
    " — run scripts/02_moorea_site_map.R once to build it."
  )
}
osm_layers <- readRDS(osm_cache)
moorea_poly <- osm_layers$island
reefs <- osm_layers$reefs
motus <- osm_layers$motus

lagoon <- NULL
if (has_sf(osm_layers$lagoon_osm)) {
  lagoon <- osm_layers$lagoon_osm
} else if (has_sf(reefs) || has_sf(motus)) {
  parts <- Filter(has_sf, list(reefs, motus))
  unions <- lapply(parts, function(x) {
    st_union(st_transform(st_make_valid(x), utm_crs))
  })
  outer_u <- Reduce(st_union, unions)
  island_u <- st_transform(st_make_valid(moorea_poly), utm_crs) |>
    st_union()
  hull <- tryCatch(
    st_concave_hull(outer_u, ratio = lagoon_hull_ratio),
    error = function(e) st_convex_hull(outer_u)
  )
  lagoon_u <- st_make_valid(st_difference(st_make_valid(hull), island_u))
  if (!all(st_is_empty(lagoon_u))) {
    lagoon <- st_sf(geometry = st_transform(lagoon_u, 4326))
  }
}

add_water_basemap <- function(p) {
  if (has_sf(lagoon)) {
    p <- p + geom_sf(data = lagoon, fill = col_lagoon, color = NA, inherit.aes = FALSE)
  }
  if (has_sf(reefs)) {
    p <- p +
      geom_sf(
        data = reefs,
        fill = col_reef,
        color = "#C4A574",
        linewidth = 0.15,
        inherit.aes = FALSE
      )
  }
  p
}

add_island_overlay <- function(p) {
  if (has_sf(motus)) {
    p <- p +
      geom_sf(
        data = motus,
        fill = col_island,
        color = "grey25",
        linewidth = 0.2,
        inherit.aes = FALSE
      )
  }
  p +
    geom_sf(
      data = moorea_poly,
      fill = col_island,
      color = "grey20",
      linewidth = 0.4,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = moorea_label,
      aes(x = Long, y = Lat, label = name),
      fontface = "bold.italic",
      size = 4,
      color = "grey20",
      alpha = 0.9,
      inherit.aes = FALSE
    )
}

sp_pts <- analysis |>
  select(Long, Lat, fine_fraction) |>
  as.data.frame()
coordinates(sp_pts) <- ~ Long + Lat
proj4string(sp_pts) <- CRS("+proj=longlat +datum=WGS84")
sp_utm <- spTransform(
  sp_pts,
  CRS("+proj=utm +zone=6 +south +datum=WGS84 +units=m +no_defs")
)

vgm_emp <- gstat::variogram(fine_fraction ~ 1, sp_utm)
vgm_fit <- fit.variogram(vgm_emp, model = vgm("Sph"))

map_corners <- data.frame(
  Long = c(map_xlim[1], map_xlim[2], map_xlim[1], map_xlim[2]),
  Lat = c(map_ylim[1], map_ylim[1], map_ylim[2], map_ylim[2])
)
coordinates(map_corners) <- ~ Long + Lat
proj4string(map_corners) <- CRS("+proj=longlat +datum=WGS84")
map_corners_utm <- spTransform(
  map_corners,
  CRS("+proj=utm +zone=6 +south +datum=WGS84 +units=m +no_defs")
)
bb <- bbox(map_corners_utm)
grid_utm <- expand.grid(
  x = seq(bb[1, 1], bb[1, 2], length.out = grid_n),
  y = seq(bb[2, 1], bb[2, 2], length.out = grid_n)
)
coordinates(grid_utm) <- ~ x + y
gridded(grid_utm) <- TRUE
proj4string(grid_utm) <- CRS("+proj=utm +zone=6 +south +datum=WGS84 +units=m +no_defs")

ok_k <- krige(fine_fraction ~ 1, sp_utm, grid_utm, model = vgm_fit)
ok_ll <- spTransform(ok_k, CRS("+proj=longlat +datum=WGS84"))
ok_plot <- as.data.frame(ok_ll) |>
  rename(Long = coords.x1, Lat = coords.x2)

island_union <- st_make_valid(st_union(moorea_poly))
ok_pts <- st_as_sf(ok_plot, coords = c("Long", "Lat"), crs = 4326, remove = FALSE)
sample_pts <- st_as_sf(analysis, coords = c("Long", "Lat"), crs = 4326)
ok_plot$on_land <- as.logical(
  st_intersects(ok_pts, island_union, sparse = FALSE)[, 1]
)
# Keep predictions near samples (~3 km; sf/s2 distances are metres)
ok_near <- lengths(st_is_within_distance(ok_pts, sample_pts, dist = 3000)) > 0
fine_max <- max(analysis$fine_fraction, na.rm = TRUE)
ok_plot$var1.pred <- pmin(pmax(ok_plot$var1.pred, 0), fine_max)
ok_water <- ok_plot |>
  filter(!on_land, ok_near)

p_fig1 <- add_water_basemap(ggplot())
p_fig1 <- p_fig1 +
  geom_tile(
    data = ok_water,
    aes(x = Long, y = Lat, fill = var1.pred),
    alpha = 0.88
  )
p_fig1 <- add_island_overlay(p_fig1) +
  geom_point(
    data = analysis,
    aes(x = Long, y = Lat, fill = fine_fraction, shape = bay_group),
    color = "grey15",
    size = 2.6,
    stroke = 0.45
  ) +
  scale_shape_manual(
    name = NULL,
    values = c(
      "In bay (Cook's / Opunohu)" = 21,
      "Out of bay" = 22
    )
  ) +
  scale_fill_gradient(
    name = "Fine Fraction",
    low = "#F7F1E8",
    high = "#8B5A2B",
    limits = c(0, fine_max),
    trans = "sqrt",
    oob = scales::squish
  ) +
  coord_sf(
    xlim = map_xlim,
    ylim = map_ylim,
    expand = FALSE,
    default_crs = sf::st_crs(4326),
    crs = sf::st_crs(4326)
  ) +
  labs(
    title = "Fine Fraction — ordinary kriging",
    subtitle = "Circles = in-bay (Cook's / Opunohu); squares = out of bay",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.background = element_rect(fill = col_ocean, color = NA),
    panel.grid = element_blank()
  )

ggsave(
  file.path(out_dir, "04_fig1_kriging_fine_fraction.png"),
  p_fig1,
  width = 9,
  height = 7,
  dpi = 200,
  bg = "white"
)

# 3️⃣ Fig 2 — in-bay vs out-of-bay (Welch t-test) ----
# =============================================================================
# Fine Fraction is right-skewed → test on log1p; also report raw means.

tt_log <- t.test(log1p_fine ~ bay_group, data = analysis, var.equal = FALSE)
tt_raw <- t.test(fine_fraction ~ bay_group, data = analysis, var.equal = FALSE)

group_summ <- analysis |>
  group_by(bay_group) |>
  summarise(
    n = n(),
    mean_fine = mean(fine_fraction),
    sd_fine = sd(fine_fraction),
    median_fine = median(fine_fraction),
    mean_log1p = mean(log1p_fine),
    .groups = "drop"
  )

ttest_out <- tibble(
  test = c("Welch_t_log1p_FineFraction", "Welch_t_raw_FineFraction"),
  n_in_bay = sum(analysis$bay %in% in_bay_prefixes),
  n_out_bay = sum(!(analysis$bay %in% in_bay_prefixes)),
  mean_in_bay_raw = group_summ$mean_fine[group_summ$bay_group == "In bay (Cook's / Opunohu)"],
  mean_out_bay_raw = group_summ$mean_fine[group_summ$bay_group == "Out of bay"],
  estimate_diff = c(unname(tt_log$estimate[1] - tt_log$estimate[2]),
                    unname(tt_raw$estimate[1] - tt_raw$estimate[2])),
  t = c(unname(tt_log$statistic), unname(tt_raw$statistic)),
  df = c(unname(tt_log$parameter), unname(tt_raw$parameter)),
  p_value = c(tt_log$p.value, tt_raw$p.value),
  note = c(
    "Primary test (skewed response: log1p)",
    "Same contrast on raw Fine Fraction (supplementary)"
  )
)
write_csv(ttest_out, file.path(out_dir, "04_fig2_ttest.csv"))

p_fig2 <- ggplot(
  analysis,
  aes(x = bay_group, y = fine_fraction, fill = bay_group)
) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, height = 0, size = 2.2, alpha = 0.85, shape = 21) +
  scale_fill_manual(values = c(col_in, col_out), guide = "none") +
  labs(
    title = "Fine Fraction — in bay vs out of bay",
    subtitle = paste0(
      "Welch t-test on log1p(Fine Fraction): p = ",
      signif(tt_log$p.value, 3),
      "  |  in-bay n = ",
      sum(analysis$bay %in% in_bay_prefixes),
      ", out-of-bay n = ",
      sum(!(analysis$bay %in% in_bay_prefixes))
    ),
    x = NULL,
    y = "Fine Fraction"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.x = element_blank())

ggsave(
  file.path(out_dir, "04_fig2_in_vs_out_bay.png"),
  p_fig2,
  width = 7,
  height = 5.5,
  dpi = 200,
  bg = "white"
)

# 4️⃣ Fig 3 — in-bay linear model (outlet only; no grain size) ----
# =============================================================================
# Default: outlet only. Set fig3_include_shore <- TRUE in Settings to also
# fit/plot the outlet + shore candidate if PIs ask later.

in_bay <- analysis |>
  filter(bay %in% in_bay_prefixes)

m_outlet <- lm(log1p_fine ~ dist_outlet_m, data = in_bay)
m_outlet_shore <- if (isTRUE(fig3_include_shore)) {
  lm(log1p_fine ~ dist_outlet_m + dist_shore_m, data = in_bay)
} else {
  NULL
}

lm_row <- function(fit, label) {
  s <- summary(fit)
  co <- s$coefficients
  tibble(
    model = label,
    n = nobs(fit),
    aic = AIC(fit),
    adj_r2 = s$adj.r.squared,
    term = rownames(co),
    estimate = co[, 1],
    se = co[, 2],
    t = co[, 3],
    p_value = co[, 4]
  )
}

lm_results <- lm_row(m_outlet, "log1p(Fine Fraction) ~ outlet")
if (!is.null(m_outlet_shore)) {
  lm_results <- bind_rows(
    lm_results,
    lm_row(m_outlet_shore, "log1p(Fine Fraction) ~ outlet + shore")
  )
}
write_csv(lm_results, file.path(out_dir, "04_fig3_lm_results.csv"))

# Main-text Fig 3: outlet only, on log1p scale to match the LM
p_fig3 <- ggplot(in_bay, aes(x = dist_outlet_m, y = log1p_fine)) +
  geom_point(aes(color = bay), size = 2.8, alpha = 0.9) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "grey20",
    linewidth = 0.7
  ) +
  scale_color_manual(values = c("C" = "#1B4F72", "O" = "#117A65"), name = "Bay") +
  labs(
    title = "In-bay only — Fine Fraction vs outlet distance",
    subtitle = paste0(
      "LM: log1p(Fine Fraction) ~ outlet  |  ",
      "adj R² = ",
      signif(summary(m_outlet)$adj.r.squared, 3),
      ", outlet p = ",
      signif(summary(m_outlet)$coefficients["dist_outlet_m", 4], 3),
      "  |  n = ",
      nrow(in_bay)
    ),
    x = "Distance from nearest major outlet (m)",
    y = "log1p(Fine Fraction)"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(out_dir, "04_fig3_inbay_lm.png"),
  p_fig3,
  width = 7.5,
  height = 5.5,
  dpi = 200,
  bg = "white"
)

# Optional shore model plot (off unless Settings flag is TRUE)
if (!is.null(m_outlet_shore)) {
  p_fig3_shore <- ggplot(in_bay, aes(x = dist_outlet_m, y = log1p_fine)) +
    geom_point(aes(color = bay, size = dist_shore_m), alpha = 0.9) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      color = "grey20",
      linewidth = 0.7
    ) +
    scale_color_manual(values = c("C" = "#1B4F72", "O" = "#117A65"), name = "Bay") +
    scale_size_continuous(name = "Shore distance (m)", range = c(2, 7)) +
    labs(
      title = "In-bay only — Fine Fraction vs outlet (+ shore in model)",
      subtitle = paste0(
        "LM: log1p(Fine Fraction) ~ outlet + shore  |  ",
        "adj R² = ",
        signif(summary(m_outlet_shore)$adj.r.squared, 3),
        ", outlet p = ",
        signif(summary(m_outlet_shore)$coefficients["dist_outlet_m", 4], 3),
        ", shore p = ",
        signif(summary(m_outlet_shore)$coefficients["dist_shore_m", 4], 3),
        "  |  n = ",
        nrow(in_bay)
      ),
      x = "Distance from nearest major outlet (m)",
      y = "log1p(Fine Fraction)"
    ) +
    theme_bw(base_size = 12)

  ggsave(
    file.path(out_dir, "04_fig3_inbay_lm_outlet_shore.png"),
    p_fig3_shore,
    width = 8,
    height = 5.5,
    dpi = 200,
    bg = "white"
  )
}

# 5️⃣ Short summary for the paper path ----
# =============================================================================

syn <- c(
  paste0("# Paper figures summary (script 04)"),
  paste0(""),
  paste0("Run: ", stamp),
  paste0("Response: Fine Fraction"),
  paste0(
    "In-bay prefixes: ",
    paste(in_bay_prefixes, collapse = ", "),
    " (Cook's + Opunohu); all other sites = out of bay."
  ),
  paste0("n analysis-ready: ", nrow(analysis)),
  paste0(
    "n in-bay / out-of-bay: ",
    sum(analysis$bay %in% in_bay_prefixes),
    " / ",
    sum(!(analysis$bay %in% in_bay_prefixes))
  ),
  paste0(""),
  paste0("## Fig 1"),
  paste0("- Ordinary kriging of Fine Fraction: 04_fig1_kriging_fine_fraction.png"),
  paste0(""),
  paste0("## Fig 2"),
  paste0(
    "- Welch t-test on log1p(Fine Fraction), in vs out bay: p = ",
    signif(tt_log$p.value, 3)
  ),
  paste0(
    "- Mean Fine Fraction in-bay = ",
    signif(group_summ$mean_fine[1], 3),
    "; out-of-bay = ",
    signif(group_summ$mean_fine[2], 3)
  ),
  paste0("- Plot: 04_fig2_in_vs_out_bay.png; table: 04_fig2_ttest.csv"),
  paste0(""),
  paste0("## Fig 3 (in-bay only; no grain size)"),
  paste0("- Main model: log1p(Fine Fraction) ~ outlet (fig3_include_shore = ",
         fig3_include_shore, ")"),
  paste0(
    "- log1p(Fine Fraction) ~ outlet: adj R² = ",
    signif(summary(m_outlet)$adj.r.squared, 3),
    ", outlet p = ",
    signif(summary(m_outlet)$coefficients["dist_outlet_m", 4], 3),
    ", AIC = ",
    signif(AIC(m_outlet), 5)
  ),
  if (!is.null(m_outlet_shore)) {
    paste0(
      "- Also fit outlet + shore: adj R² = ",
      signif(summary(m_outlet_shore)$adj.r.squared, 3),
      ", outlet p = ",
      signif(summary(m_outlet_shore)$coefficients["dist_outlet_m", 4], 3),
      ", shore p = ",
      signif(summary(m_outlet_shore)$coefficients["dist_shore_m", 4], 3)
    )
  } else {
    "- Outlet + shore model not run (set fig3_include_shore <- TRUE to enable)."
  },
  paste0("- Plot: 04_fig3_inbay_lm.png"),
  paste0("- Coefficients: 04_fig3_lm_results.csv"),
  paste0(""),
  paste0("## Note"),
  paste0(
    "- Script 03 exploratory analyses are retained as archive context; ",
    "script 04 is the lean paper figure path."
  )
)
writeLines(syn, file.path(out_dir, "04_summary.md"))

meta <- c(
  paste0("RUN: ", stamp),
  paste0("Input: ", snapshot_file),
  paste0("n: ", nrow(analysis)),
  paste0("In-bay prefixes: ", paste(in_bay_prefixes, collapse = ",")),
  paste0("Welch t (log1p) p: ", signif(tt_log$p.value, 4)),
  paste0("Fig3 model: outlet only (fig3_include_shore=", fig3_include_shore, ")"),
  paste0("Fig3 adj R2: ", signif(summary(m_outlet)$adj.r.squared, 4)),
  paste0("Fig3 outlet p: ",
         signif(summary(m_outlet)$coefficients["dist_outlet_m", 4], 4)),
  "Outputs: outputs/04_*"
)
writeLines(meta, file.path(out_dir, "04_run_metadata.txt"))

message(stamp, " — paper figures written to outputs/04_*")
message(
  "Fig2 Welch t (log1p) p = ", signif(tt_log$p.value, 3),
  " | Fig3 outlet-only adj R2 = ",
  signif(summary(m_outlet)$adj.r.squared, 3),
  ", p = ",
  signif(summary(m_outlet)$coefficients["dist_outlet_m", 4], 3)
)
