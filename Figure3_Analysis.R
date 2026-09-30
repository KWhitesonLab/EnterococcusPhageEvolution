## ============================================================================
## Figure 3 -- Host range comparison between ancestral and evolved phages
## against the three training strains (DP05, DP12, DP15)
##
## Accelerated evolution of Enterococcus bacteriophage cocktails
## Dunham et al.
##
## This script reproduces all three panels of Figure 3 from a single input
## file, Figure3_RawData.xlsx:
##   Panel a) Ancestral phages/cocktails vs. their three training strains
##   Panel b) Evolved phages/cocktails, Parallel-host condition
##   Panel c) Evolved phages/cocktails, Mixed-host condition
##
## HOW TO USE
##   Put this script in the same folder as Figure3_RawData.xlsx and run it
##   -- no setwd() needed, from any working directory (see get_script_dir()
##   below). It writes six files into that same folder: Figure3_panel_a/b/c
##   as both .svg (vector, for Illustrator) and _preview.png (raster).
##
## INPUT DATA (Figure3_RawData.xlsx)
##   Evolved_Index     -- the 24 directed-evolution conditions (12 single
##                         phage/cocktail isolates + 12 twenty-cycle-evolved
##                         cocktails), each with its evolution condition
##                         (Parallel/Mixed), lineage (Similar/Dissimilar),
##                         and biological replicate number.
##   Ancestral_Index    -- the 8 ancestral phage/cocktail conditions (5
##                         single phages + a Similar cocktail + a
##                         Dissimilar cocktail, with the shared phage "Bob"
##                         appearing under both lineage tags).
##   Plate_Layout        -- the well-to-host map (which of the 96 wells on
##                         every plate contains which of the three training
##                         strains, and which single well per host is the
##                         untreated growth control) -- one shared layout,
##   One sheet per condition (named by its code, e.g. "C1M", "S2I_T20",
##   "Anc_Bill") -- the raw LogPhase600 plate-reader export for that
##   condition, copied in verbatim: instrument metadata, the "Time" header
##   row, and 289 timepoint readings (0-48h, every 10 minutes) across all
##   96 wells.
## ============================================================================

library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(svglite)
library(patchwork)
library(stringr)

## ----------------------------------------------------------------------
## 0. Locate this script's own folder, so the data file and outputs are
##    found/written relative to where the script lives rather than
##    getwd() -- this makes the script runnable via Rscript, source(), or
##    RStudio's Run/Source from any working directory, as long as
##    Figure3_RawData.xlsx sits next to it.
## ----------------------------------------------------------------------

get_script_dir <- function() {
  # Rscript (command-line) invocation: --file= is in the raw command args.
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  # source("this_script.R"): the sourcing frame records its own file path.
  for (fr in sys.frames()) {
    if (!is.null(fr$ofile)) return(dirname(normalizePath(fr$ofile)))
  }
  # RStudio, run interactively: ask the IDE for the active document's path.
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) return(dirname(normalizePath(ctx$path)))
  }
  # Last resort: current working directory, with a clear warning so a
  # wrong-folder failure never happens silently.
  warning(
    "Could not auto-detect this script's folder (unusual execution method) -- ",
    "falling back to the current working directory (", getwd(), "). ",
    "If Figure3_RawData.xlsx isn't found below, either run this script with ",
    "Rscript, use source(), or setwd() to the folder containing both files first."
  )
  getwd()
}

script_dir <- get_script_dir()
data_file  <- file.path(script_dir, "Figure3_RawData.xlsx")
out_dir    <- script_dir
stopifnot(file.exists(data_file))

time_vec <- seq(0, 48, by = 1/6) %>% round(digits = 3)   # 289 points, every 10 min over 48h

## ----------------------------------------------------------------------
## 1. Read the plate layout and the two condition indexes.
## ----------------------------------------------------------------------

plate_layout <- read.xlsx(data_file, sheet = "Plate_Layout", skipEmptyRows = FALSE)
plate_layout <- plate_layout[, c("well", "host", "lat", "lon", "is_control")]
plate_layout$is_control <- as.logical(plate_layout$is_control)

registry           <- read.xlsx(data_file, sheet = "Evolved_Index")
ancestral_registry <- read.xlsx(data_file, sheet = "Ancestral_Index")

## ----------------------------------------------------------------------
## 2. Read one condition's raw sheet and reshape it to long format
##    (time, well, OD600).
##
##    Each raw sheet is the instrument export copied in verbatim: a
##    metadata preamble of variable length, then the "Time" header row,
##    then 289 timepoint rows. The header row is located by searching
##    column A for the literal value "Time" rather than assuming a fixed
##    number of preamble lines, since that length isn't consistent across
##    export files. The instrument's own Time column (HH:MM:SS text) is
##    discarded and replaced by position with the precomputed time_vec,
##    since row order is what actually encodes elapsed time.
## ----------------------------------------------------------------------

read_condition_sheet <- function(sheet_name) {
  raw <- read.xlsx(data_file, sheet = sheet_name, colNames = FALSE, skipEmptyRows = FALSE)
  hdr_row <- which(raw[[1]] == "Time")[1]
  stopifnot(!is.na(hdr_row))
  well_cols <- as.character(unlist(raw[hdr_row, -1]))
  data_rows <- raw[(hdr_row + 1):(hdr_row + length(time_vec)), -1, drop = FALSE]
  colnames(data_rows) <- well_cols
  data_rows[] <- lapply(data_rows, function(x) suppressWarnings(as.numeric(x)))
  data_rows$time <- time_vec
  pivot_longer(data_rows, -time, names_to = "well", values_to = "OD600")
}

## ----------------------------------------------------------------------
## 3. Per-condition loader: merge the raw long-format readings with the
##    plate layout, keep only the three training-strain hosts, and
##    compute the technical-replicate mean/SD at each timepoint. Wells
##    are grouped by (host, is_control, time): is_control separates each
##    host's single untreated growth-control well from the ~10 wells
##    receiving the actual phage/cocktail treatment, so this reproduces
##    exactly the control-vs-treatment split the figure plots.
## ----------------------------------------------------------------------

load_condition <- function(sheet_name) {
  raw_long <- read_condition_sheet(sheet_name)
  merge(plate_layout, raw_long, by = "well") %>%
    filter(host %in% c("DP05", "DP12", "DP15")) %>%
    group_by(host, is_control, time) %>%
    mutate(OD600.avg = mean(OD600), OD600.sd = sd(OD600)) %>%
    ungroup()
}

## ----------------------------------------------------------------------
## 4. Build the long-format dataset behind panels b) and c): every
##    evolved condition's per-well readings, tagged with its evolution
##    condition, lineage, sub-row (isolate vs. cocktail), and biological
##    replicate from Evolved_Index.
## ----------------------------------------------------------------------

all_data <- bind_rows(lapply(seq_len(nrow(registry)), function(i) {
  reg_row <- registry[i, ]
  merged  <- load_condition(reg_row$sheet)
  merged$cond_code <- reg_row$sheet
  merged$sub_row   <- reg_row$sub_row
  merged$evolution <- reg_row$evolution
  merged$lineage   <- reg_row$lineage
  merged$replicate <- reg_row$replicate
  merged
}))

host_labels <- c(
  DP05 = 'italic("E. faecalis")~"DP05"',
  DP12 = 'italic("E. faecalis")~"DP12"',
  DP15 = 'italic("E. faecalis")~"DP15"'
)

all_data <- all_data %>%
  mutate(
    series  = if_else(is_control, "Control", paste0("Exp. ", replicate)),
    host    = factor(host, levels = c("DP05", "DP12", "DP15"), labels = host_labels[c("DP05","DP12","DP15")]),
    lineage = factor(lineage, levels = c("Similar", "Dissimilar")),
    sub_row = factor(sub_row, levels = c("Cocktail", "Isolate"), labels = c("Cocktails", "Isolates"))
  )

series_colors <- c(
  "Exp. 1"  = "#fc8d62",
  "Exp. 2"  = "#1eb735",
  "Exp. 3"  = "#8da0cb",
  "Control" = "#d1d1d1"
)

## ----------------------------------------------------------------------
## 5. Shared plot theme -- square panels, no gridlines, black text,
##    thicker mean lines over fainter individual-well lines.
## ----------------------------------------------------------------------

fig3_theme <- theme_bw(base_size = 7) +
  theme(
    aspect.ratio      = 1,
    panel.grid        = element_blank(),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.35),
    axis.line         = element_blank(),
    strip.background  = element_blank(),
    strip.text        = element_blank(),
    panel.spacing     = unit(0.12, "cm"),
    axis.text         = element_text(colour = "black", size = 6),
    axis.title        = element_text(colour = "black", size = 7, face = "bold"),
    plot.title        = element_blank(),
    legend.title      = element_text(colour = "black", size = 7, face = "bold"),
    legend.text       = element_text(colour = "black", size = 6.5),
    legend.key.size   = unit(0.3, "cm"),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    text              = element_text(colour = "black")
  )

IND_ALPHA <- 0.3
IND_LWD   <- 0.2
MEAN_LWD  <- 1.15

## ----------------------------------------------------------------------
## 6. Panels b) and c): one sub-plot per host x lineage combination, with
##    Cocktails stacked above Isolates.
## ----------------------------------------------------------------------

build_subrow_plot <- function(data, evolution_filter, subrow_filter, legend_name) {
  d <- filter(data, evolution == evolution_filter, sub_row == subrow_filter)
  d_mean <- distinct(d, sub_row, host, lineage, cond_code, series, time, OD600.avg)

  ggplot() +
    geom_line(data = d, aes(x = time, y = OD600, group = interaction(cond_code, well), color = series),
              linewidth = IND_LWD, alpha = IND_ALPHA) +
    geom_line(data = d_mean, aes(x = time, y = OD600.avg, group = interaction(cond_code, series), color = series),
              linewidth = MEAN_LWD) +
    facet_grid(sub_row ~ host + lineage, labeller = label_parsed) +
    scale_color_manual(values = series_colors, breaks = names(series_colors), name = legend_name) +
    scale_x_continuous(breaks = seq(0, 50, 10)) +
    coord_cartesian(ylim = c(-0.1, 2.0), xlim = c(0, 48)) +
    labs(x = "Time (h)", y = expression(OD[600])) +
    fig3_theme
}

build_panel_bc <- function(evolution_filter) {
  p_top <- build_subrow_plot(all_data, evolution_filter, "Cocktails", "Cocktails") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p_bot <- build_subrow_plot(all_data, evolution_filter, "Isolates", "Isolates")
  (p_top / p_bot) + plot_layout(heights = c(1, 1))
}

p_b <- build_panel_bc("Parallel")
p_c <- build_panel_bc("Mixed")

## ----------------------------------------------------------------------
## 7. Panel a): ancestral phages and cocktails. phage_code comes straight
##    from Ancestral_Index ("Bill" / "Car" / "SDS1" / "CCS4" / "Bob" /
##    "Cocktail"), which is also what isolate_colors/cocktail_colors below
##    are keyed on.
## ----------------------------------------------------------------------

ancestral_data <- bind_rows(lapply(seq_len(nrow(ancestral_registry)), function(i) {
  reg_row <- ancestral_registry[i, ]
  merged  <- load_condition(reg_row$sheet)
  merged$phage_code <- reg_row$phage_code
  merged$sub_row     <- reg_row$sub_row
  merged$lineage      <- reg_row$lineage
  merged$series      <- if_else(merged$is_control, "Control", reg_row$phage_code)
  merged
})) %>%
  mutate(
    host    = factor(host, levels = c("DP05", "DP12", "DP15"), labels = host_labels[c("DP05","DP12","DP15")]),
    lineage = factor(lineage, levels = c("Similar", "Dissimilar")),
    sub_row = factor(sub_row, levels = c("Cocktail", "Isolate"), labels = c("Cocktails", "Isolates"))
  )

cocktail_colors <- c("Cocktail" = "#000000", "Control" = "#d1d1d1")
isolate_colors  <- c(
  "Bill" = "#c400c4", "Bob" = "#8ea0e0", "Car" = "#fc8d62",
  "CCS4" = "#b5121b", "SDS1" = "#3fae49", "Control" = "#d1d1d1"
)

build_ancestral_subplot <- function(data, subrow_filter, colors, legend_name) {
  d <- filter(data, sub_row == subrow_filter)
  d_mean <- distinct(d, sub_row, host, lineage, phage_code, series, time, OD600.avg)

  ggplot() +
    geom_line(data = d, aes(x = time, y = OD600, group = interaction(phage_code, well), color = series),
              linewidth = IND_LWD, alpha = IND_ALPHA) +
    geom_line(data = d_mean, aes(x = time, y = OD600.avg, group = interaction(phage_code, series), color = series),
              linewidth = MEAN_LWD) +
    facet_grid(sub_row ~ host + lineage, labeller = label_parsed) +
    scale_color_manual(values = colors, breaks = names(colors), name = legend_name) +
    scale_x_continuous(breaks = seq(0, 50, 10)) +
    coord_cartesian(ylim = c(-0.1, 2.0), xlim = c(0, 48)) +
    labs(x = "Time (h)", y = expression(OD[600])) +
    fig3_theme
}

p_a_top <- build_ancestral_subplot(ancestral_data, "Cocktails", cocktail_colors, "Cocktails") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
p_a_bot <- build_ancestral_subplot(ancestral_data, "Isolates", isolate_colors, "Isolates")

p_a <- (p_a_top / p_a_bot) + plot_layout(heights = c(1, 1))

## Sanity check: ancestral phage Bill should show little to no growth
## suppression against DP15, unlike DP05/DP12.
sanity_check <- ancestral_data %>%
  filter(phage_code %in% c("Bill", "Control"), sub_row == "Isolates", lineage == "Similar", time == 48) %>%
  distinct(host, series, OD600.avg) %>%
  arrange(host, series)
cat("\n=== Panel a) sanity check (Bill, t=48h) ===\n"); print(sanity_check); cat("\n")

## ----------------------------------------------------------------------
## 8. Export panels a), b), c) as separate files, 7.25in wide (Science
##    Advances full-text-width figure), as both editable SVG and a PNG
##    preview.
## ----------------------------------------------------------------------

fig_width <- 7.25

save_panel <- function(p, name, height) {
  ggsave(file.path(out_dir, paste0(name, ".svg")), plot = p, width = fig_width, height = height, units = "in", device = svglite)
  ggsave(file.path(out_dir, paste0(name, "_preview.png")), plot = p, width = fig_width, height = height, units = "in", dpi = 300)
}

save_panel(p_a, "Figure3_panel_a", 2.85)
save_panel(p_b, "Figure3_panel_b", 2.85)
save_panel(p_c, "Figure3_panel_c", 2.85)

cat("Done.\n")
cat("Panels b/c conditions loaded:", paste(unique(all_data$cond_code), collapse = ", "), "\n")
cat("Panel a conditions loaded:", paste(unique(ancestral_data$phage_code), collapse = ", "), "\n")
