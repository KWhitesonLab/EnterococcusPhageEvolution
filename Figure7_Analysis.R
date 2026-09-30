## ============================================================================
## Figure 7 -- Host range heatmap, 72 naive E. faecalis strains + training
## strains/media controls, ancestral and evolved phages/cocktails
##
## Accelerated evolution of Enterococcus bacteriophage cocktails
## Dunham et al.
##
## Reproduces Figure 7 (heatmap panel + Summary bar panel) from Figure7_RawData.xlsx.
##
## HOW TO USE
##   Put this script in the same folder as Figure7_RawData.xlsx and run it
##   -- no setwd() needed, from any working directory (see get_script_dir()
##   below). It writes Figure7_heatmap.svg (vector, for Illustrator) and
##   Figure7_heatmap_preview.png into that same folder.
##
## INPUT DATA (Figure7_RawData.xlsx)
##   HRScoring_Naive  -- 72 naive E. faecalis strains x 35 phage/cocktail
##                        conditions, scored 1 = Resistant, 2 = Partial
##                        Suppression, 3 = Full Suppression.
##   HRScoring_TSBHI   -- the same 35 conditions scored against the 8
##                        training-strain/media-control rows (BHI1-5,
##                        DP05, DP12, DP15).
##   Strain_Order      -- the 72 naive strains in the exact order used on
##                        the figure's x-axis.
##   Condition_Index   -- the 35 condition codes with their short/full
##                        display labels and which of the four category
##                        blocks (Ancestral Phages / Ancestral Cocktails /
##                        Evolved Isolates / Evolved Cocktails) each
##                        belongs to.
## ============================================================================

library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(svglite)
library(patchwork)

## ----------------------------------------------------------------------
## 0. Locate this script's own folder, so the data file and outputs are
##    found/written relative to where the script lives rather than
##    getwd() -- this makes the script runnable via Rscript, source(), or
##    RStudio's Run/Source from any working directory, as long as
##    Figure7_RawData.xlsx sits next to it.
## ----------------------------------------------------------------------

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  for (fr in sys.frames()) {
    if (!is.null(fr$ofile)) return(dirname(normalizePath(fr$ofile)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) return(dirname(normalizePath(ctx$path)))
  }
  warning(
    "Could not auto-detect this script's folder (unusual execution method) -- ",
    "falling back to the current working directory (", getwd(), "). ",
    "If Figure7_RawData.xlsx isn't found below, either run this script with ",
    "Rscript, use source(), or setwd() to the folder containing both files first."
  )
  getwd()
}

script_dir <- get_script_dir()
data_file  <- file.path(script_dir, "Figure7_RawData.xlsx")
out_dir    <- script_dir
stopifnot(file.exists(data_file))

## ----------------------------------------------------------------------
## 1. Load scoring data and condition metadata.
## ----------------------------------------------------------------------

hr_naive_raw <- read.xlsx(data_file, sheet = "HRScoring_Naive")
hr_tsbhi_raw <- read.xlsx(data_file, sheet = "HRScoring_TSBHI")
strain_order_naive <- read.xlsx(data_file, sheet = "Strain_Order")$strain
condition_index <- read.xlsx(data_file, sheet = "Condition_Index")

rownames(hr_naive_raw) <- hr_naive_raw$strain
hr_naive <- hr_naive_raw[, -1, drop = FALSE]
rownames(hr_tsbhi_raw) <- hr_tsbhi_raw$strain
hr_tsbhi <- hr_tsbhi_raw[, -1, drop = FALSE]

strain_order <- c("BHI1", "BHI2", "BHI3", "BHI4", "BHI5", "DP05", "DP12", "DP15", strain_order_naive)

hr_all <- rbind(hr_tsbhi, hr_naive)
hr_all <- hr_all[strain_order, ]

## ----------------------------------------------------------------------
## 2. Condition order and display metadata, from Condition_Index (kept in
##    the workbook's own row order -- Similar-before-Dissimilar throughout).
## ----------------------------------------------------------------------

condition_index <- condition_index[order(condition_index$order_index), ]
cond_order        <- condition_index$condition_code
cond_labels       <- setNames(condition_index$full_label, cond_order)
cond_labels_short <- setNames(condition_index$short_label, cond_order)
cond_block        <- setNames(condition_index$block, cond_order)

stopifnot(setequal(cond_order, colnames(hr_all)))

strain_group <- setNames(
  ifelse(grepl("^BHI", strain_order), "Media Control",
         ifelse(strain_order %in% c("DP05", "DP12", "DP15"), "Training Strain", "Naive Strain")),
  strain_order
)

## ----------------------------------------------------------------------
## 3. Manual row (condition) y-positions: 3-tier gap hierarchy -- big gaps
##    between primary categories (Ancestral Phages / Ancestral Cocktails /
##    Evolved Isolates / Evolved Cocktails), small gaps between secondary
##    categories (Similar / Dissimilar), no gap between tertiary
##    categories (Mixed / Parallel). Done with an explicit numeric y scale
##    since facet_grid only supports one uniform panel.spacing value.
## ----------------------------------------------------------------------

n_cond <- length(cond_order)
big_idx   <- c(6, 12, 24)
small_idx <- c(9, 18, 30)
step      <- 1.0
small_gap <- 1.5
big_gap   <- 2.6

y_raw <- numeric(n_cond)
for (i in 2:n_cond) {
  inc <- if (i %in% big_idx) big_gap else if (i %in% small_idx) small_gap else step
  y_raw[i] <- y_raw[i - 1] + inc
}
y_final <- max(y_raw) + min(y_raw) - y_raw   # flip so row 1 (CCS4_C) lands at the top
cond_y <- setNames(y_final, cond_order)
y_limits <- range(y_final) + c(-step / 2 - 0.02, step / 2 + 0.02)

## ----------------------------------------------------------------------
## 4. Long-format heatmap data.
## ----------------------------------------------------------------------

hr_ordered <- hr_all[, cond_order]
mat <- as.matrix(hr_ordered)

df <- as.data.frame(mat) %>%
  mutate(strain = rownames(mat)) %>%
  pivot_longer(-strain, names_to = "cond", values_to = "score") %>%
  mutate(
    cond_y = cond_y[cond],
    strain = factor(strain, levels = strain_order),
    strain_group = factor(strain_group[as.character(strain)], levels = c("Media Control", "Training Strain", "Naive Strain")),
    score = factor(score, levels = c("1", "2", "3"))
  )

## ----------------------------------------------------------------------
## 5. Palette and column-label handling.
## ----------------------------------------------------------------------

score_colors <- c("1" = "#eff3ff", "2" = "#6baed6", "3" = "#0d2c39")
score_labels <- c("1" = "Resistant", "2" = "Partial Suppression", "3" = "Full Suppression")

## "ESKAPE E. faecalis" is far longer than the other strain IDs (which are
## all ~7-character codes); left at full length it forces the x-axis text
## row much taller than every other column needs, stealing panel height and
## breaking square-cell sizing. Abbreviated on the axis only.
strain_tick_labels <- setNames(strain_order, strain_order)
strain_tick_labels["ESKAPE E. faecalis"] <- "ESKAPE"

## ----------------------------------------------------------------------
## 6. Heatmap panel.
## ----------------------------------------------------------------------

p_heat <- ggplot(df, aes(x = strain, y = cond_y, fill = score)) +
  geom_tile(color = NA, height = step, width = 1) +
  scale_fill_manual(values = score_colors, labels = score_labels, name = "Phage Susceptibility", na.value = "white") +
  scale_x_discrete(expand = c(0, 0), position = "bottom", labels = strain_tick_labels) +
  scale_y_continuous(
    breaks = cond_y[cond_order], labels = cond_labels_short[cond_order],
    limits = y_limits, expand = c(0, 0)
  ) +
  facet_grid(. ~ strain_group, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3.6, family = "mono", lineheight = 0.7),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 4.2, family = "mono", lineheight = 0.7),
    axis.ticks.x = element_line(color = "grey50", linewidth = 0.1),
    axis.ticks.length.x = unit(1, "pt"),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 6.5),
    legend.text = element_text(size = 6),
    legend.key.size = unit(7, "pt"),
    legend.margin = margin(0, 0, 0, 2),
    plot.margin = margin(2, 2, 2, 2)
  )

## ----------------------------------------------------------------------
## 7. Summary panel: % of the 72 naive strains showing full suppression
##    and at-least-partial suppression, per condition -- computed over the
##    72 naive-strain columns only (training strains and media controls
##    excluded), per the manuscript's Host Range Analysis methods text.
## ----------------------------------------------------------------------

naive_mat <- hr_ordered[strain_order_naive, ]
stopifnot(nrow(naive_mat) == 72, setequal(colnames(naive_mat), cond_order))

summary_df <- data.frame(
  cond = cond_order,
  pct_full    = sapply(cond_order, function(cn) mean(naive_mat[[cn]] == "3") * 100),
  pct_partial = sapply(cond_order, function(cn) mean(naive_mat[[cn]] == "2") * 100)
) %>%
  mutate(
    cond_y = cond_y[cond],
    full_xmax = pct_full,
    partial_xmax = pct_full + pct_partial
  )

bar_half_height <- step / 2 * 0.8

p_summary <- ggplot(summary_df) +
  geom_rect(aes(xmin = 0, xmax = full_xmax, ymin = cond_y - bar_half_height, ymax = cond_y + bar_half_height),
            fill = score_colors[["3"]]) +
  geom_rect(aes(xmin = full_xmax, xmax = partial_xmax, ymin = cond_y - bar_half_height, ymax = cond_y + bar_half_height),
            fill = score_colors[["2"]]) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0), breaks = c(0, 50, 100)) +
  scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
  labs(x = "% of strains", y = NULL, title = "Summary") +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.title.x = element_text(size = 6),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.15),
    plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
    plot.margin = margin(2, 2, 2, 1)
  )

## ----------------------------------------------------------------------
## 8. Combine heatmap + summary via patchwork and export.
##    Target canvas: 7.25 x 3in (Science Advances full-text-width figure),
##    every heatmap cell square. The 5.5:1 width split between the two
##    panels is what makes the tiles render square at this canvas size;
##    it isn't a theme setting ggplot exposes directly.
## ----------------------------------------------------------------------

p_combined <- p_heat + p_summary + plot_layout(widths = c(5.5, 1))

fig_width  <- 7.25
fig_height <- 3.0

ggsave(file.path(out_dir, "Figure7_heatmap.svg"), p_combined, device = svglite, width = fig_width, height = fig_height, units = "in")
ggsave(file.path(out_dir, "Figure7_heatmap_preview.png"), p_combined, width = fig_width, height = fig_height, units = "in", dpi = 500, bg = "white")

cat("Done.\n")
cat("Conditions loaded:", length(cond_order), "\n")
cat("Naive strains:", length(strain_order_naive), " | Training/media rows:", nrow(hr_tsbhi), "\n")
