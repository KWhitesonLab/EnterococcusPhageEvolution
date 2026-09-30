## ============================================================================
## Figure 2A-B -- Cocktail suppression dynamics across directed evolution
## cycles (dilution-factor trajectory)
##
## Accelerated evolution of Enterococcus bacteriophage cocktails
## Dunham et al.
##
## Reproduces Figure 2A and Figure 2B from a single input file,
## Figure2AB_RawData.xlsx.
##
## HOW TO USE
##   Put this script in the same folder as Figure2AB_RawData.xlsx and run it
##   -- no setwd() needed, from any working directory (see get_script_dir()
##   below). It writes Figure2A.svg/.png (Mixed-host condition, single
##   panel) and Figure2B.svg/.png (Parallel-host condition, one facet per
##   host culture) into that same folder.
##
## INPUT DATA (Figure2AB_RawData.xlsx)
##   Collection_Results_Sister -- the Similar cocktail (Bill, Bob, Car):
##     one row per (experimental cycle, biological replicate, host
##     culture), recording the most dilute well of that cycle's serial
##     dilution series that still showed visible suppression.
##   Collection_Results_Cousin -- the same, for the Dissimilar cocktail
##     (Bob, CCS4, SDS1).
##   well.1 is a log-dilution index (0 = undiluted, -1 = 1:10, -2 = 1:100,
##   -3 = 1:1,000, -4 = 1:10,000); well.2 is a secondary well reading, not
##   used in this figure.
##   Figure 2A covers the Mixed-host condition (host == "Mix"), where all
##   three co-cultured strains of a host set were challenged together each
##   cycle. Figure 2B covers the Parallel-host condition -- the same three
##   host cultures (DP05, DP12, DP15) challenged separately.
## ============================================================================

library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(svglite)

## ----------------------------------------------------------------------
## 0. Locate this script's own folder, so the data file and outputs are
##    found/written relative to where the script lives rather than
##    getwd() -- this makes the script runnable via Rscript, source(), or
##    RStudio's Run/Source from any working directory, as long as
##    Figure2AB_RawData.xlsx sits next to it.
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
    "If Figure2AB_RawData.xlsx isn't found below, either run this script with ",
    "Rscript, use source(), or setwd() to the folder containing both files first."
  )
  getwd()
}

script_dir <- get_script_dir()
data_file  <- file.path(script_dir, "Figure2AB_RawData.xlsx")
out_dir    <- script_dir
stopifnot(file.exists(data_file))

## ----------------------------------------------------------------------
## 1. Load and combine the two cocktails' raw dilution-tracking data.
## ----------------------------------------------------------------------

sim <- read.xlsx(data_file, sheet = "Collection_Results_Sister")
dis <- read.xlsx(data_file, sheet = "Collection_Results_Cousin")

sim$cocktail <- "Similar"
dis$cocktail <- "Dissimilar"

raw <- bind_rows(sim, dis) %>%
  mutate(
    timepoint = as.numeric(timepoint),
    replicate = as.numeric(replicate),
    well.1    = as.numeric(well.1),
    host      = recode(host, "Mix" = "Mixed"),
    host      = factor(host, levels = c("Mixed", "DP05", "DP12", "DP15")),
    cocktail  = factor(cocktail, levels = c("Similar", "Dissimilar"))
  )

## ----------------------------------------------------------------------
## 2. Per-cycle, per-host mean and SD across the three biological
##    replicates. titer_proxy is the negated dilution index, so it
##    increases with suppression strength (higher = more dilute cocktail
##    still suppressed growth = more effective).
## ----------------------------------------------------------------------

summ <- raw %>%
  group_by(cocktail, host, timepoint) %>%
  summarise(
    mean_dilution = mean(well.1),
    sd_dilution   = sd(well.1),
    .groups = "drop"
  ) %>%
  mutate(
    titer_proxy = -mean_dilution,
    titer_lo    = pmax(titer_proxy - sd_dilution, 0),
    titer_hi    = titer_proxy + sd_dilution
  )

## ----------------------------------------------------------------------
## 3. Palette and axis labels, shared by both panels.
## ----------------------------------------------------------------------

cocktail_colors <- c("Similar" = "#0072B2", "Dissimilar" = "#D55E00")
cocktail_legend_labels <- c("Similar" = "Similar Phages", "Dissimilar" = "Dissimilar Phages")
dilution_labels <- c("1", "1:10", "1:100", "1:1,000", "1:10,000")

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, color = "grey85"),
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_blank()
  )

## ----------------------------------------------------------------------
## 4. Panel A: Mixed-host condition, single panel.
## ----------------------------------------------------------------------

summ_A <- summ %>% filter(host == "Mixed")

p_A <- ggplot(summ_A, aes(x = timepoint, y = titer_proxy, color = cocktail, fill = cocktail)) +
  geom_ribbon(aes(ymin = titer_lo, ymax = titer_hi), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.0) +
  scale_color_manual(values = cocktail_colors, labels = cocktail_legend_labels, name = NULL) +
  scale_fill_manual(values = cocktail_colors, guide = "none") +
  scale_x_log10(breaks = c(1, 10, 20)) +
  scale_y_continuous(breaks = 0:4, labels = dilution_labels, limits = c(0, 4.1)) +
  labs(
    title = "Mixed Hosts", subtitle = "DP05, DP12, & DP15",
    x = "Experimental Cycle", y = "Effective Cocktail Dilution"
  ) +
  base_theme

ggsave(file.path(out_dir, "Figure2A.svg"), p_A, device = svglite, width = 4.2, height = 4.3, units = "in")
ggsave(file.path(out_dir, "Figure2A_preview.png"), p_A, width = 4.2, height = 4.3, units = "in", dpi = 300, bg = "white")

## ----------------------------------------------------------------------
## 5. Panel B: Parallel-host condition, one facet per host culture.
## ----------------------------------------------------------------------

summ_B <- summ %>% filter(host %in% c("DP05", "DP12", "DP15")) %>% mutate(host = droplevels(host))

p_B <- ggplot(summ_B, aes(x = timepoint, y = titer_proxy, color = cocktail, fill = cocktail)) +
  geom_ribbon(aes(ymin = titer_lo, ymax = titer_hi), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.0) +
  facet_wrap(~host, nrow = 1, labeller = label_bquote(italic("E. faecalis") ~ .(as.character(host)))) +
  scale_color_manual(values = cocktail_colors, labels = cocktail_legend_labels, name = NULL) +
  scale_fill_manual(values = cocktail_colors, guide = "none") +
  scale_x_log10(breaks = c(1, 10, 20)) +
  scale_y_continuous(breaks = 0:4, labels = dilution_labels, limits = c(0, 4.1)) +
  labs(
    title = "Parallel Hosts",
    x = "Experimental Cycle", y = "Effective Cocktail Dilution"
  ) +
  base_theme

ggsave(file.path(out_dir, "Figure2B.svg"), p_B, device = svglite, width = 9.0, height = 4.3, units = "in")
ggsave(file.path(out_dir, "Figure2B_preview.png"), p_B, width = 9.0, height = 4.3, units = "in", dpi = 300, bg = "white")

cat("Done.\n")
