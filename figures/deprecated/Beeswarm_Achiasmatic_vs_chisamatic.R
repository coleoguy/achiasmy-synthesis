# ggbeeswarm density-shaped swarms (varwidth) in 1 panel:
# 1) Mammalia: x = Subclass (Theria vs Metatheria), color = Meiosis type (NO legend)
# 2) Coleoptera: x = Suborder (Polyphaga left, Adephaga right), color = Meiosis type (NO legend)
# 3) Diptera: x = Suborder, color = Meiosis type (KEEP legend)
# Mean line per x group in each plot
# Add "n=" under each x-axis tick label (sample size per x group)

setwd("~/GitHub/achiasmy-synthesis/data")
set.seed(1)

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(ggbeeswarm)

# -------------------------
# Common meiosis recode
# -------------------------
recode_meiosis <- function(x) {
  x <- str_to_lower(str_squish(as.character(x)))
  dplyr::case_when(
    x %in% c("achiasmatic", "asynaptic", "distance pairing") ~ "Achiasmatic",
    is.na(x) | x == ""                                       ~ "Chiasmatic",
    TRUE                                                     ~ NA_character_
  )
}

# -------------------------
# Mammalia prep
# -------------------------
mamm <- read_csv("Mammalia achiasmy.csv", show_col_types = FALSE) %>%
  mutate(across(any_of(c("Class","Subclass","Order","Family","Genus","Species")),
                ~ str_squish(as.character(.x)))) %>%
  mutate(
    Meiosis_plot     = recode_meiosis(`Meiosis Type`),
    `Diploid Number` = suppressWarnings(as.numeric(`Diploid Number`)),
    Haploid          = `Diploid Number` / 2,
    organism_id      = paste(Class, Subclass, Order, Family, Genus, Species, sep="|")
  ) %>%
  filter(!is.na(Haploid), !is.na(Meiosis_plot),
         !is.na(Subclass), Subclass %in% c("Theria","Metatheria")) %>%
  group_by(organism_id) %>%
  summarise(
    Haploid      = sample(Haploid, size = 1),
    Subclass     = first(Subclass),
    Meiosis_plot = if (any(Meiosis_plot == "Achiasmatic")) "Achiasmatic" else "Chiasmatic",
    .groups = "drop"
  ) %>%
  mutate(
    x_group      = factor(Subclass, levels = c("Metatheria", "Theria")),
    Meiosis_plot = factor(Meiosis_plot, levels = c("Chiasmatic", "Achiasmatic"))
  )

mamm_mean <- mamm %>%
  group_by(x_group) %>%
  summarise(mean_hap = mean(Haploid, na.rm = TRUE), .groups = "drop")

mamm_n <- mamm %>% count(x_group, name = "n")

# -------------------------
# Insect prep (Coleoptera/Diptera)
# -------------------------
prep_insect <- function(file, keep_suborders = NULL, suborder_levels = NULL) {
  dat <- read_csv(file, show_col_types = FALSE) %>%
    mutate(across(any_of(c("Class","Order","Suborder","Family","Genus","Species")),
                  ~ str_squish(as.character(.x)))) %>%
    mutate(
      Meiosis_plot     = recode_meiosis(`Meiosis Type`),
      `Diploid Number` = suppressWarnings(as.numeric(`Diploid Number`)),
      Haploid          = `Diploid Number` / 2,
      Suborder_clean   = str_to_lower(str_squish(as.character(Suborder)))
    ) %>%
    filter(!is.na(Haploid), Haploid > 0,
           !is.na(Suborder_clean), Suborder_clean != "",
           !is.na(Meiosis_plot))
  
  if (!is.null(keep_suborders)) {
    keep_low <- str_to_lower(keep_suborders)
    dat <- dat %>% filter(Suborder_clean %in% keep_low)
  }
  
  dat <- dat %>%
    mutate(Suborder = case_when(
      Suborder_clean == "adephaga"  ~ "Adephaga",
      Suborder_clean == "polyphaga" ~ "Polyphaga",
      TRUE                          ~ Suborder
    ))
  
  id_cols <- intersect(c("Class","Order","Suborder","Family","Genus","Species"), names(dat))
  
  out <- dat %>%
    mutate(organism_id = do.call(paste, c(across(all_of(id_cols)), sep="|"))) %>%
    group_by(organism_id) %>%
    summarise(
      Haploid      = sample(Haploid, size = 1),
      Suborder     = first(Suborder),
      Meiosis_plot = if (any(Meiosis_plot == "Achiasmatic")) "Achiasmatic" else "Chiasmatic",
      .groups = "drop"
    ) %>%
    mutate(
      Meiosis_plot = factor(Meiosis_plot, levels = c("Chiasmatic", "Achiasmatic")),
      x_group = if (is.null(suborder_levels)) {
        factor(Suborder, levels = sort(unique(Suborder)))
      } else {
        factor(Suborder, levels = suborder_levels)
      }
    )
  
  out
}

# Coleoptera: ONLY Adephaga + Polyphaga, plotted as Polyphaga then Adephaga
cole <- prep_insect(
  "Coleoptera achiasmy.csv",
  keep_suborders  = c("adephaga", "polyphaga"),
  suborder_levels = c("Polyphaga", "Adephaga")
)

# Diptera: all suborders
dipt <- prep_insect("Diptera achiasmy.csv")

cole_mean <- cole %>% group_by(x_group) %>% summarise(mean_hap = mean(Haploid, na.rm = TRUE), .groups="drop")
dipt_mean <- dipt %>% group_by(x_group) %>% summarise(mean_hap = mean(Haploid, na.rm = TRUE), .groups="drop")

cole_n <- cole %>% count(x_group, name = "n")
dipt_n <- dipt %>% count(x_group, name = "n")

# -------------------------
# Plot helper (no overlap across x groups + add n= to tick labels)
# -------------------------
make_plot <- function(df, mean_df, n_df, title_text, xlab, x_angle = 45, show_legend = TRUE) {
  
  # build per-plot tick labels with n
  lev <- levels(df$x_group)
  n_map <- setNames(n_df$n, as.character(n_df$x_group))
  tick_labels <- setNames(
    paste0(lev, "\n(n=", n_map[lev], ")"),
    lev
  )
  
  ggplot(df, aes(x = x_group, y = Haploid, color = Meiosis_plot)) +
    ggbeeswarm::geom_quasirandom(
      method    = "quasirandom",
      groupOnX  = TRUE,
      varwidth  = TRUE,
      width     = 0.45,
      alpha     = 0.30,
      size      = 1.8
    ) +
    geom_crossbar(
      data = mean_df,
      aes(x = x_group, y = mean_hap, ymin = mean_hap, ymax = mean_hap),
      inherit.aes = FALSE,
      width = 0.85
    ) +
    scale_x_discrete(labels = tick_labels) +
    labs(
      x = xlab,
      y = "Haploid chromosome number",
      color = "Meiosis type",
      title = title_text
    ) +
    theme_classic() +
    theme(
      axis.text.x  = element_text(angle = x_angle, hjust = 1),
      axis.title.x = element_text(margin = margin(t = 10)),
      plot.margin  = margin(t = 6, r = 6, b = 18, l = 6),
      legend.position = if (show_legend) "right" else "none"
    )
}

# -------------------------
# Build plots
# -------------------------
p_mamm <- make_plot(mamm, mamm_mean, mamm_n, "Mammalia", "Subclass", x_angle = 45, show_legend = FALSE)
p_cole <- make_plot(cole, cole_mean, cole_n, "Coleoptera", "Suborder", x_angle = 45, show_legend = FALSE)
p_dipt <- make_plot(dipt, dipt_mean, dipt_n, "Diptera", "Suborder", x_angle = 45, show_legend = TRUE)

# -------------------------
# Combine into 1 panel (3 columns)
# -------------------------
panel <- (p_mamm | p_cole | p_dipt)

# Optional save:

# =========================
# SAVE AS EMF
# =========================

devEMF::emf(
  file = "Beeswarm_mamm_cole_dipt.emf",
  width = 11, height = 8,
  bg = "white",
  coordDPI = 300
)
print(panel)
dev.off()

# Also display in the R plotting window
print(panel)
