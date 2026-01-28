#SEGREGATION TYPE VIOLIN PLOT
#MEGHANN MCCONNELL
#BLACKMON LAB
#JAN 2026

setwd("~/GitHub/achiasmy-synthesis/data")
set.seed(1)

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(viridisLite)
library(scales)

# -------------------------
# viridis colors (3 categories)
# -------------------------
meiosis_levels <- c("Chiasmatic","Achiasmatic","Asynaptic")
meiosis_cols   <- setNames(viridisLite::viridis(3), meiosis_levels)

# -------------------------
# color helpers
# -------------------------
blend_hex <- function(c1, c2, w = 0.5) {
  r1 <- grDevices::col2rgb(c1)
  r2 <- grDevices::col2rgb(c2)
  r  <- (1 - w) * r1 + w * r2
  grDevices::rgb(r[1], r[2], r[3], maxColorValue = 255)
}
darken_hex <- function(col, f = 0.7) {
  r <- grDevices::col2rgb(col)
  r <- pmax(0, pmin(255, r * f))
  grDevices::rgb(r[1], r[2], r[3], maxColorValue = 255)
}

# shift ONLY Chiasmatic toward a blue hue
meiosis_cols["Chiasmatic"] <- blend_hex(meiosis_cols["Chiasmatic"], "#0072B2", w = 0.65)

# -------------------------
# Meiosis recode (3-category)
# -------------------------
recode_meiosis3 <- function(x) {
  x <- str_to_lower(str_squish(as.character(x)))
  dplyr::case_when(
    is.na(x) | x == ""                          ~ "Chiasmatic",
    x %in% c("asynaptic")                       ~ "Asynaptic",
    x %in% c("achiasmatic", "distance pairing") ~ "Achiasmatic",
    TRUE                                        ~ NA_character_
  )
}

achiasmatic_tokens <- c("achiasmatic", "distance pairing")
message("Raw 'Meiosis Type' values categorized as Achiasmatic: ",
        paste(shQuote(achiasmatic_tokens), collapse = ", "))

collapse_meiosis3 <- function(v) {
  v <- v[!is.na(v)]
  if (!length(v)) return(NA_character_)
  if (any(v == "Asynaptic"))   return("Asynaptic")
  if (any(v == "Achiasmatic")) return("Achiasmatic")
  if (any(v == "Chiasmatic"))  return("Chiasmatic")
  NA_character_
}

# -------------------------
# Stats for custom box
# -------------------------
calc_box_stats <- function(df) {
  df %>%
    group_by(x_group, Meiosis_plot) %>%
    summarise(
      ymin = min(Haploid, na.rm = TRUE),
      ymax = max(Haploid, na.rm = TRUE),
      q1   = as.numeric(quantile(Haploid, 0.25, na.rm = TRUE, type = 7)),
      med  = as.numeric(quantile(Haploid, 0.50, na.rm = TRUE, type = 7)),
      q3   = as.numeric(quantile(Haploid, 0.75, na.rm = TRUE, type = 7)),
      .groups = "drop"
    )
}

# -------------------------
# Violin + custom IQR box + median dot + min/max whiskers
# NOW includes violin_adjust (bandwidth multiplier) per plot
# -------------------------
make_violin_custombox <- function(df, n_df,
                                  title_text, xlab,
                                  x_angle = 45,
                                  show_legend = TRUE,
                                  overlay_dodge = FALSE,
                                  dodge_width = 0.85,
                                  violin_alpha = 0.80,
                                  violin_adjust = 1.0,
                                  box_width = 0.18,
                                  whisker_lwd = 0.6,
                                  dot_size = 2.0) {
  
  df <- df %>%
    mutate(
      Meiosis_plot = factor(as.character(Meiosis_plot), levels = meiosis_levels),
      x_group      = factor(x_group, levels = levels(x_group))
    )
  
  lev <- levels(df$x_group)
  n_map <- setNames(n_df$n, as.character(n_df$x_group))
  tick_labels <- setNames(paste0(lev, "\n(n=", n_map[lev], ")"), lev)
  
  stat_df <- calc_box_stats(df) %>%
    mutate(Meiosis_plot = factor(as.character(Meiosis_plot), levels = meiosis_levels))
  
  box_fill_map <- c(
    "Chiasmatic"  = darken_hex(meiosis_cols["Chiasmatic"],  f = 0.65),
    "Achiasmatic" = darken_hex(meiosis_cols["Achiasmatic"], f = 0.65),
    "Asynaptic"   = darken_hex(blend_hex(meiosis_cols["Asynaptic"], "#E69F00", w = 0.60), f = 0.75)
  )
  
  p <- ggplot(df, aes(x = x_group, y = Haploid))
  
  if (overlay_dodge) {
    p <- p +
      geom_violin(
        aes(fill = Meiosis_plot, color = Meiosis_plot,
            group = interaction(x_group, Meiosis_plot)),
        position = position_dodge(width = dodge_width),
        trim = FALSE,
        scale = "width",
        adjust = violin_adjust,
        alpha = violin_alpha,
        linewidth = 0.35
      ) +
      geom_linerange(
        data = stat_df,
        aes(x = x_group, ymin = ymin, ymax = ymax,
            group = interaction(x_group, Meiosis_plot)),
        position = position_dodge(width = dodge_width),
        inherit.aes = FALSE,
        linewidth = whisker_lwd,
        color = "black"
      ) +
      geom_rect(
        data = dplyr::filter(stat_df, Meiosis_plot == "Chiasmatic"),
        aes(
          xmin = as.numeric(x_group) - box_width/2,
          xmax = as.numeric(x_group) + box_width/2,
          ymin = q1, ymax = q3,
          group = interaction(x_group, Meiosis_plot)
        ),
        position = position_dodge(width = dodge_width),
        inherit.aes = FALSE,
        fill = box_fill_map["Chiasmatic"],
        color = NA,
        alpha = 1
      ) +
      geom_rect(
        data = dplyr::filter(stat_df, Meiosis_plot == "Achiasmatic"),
        aes(
          xmin = as.numeric(x_group) - box_width/2,
          xmax = as.numeric(x_group) + box_width/2,
          ymin = q1, ymax = q3,
          group = interaction(x_group, Meiosis_plot)
        ),
        position = position_dodge(width = dodge_width),
        inherit.aes = FALSE,
        fill = box_fill_map["Achiasmatic"],
        color = NA,
        alpha = 1
      ) +
      geom_rect(
        data = dplyr::filter(stat_df, Meiosis_plot == "Asynaptic"),
        aes(
          xmin = as.numeric(x_group) - box_width/2,
          xmax = as.numeric(x_group) + box_width/2,
          ymin = q1, ymax = q3,
          group = interaction(x_group, Meiosis_plot)
        ),
        position = position_dodge(width = dodge_width),
        inherit.aes = FALSE,
        fill = box_fill_map["Asynaptic"],
        color = NA,
        alpha = 1
      ) +
      geom_point(
        data = stat_df,
        aes(x = x_group, y = med,
            group = interaction(x_group, Meiosis_plot)),
        position = position_dodge(width = dodge_width),
        inherit.aes = FALSE,
        color = "white",
        size = dot_size
      )
  } else {
    p <- p +
      geom_violin(
        aes(fill = Meiosis_plot, color = Meiosis_plot),
        trim = FALSE,
        scale = "width",
        adjust = violin_adjust,
        alpha = violin_alpha,
        linewidth = 0.35
      ) +
      geom_linerange(
        data = stat_df,
        aes(x = x_group, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        linewidth = whisker_lwd,
        color = "black"
      ) +
      geom_rect(
        data = dplyr::filter(stat_df, Meiosis_plot == "Chiasmatic"),
        aes(
          xmin = as.numeric(x_group) - box_width/2,
          xmax = as.numeric(x_group) + box_width/2,
          ymin = q1, ymax = q3
        ),
        inherit.aes = FALSE,
        fill = box_fill_map["Chiasmatic"],
        color = NA,
        alpha = 1
      ) +
      geom_rect(
        data = dplyr::filter(stat_df, Meiosis_plot == "Achiasmatic"),
        aes(
          xmin = as.numeric(x_group) - box_width/2,
          xmax = as.numeric(x_group) + box_width/2,
          ymin = q1, ymax = q3
        ),
        inherit.aes = FALSE,
        fill = box_fill_map["Achiasmatic"],
        color = NA,
        alpha = 1
      ) +
      geom_rect(
        data = dplyr::filter(stat_df, Meiosis_plot == "Asynaptic"),
        aes(
          xmin = as.numeric(x_group) - box_width/2,
          xmax = as.numeric(x_group) + box_width/2,
          ymin = q1, ymax = q3
        ),
        inherit.aes = FALSE,
        fill = box_fill_map["Asynaptic"],
        color = NA,
        alpha = 1
      ) +
      geom_point(
        data = stat_df,
        aes(x = x_group, y = med),
        inherit.aes = FALSE,
        color = "white",
        size = dot_size
      )
  }
  
  p +
    scale_x_discrete(labels = tick_labels) +
    scale_fill_manual(values = meiosis_cols, drop = FALSE) +
    scale_color_manual(values = meiosis_cols, drop = FALSE) +
    labs(title = title_text, x = xlab, y = "Haploid chromosome number",
         fill = "Segregation type", color = "Segregation type") +
    coord_cartesian(ylim = c(0, NA)) +
    theme_classic() +
    theme(
      axis.text.x  = element_text(angle = x_angle, hjust = 1),
      axis.title.x = element_text(margin = margin(t = 10)),
      plot.margin  = margin(t = 6, r = 6, b = 18, l = 6),
      legend.position = if (show_legend) "right" else "none"
    )
}

# =========================
# DATA PREP
# =========================

# Mammalia
mamm_raw <- read_csv("Mammalia achiasmy.csv", show_col_types = FALSE) %>%
  mutate(across(any_of(c("Class","Subclass","Order","Family","Genus","Species")),
                ~ str_squish(as.character(.x)))) %>%
  mutate(
    Subclass2 = case_when(Subclass == "Theria" ~ "Eutheria", TRUE ~ Subclass),
    Meiosis3         = recode_meiosis3(`Meiosis Type`),
    `Diploid Number` = suppressWarnings(as.numeric(`Diploid Number`)),
    Haploid          = `Diploid Number` / 2,
    organism_id      = paste(Class, Subclass2, Order, Family, Genus, Species, sep="|")
  ) %>%
  filter(is.finite(Haploid), Haploid > 0,
         !is.na(Subclass2), Subclass2 %in% c("Metatheria","Eutheria"),
         !is.na(Meiosis3))

mamm_org <- mamm_raw %>%
  group_by(organism_id) %>%
  summarise(
    Haploid      = sample(Haploid, size = 1),
    Subclass2    = first(Subclass2),
    Meiosis_plot = collapse_meiosis3(Meiosis3),
    .groups = "drop"
  ) %>%
  filter(!is.na(Meiosis_plot))

# PANEL 1: switch so Chiasmatic group is LEFTMOST:
# order Eutheria (chiasmatic-only) then Metatheria
mamm_panel1 <- mamm_org %>%
  filter(Subclass2 == "Metatheria" | (Subclass2 == "Eutheria" & Meiosis_plot == "Chiasmatic")) %>%
  mutate(x_group = factor(Subclass2, levels = c("Eutheria","Metatheria")))
mamm1_n <- mamm_panel1 %>% count(x_group, name="n")

# PANEL 2 (already has Chiasmatic first)
euth_panel2 <- mamm_org %>%
  filter(Subclass2 == "Eutheria") %>%
  mutate(x_group = factor(as.character(Meiosis_plot),
                          levels = c("Chiasmatic","Achiasmatic","Asynaptic")))
euth2_n <- euth_panel2 %>% count(x_group, name="n")

# Insects
prep_insect <- function(file) {
  dat <- read_csv(file, show_col_types = FALSE) %>%
    mutate(across(any_of(c("Class","Order","Suborder","Family","Genus","Species")),
                  ~ str_squish(as.character(.x)))) %>%
    mutate(
      Meiosis3         = recode_meiosis3(`Meiosis Type`),
      `Diploid Number` = suppressWarnings(as.numeric(`Diploid Number`)),
      Haploid          = `Diploid Number` / 2,
      Suborder_clean   = str_to_lower(str_squish(as.character(Suborder)))
    ) %>%
    filter(is.finite(Haploid), Haploid > 0,
           !is.na(Suborder_clean), Suborder_clean != "",
           !is.na(Meiosis3)) %>%
    mutate(Suborder = case_when(
      Suborder_clean == "adephaga"  ~ "Adephaga",
      Suborder_clean == "polyphaga" ~ "Polyphaga",
      TRUE                          ~ Suborder
    ))
  
  id_cols <- intersect(c("Class","Order","Suborder","Family","Genus","Species"), names(dat))
  
  dat %>%
    mutate(organism_id = do.call(paste, c(across(all_of(id_cols)), sep="|"))) %>%
    group_by(organism_id) %>%
    summarise(
      Haploid      = sample(Haploid, size = 1),
      Suborder     = first(Suborder),
      Meiosis_plot = collapse_meiosis3(Meiosis3),
      .groups = "drop"
    ) %>%
    filter(!is.na(Meiosis_plot))
}

# Coleoptera
cole_org <- prep_insect("Coleoptera achiasmy.csv") %>%
  filter(Suborder %in% c("Polyphaga","Adephaga"))

# PANEL 3: switch so Chiasmatic (Adephaga) is LEFTMOST:
cole_panel3 <- cole_org %>%
  filter(
    (Suborder == "Polyphaga" & Meiosis_plot == "Asynaptic") |
      (Suborder == "Adephaga" & Meiosis_plot == "Chiasmatic")
  ) %>%
  mutate(x_group = factor(Suborder, levels = c("Adephaga","Polyphaga")))
cole3_n <- cole_panel3 %>% count(x_group, name="n")

# PANEL 4 (already has Chiasmatic first)
ade_panel4 <- cole_org %>%
  filter(Suborder == "Adephaga") %>%
  filter(Meiosis_plot != "Achiasmatic") %>%
  mutate(x_group = factor(as.character(Meiosis_plot),
                          levels = c("Chiasmatic","Asynaptic")))
ade_present <- ade_panel4 %>% distinct(x_group) %>% pull(x_group) %>% as.character()
ade_panel4  <- ade_panel4 %>%
  mutate(x_group = factor(as.character(x_group),
                          levels = c("Chiasmatic","Asynaptic")[c("Chiasmatic","Asynaptic") %in% ade_present]))
ade4_n <- ade_panel4 %>% count(x_group, name="n")

# Diptera
dipt_org <- prep_insect("Diptera achiasmy.csv")

# PANEL 5: ensure Chiasmatic Nematocera is left of Achiasmatic-only Brachycera
# (we force order: Nematocera, then Brachycera, then everything else)
dipt_panel5 <- dipt_org %>%
  mutate(Sub_l = str_to_lower(str_squish(as.character(Suborder)))) %>%
  filter(
    !(Sub_l == "nematocera" & Meiosis_plot == "Achiasmatic"),
    !(Sub_l == "brachycera" & Meiosis_plot != "Achiasmatic")
  )

dip_levels <- unique(as.character(dipt_panel5$Suborder))
front <- c("Nematocera","Brachycera")
dip_levels <- c(front[front %in% dip_levels], setdiff(dip_levels, front))

dipt_panel5 <- dipt_panel5 %>%
  mutate(x_group = factor(Suborder, levels = dip_levels))
dipt5_n <- dipt_panel5 %>% count(x_group, name="n")

# PANEL 6 (already has Chiasmatic first)
nem_panel6 <- dipt_org %>%
  filter(str_to_lower(str_squish(as.character(Suborder))) == "nematocera") %>%
  mutate(x_group = factor(as.character(Meiosis_plot),
                          levels = c("Chiasmatic","Achiasmatic","Asynaptic")))
nem_present <- nem_panel6 %>% distinct(x_group) %>% pull(x_group) %>% as.character()
nem_panel6  <- nem_panel6 %>%
  mutate(x_group = factor(as.character(x_group),
                          levels = c("Chiasmatic","Achiasmatic","Asynaptic")[c("Chiasmatic","Achiasmatic","Asynaptic") %in% nem_present]))
nem6_n <- nem_panel6 %>% count(x_group, name="n")

# =========================
# PER-PANEL BANDWIDTH (adjust) SETTINGS
# =========================
adj_p1 <- 0.8
adj_p2 <- 0.8
adj_p3 <- 1.0
adj_p4 <- 1.0
adj_p5 <- 1.8
adj_p6 <- 1.8

# =========================
# BUILD PLOTS
# =========================
p1 <- make_violin_custombox(mamm_panel1, mamm1_n, "Mammalia", "Subclass",
                            x_angle = 0, show_legend = FALSE,
                            overlay_dodge = TRUE, violin_adjust = adj_p1) +
  theme(axis.text.x = element_text(hjust = 0.5))

p2 <- make_violin_custombox(euth_panel2, euth2_n, "Eutheria", "Segregation type",
                            x_angle = 0, show_legend = TRUE,
                            overlay_dodge = FALSE, violin_adjust = adj_p2) +
  theme(axis.text.x = element_text(hjust = 0.5))

p3 <- make_violin_custombox(cole_panel3, cole3_n, "Coleoptera", "Suborder",
                            x_angle = 0, show_legend = FALSE,
                            overlay_dodge = TRUE, violin_adjust = adj_p3) +
  theme(axis.text.x = element_text(hjust = 0.5))

p4 <- make_violin_custombox(ade_panel4, ade4_n, "Adephaga", "Segregation type",
                            x_angle = 0, show_legend = FALSE,
                            overlay_dodge = FALSE, violin_adjust = adj_p4) +
  theme(axis.text.x = element_text(hjust = 0.5))

p5 <- make_violin_custombox(dipt_panel5, dipt5_n, "Diptera", "Suborder",
                            x_angle = 0, show_legend = FALSE,
                            overlay_dodge = TRUE, violin_adjust = adj_p5) +
  scale_y_continuous(breaks = scales::breaks_width(2)) +
  theme(axis.text.x = element_text(hjust = 0.5))

p6 <- make_violin_custombox(nem_panel6, nem6_n, "Nematocera", "Segregation type",
                            x_angle = 0, show_legend = FALSE,
                            overlay_dodge = FALSE, violin_adjust = adj_p6) +
  scale_y_continuous(breaks = scales::breaks_width(2)) +
  theme(axis.text.x = element_text(hjust = 0.5))

panel <- (p1 | p2) / (p3 | p4) / (p5 | p6) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right", plot.margin = margin(2, 2, 2, 2))

print(panel)

# =========================
# SAVE AS EMF
# =========================
devEMF::emf(
  file = "Violin_plot.emf",
  width = 9, height = 13,
  bg = "white",
  coordDPI = 300
)
print(panel)
dev.off()
