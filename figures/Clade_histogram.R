# =========================
# 6-panel barplots with:
# - LESS whitespace between panels (esp. 1–4)
# - EACH panel has its OWN y-scale (so y-axis goes just above its max)
# - Layout: (p1|p2)/(p3|p4)/(p5|p6)
# =========================

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)

setwd("~/GitHub/achiasmy-synthesis/data")
Taxa <- read_csv("Achiasmy_full_data.csv", show_col_types = FALSE)

clean <- function(x) x |> str_trim() |> str_replace_all("[\\(\\),;]", "_") |> str_squish()

set.seed(1)

.expand_range <- function(a, b) {
  a <- suppressWarnings(as.integer(a)); b <- suppressWarnings(as.integer(b))
  if (is.na(a) || is.na(b)) return(integer(0))
  if (b < a) { tmp <- a; a <- b; b <- tmp }
  seq.int(a, b)
}

parse_diploid_cell <- function(x) {
  if (is.na(x) || !nzchar(x)) return(integer(0))
  s <- as.character(x)
  s <- str_replace_all(s, "\u2013|\u2014|\u2212", "-")
  s <- str_replace_all(s, "[^0-9,;\\-]", " ")
  s <- str_squish(s)
  
  out <- integer(0)
  
  rng <- str_match_all(s, "(?<!\\d)(\\d+)\\s*-\\s*(\\d+)(?!\\d)")[[1]]
  if (nrow(rng)) {
    out <- c(out, unlist(apply(rng[,2:3, drop = FALSE], 1,
                               function(v) .expand_range(v[1], v[2]))))
    s <- str_replace_all(s, "(?<!\\d)\\d+\\s*-\\s*\\d+(?!\\d)", " ")
  }
  
  singles <- str_split(s, "[,;\\s]+", simplify = TRUE)
  singles <- singles[nzchar(singles)]
  if (length(singles)) out <- c(out, suppressWarnings(as.integer(singles)))
  
  unique(out[is.finite(out)])
}

choose_one_per_species <- function(df) {
  df %>%
    select(Species_key, diploid_num_vals) %>%
    unnest_longer(diploid_num_vals, values_to = "val", keep_empty = TRUE) %>%
    filter(!is.na(val)) %>%
    count(Species_key, val, name = "n") %>%
    group_by(Species_key) %>%
    filter(n == max(n, na.rm = TRUE)) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    transmute(Species_key, diploid_num_chosen = as.integer(val))
}

if (!"Diploid Number" %in% names(Taxa)) stop("Column 'Diploid Number' not found.")

Taxa <- Taxa %>%
  mutate(
    Kingdom = clean(Kingdom),
    Class   = clean(Class),
    Order   = clean(Order),
    Family  = clean(Family),
    Genus   = clean(Genus),
    Species = clean(Species),
    Species_key = paste(Genus, sub(" ,*", "", Species), sep = "_"),
    path_key    = paste(Kingdom, Class, Order, Family, Genus, sub(" ,*", "", Species), sep = "|")
  )

diplo_choice <- Taxa %>%
  mutate(diploid_num_vals = lapply(.data[["Diploid Number"]], parse_diploid_cell)) %>%
  select(Species_key, diploid_num_vals) %>%
  choose_one_per_species()

branch_cols <- c(
  "Coleoptera" = "#37C789",
  "Diptera"    = "#27C8F5",
  "Misc. Insect Orders"      = "#8EF527",
  "Mammalia"           = "#F54927",
  "Arachnida"          = "#F5B040",
  "Maxillopoda"        = "#D658D4"
)

border_cols <- vapply(branch_cols, darken_hex, character(1))

dat <- Taxa %>%
  group_by(path_key) %>% slice_head(n = 1) %>% ungroup() %>%
  left_join(diplo_choice, by = "Species_key") %>%
  mutate(
    diploid_num = as.integer(diploid_num_chosen),
    haploid_num = diploid_num / 2
  ) %>%
  filter(is.finite(haploid_num), haploid_num > 0) %>%
  mutate(
    Class_merge = ifelse(as.character(Class) == "Copepoda", "Maxillopoda", as.character(Class)),
    is_insecta  = str_to_lower(Class_merge) == "insecta",
    clade = case_when(
      is_insecta & as.character(Order) == "Coleoptera" ~ "Coleoptera",
      is_insecta & as.character(Order) == "Diptera"    ~ "Diptera",
      is_insecta                                      ~ "Misc. Insect Orders",
      TRUE                                            ~ as.character(Class_merge)
    )
  ) %>%
  distinct(clade, Species_key, haploid_num) %>%
  filter(clade %in% names(branch_cols)) %>%
  group_by(clade) %>% filter(n_distinct(Species_key) >= 10) %>% ungroup()

bar_dat <- dat %>%
  mutate(h = as.integer(round(haploid_num))) %>%
  group_by(clade, h) %>%
  summarize(species_n = n_distinct(Species_key), .groups = "drop")

# n for titles
n_lookup <- dat %>% group_by(clade) %>% summarize(n = n_distinct(Species_key), .groups = "drop")
n_lookup <- setNames(n_lookup$n, n_lookup$clade)

# shared x range
xlims <- range(bar_dat$h, na.rm = TRUE)

# per-panel ymax so each plot's y-axis fits tightly
ymax_lookup <- bar_dat %>%
  group_by(clade) %>%
  summarize(ymax = max(species_n, na.rm = TRUE), .groups = "drop")
ymax_lookup <- setNames(ymax_lookup$ymax, ymax_lookup$clade)

theme_tight <- theme_classic(base_size = 12) +
  theme(
    plot.title  = element_text(face = "bold", size = 11, hjust = 0),
    axis.text   = element_text(size = 9),
    axis.title  = element_text(size = 11),
    # MUCH tighter margins + less panel whitespace
    plot.margin = margin(1, 1, 1, 1)
  )

make_bar <- function(cl, show_y = TRUE) {
  ymax <- ymax_lookup[[cl]]
  ggplot(filter(bar_dat, clade == cl), aes(h, species_n)) +
    geom_col(
      fill = branch_cols[[cl]],
      color = border_cols[[cl]],
      linewidth = 0.35,
      alpha = 0.75,
      width = 0.95
    ) +
    scale_x_continuous(limits = xlims, breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(
      limits = c(0, ymax * 1.05),              # <- tight per panel
      breaks = scales::pretty_breaks(n = 4),
      expand = c(0, 0)
    ) +
    labs(
      title = paste0(cl, " (n = ", n_lookup[[cl]], ")"),
      x = "Haploid chromosome number",
      y = if (show_y) "Species count" else NULL
    ) +
    theme_tight +
    theme(axis.title.y = if (show_y) element_text() else element_blank())
}

# fixed order (edit if you want)
clades6 <- c("Maxillopoda","Mammalia","Diptera","Coleoptera","Misc. Insect Orders","Arachnida")
clades6 <- clades6[clades6 %in% unique(bar_dat$clade)]

p1 <- make_bar(clades6[1], TRUE)
p2 <- make_bar(clades6[2], FALSE)
p3 <- make_bar(clades6[3], TRUE)
p4 <- make_bar(clades6[4], FALSE)
p5 <- make_bar(clades6[5], TRUE)
p6 <- make_bar(clades6[6], FALSE)

# tighter spacing between rows/cols
(p1 | p2) / (p3 | p4) / (p5 | p6) 
