# ---- Packages ----
library(ape)
library(ggtree)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggtreeExtra)
library(ggnewscale)
library(scales)
library(ggplot2)

# =========================
# 1) Load + light cleaning
# =========================
setwd("~/GitHub/achiasmy-synthesis/data")
Taxa <- read_csv("Achiasmy_full_data.csv", show_col_types = FALSE)

clean <- function(x) x |>
  str_trim() |>
  str_replace_all("[\\(\\),;]", "_") |>
  str_squish()

Taxa <- Taxa %>%
  mutate(
    Kingdom   = clean(Kingdom),
    Class   = clean(Class),
    Order   = clean(Order),
    Family  = clean(Family),
    Genus   = clean(Genus),
    Species = clean(Species),
    Species_label = paste(Genus, Species, sep = "_"),
    Species_main  = sub(" ,*", "", Species),
    Species_key   = paste(Genus, Species_main, sep = "_"),
    path_key      = paste(Kingdom, Class, Order, Family, Genus, Species_main, sep = "|")
  ) %>%
  mutate(
    across(c(Kingdom, Class, Order, Family, Genus), as.factor),
    Species_key = factor(Species_key)
  )

# =========================
# 2) Choose ONE diploid per Species_key (before de-dup) + build path-level tree
# =========================
set.seed(1)  # reproducible tie-breaks

# ---- helpers (local) ----
.expand_range <- function(a, b) {
  a <- suppressWarnings(as.integer(a)); b <- suppressWarnings(as.integer(b))
  if (is.na(a) || is.na(b)) return(integer(0))
  if (b < a) { tmp <- a; a <- b; b <- tmp }
  seq.int(a, b)
}

parse_diploid_cell <- function(x) {
  if (is.na(x) || !nzchar(x)) return(integer(0))
  s <- x
  # normalize dashes; keep digits and separators
  s <- stringr::str_replace_all(s, "\u2013|\u2014|\u2212", "-")
  s <- stringr::str_replace_all(s, "[^0-9,;\\-]", " ")
  s <- stringr::str_squish(s)
  
  out <- integer(0)
  
  # expand ranges (e.g., 14-17 -> 14,15,16,17)
  rng <- stringr::str_match_all(s, "(?<!\\d)(\\d+)\\s*-\\s*(\\d+)(?!\\d)")[[1]]
  if (nrow(rng)) {
    rng_vals <- unlist(apply(rng[,2:3, drop=FALSE], 1, function(v) .expand_range(v[1], v[2])))
    out <- c(out, rng_vals)
    s <- stringr::str_replace_all(s, "(?<!\\d)\\d+\\s*-\\s*\\d+(?!\\d)", " ")
  }
  
  singles <- stringr::str_split(s, "[,;\\s]+", simplify = TRUE)
  singles <- singles[nzchar(singles)]
  if (length(singles)) {
    nums <- suppressWarnings(as.integer(singles))
    out  <- c(out, nums[is.finite(nums)])
  }
  unique(out)
}

choose_one_per_species <- function(df) {
  df %>%
    dplyr::select(Species_key, diploid_num_vals) %>%
    tidyr::unnest_longer(diploid_num_vals, values_to = "val", keep_empty = TRUE) %>%
    dplyr::filter(!is.na(val)) %>%
    dplyr::count(Species_key, val, name = "n") %>%
    dplyr::group_by(Species_key) %>%
    dplyr::mutate(n_max = max(n, na.rm = TRUE)) %>%
    dplyr::filter(n == n_max) %>%
    dplyr::slice_sample(n = 1) %>%   # random tie among modes
    dplyr::ungroup() %>%
    dplyr::transmute(Species_key, diploid_num_chosen = as.integer(val))
}

# ---- parse the *original* diploid column and pick one per species ----
diploid_col <- names(dplyr::select(
  Taxa, tidyselect::matches("^\\s*Diploid\\s+Number\\s*$")
))

diplo_candidates <- Taxa %>%
  dplyr::mutate(diploid_num_vals = lapply(.data[[diploid_col]], parse_diploid_cell)) %>%
  dplyr::select(Species_key, diploid_num_vals)

diplo_choice <- choose_one_per_species(diplo_candidates)  # Species_key, diploid_num_chosen

# ---- build the tree with *path_key* (keeps duplicates that differ up-tree) ----
Taxa_paths <- Taxa %>%
  dplyr::filter(!if_any(c(Kingdom, Class, Order, Family, Genus, Species_main), is.na)) %>%
  dplyr::arrange(Kingdom, Class, Order, Family, Genus, Species_key, path_key) %>%
  dplyr::distinct(path_key, .keep_all = TRUE) %>%
  dplyr::select(Kingdom, Class, Order, Family, Genus, Species_key, path_key) %>%
  dplyr::mutate(across(c(Kingdom, Class, Order, Family, Genus, path_key), as.factor))

phylo_tree <- as.phylo(~ Kingdom / Class / Order / Family / Genus / path_key, data = Taxa_paths)
phylo_tree$tip.label <- vapply(strsplit(phylo_tree$tip.label, "/"), function(z) z[length(z)], "")

phylo_tree <- as.phylo(~ Kingdom / Class / Order / Family / Genus / path_key, data = Taxa_paths)
phylo_tree$tip.label <- vapply(strsplit(phylo_tree$tip.label, "/"),
                               function(z) z[length(z)], "")
# NOTE: tip labels are now path_key

# =========================
# 3) Rename key columns (pattern-based; no helper fn)
# =========================
Taxa <- Taxa %>%
  rename(
    type_meiosis  = `Meiosis Type`,
    sex_no_recomb = `Achiasmatic Sex`,
    diploid_num   = `Diploid Number`,
    sex_chr_sys   = `Simplified SCS`
  )
# =========================
# 4) One row per *path* + join to tips (use chosen diploid per Species_key)
# =========================
meta_1path <- Taxa %>%
  dplyr::group_by(path_key) %>%
  dplyr::slice_head(n = 1) %>%                   # one metadata row per distinct path
  dplyr::ungroup() %>%
  dplyr::mutate(
    type_meiosis = dplyr::na_if(stringr::str_trim(type_meiosis), "")
  ) %>%
  dplyr::left_join(diplo_choice, by = "Species_key") %>%   # chosen per species
  dplyr::mutate(diploid_num = diploid_num_chosen) %>%
  dplyr::select(-diploid_num_chosen) %>%
  dplyr::filter(!is.na(type_meiosis))                      # drop NA if desired

keep_tips <- intersect(phylo_tree$tip.label, meta_1path$path_key)
drop_tips <- setdiff(phylo_tree$tip.label, keep_tips)
if (length(drop_tips)) phylo_tree <- drop.tip(phylo_tree, drop_tips)

tip_metadata <- tibble(
  label = factor(phylo_tree$tip.label, levels = phylo_tree$tip.label)  # label == path_key
) %>%
  dplyr::left_join(meta_1path, by = c("label" = "path_key")) %>%
  dplyr::mutate(diploid_num = tidyr::replace_na(as.integer(diploid_num), 0L))

# =========================
# 5) Recodes for rings
# =========================
tip_metadata <- tip_metadata %>%
  mutate(
    sex_no_recomb = case_when(
      str_detect(str_to_lower(sex_no_recomb), "^both(\\s+male\\s+and\\s+female)?\\??$") ~ "Both",
      TRUE ~ sex_no_recomb
    ),
    sex_no_recomb_ring = case_when(
      !is.na(sex_no_recomb) & str_to_lower(sex_no_recomb) %in% c("male","m","female","f") ~ "Heterogametic",
      !is.na(sex_no_recomb) & str_detect(str_to_lower(sex_no_recomb), "^both(\\s+male\\s+and\\s+female)?\\??$") ~ "Both",
      !is.na(sex_no_recomb) & str_detect(str_to_lower(sex_no_recomb), "hermaphrod") ~ "Hermaphroditic",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Both","Heterogametic","Hermaphroditic")),
    type_meiosis_ring = case_when(
      !is.na(type_meiosis) & str_detect(str_to_lower(type_meiosis), "^\\s*achiasmatic\\s*$") ~ "Achiasmatic",
      !is.na(type_meiosis) & str_detect(str_to_lower(type_meiosis),
                                        "achiasmatic\\s*\\(\\s*meiosis\\s*1\\s*\\).*distance\\s*pairing\\s*\\(\\s*meiosis\\s*2\\s*\\)") ~ "Achiasmatic",
      !is.na(type_meiosis) & str_detect(str_to_lower(type_meiosis), "distance\\s*pairing\\s*/\\s*achiasmatic") ~ "Achiasmatic",
      !is.na(type_meiosis) & str_detect(str_to_lower(type_meiosis), "achiasm") ~ "Achiasmatic",
      !is.na(type_meiosis) & str_detect(str_to_lower(type_meiosis), "distance\\s*pairing") ~ "Asynaptic",  # <- MERGE
      TRUE ~ type_meiosis
    ) %>% factor(levels = c("Achiasmatic", "Asynaptic"))
    )

# =========================
# 6) Tree grouping by big classes (+ colored branches)
#     Merge Copepoda into Maxillopoda for grouping/coloring
# =========================
grouped_tree <- ladderize(phylo_tree, right = TRUE)

# 6a) Create a merged class for grouping: Copepoda -> Maxillopoda
tip_metadata <- tip_metadata %>%
  mutate(
    Class_merge = ifelse(as.character(Class) == "Copepoda", "Maxillopoda", as.character(Class)),
    Class_merge = factor(Class_merge)
  )

# palettes computed on merged classes
classes_all <- sort(unique(as.character(tip_metadata$Class_merge)))
class_cols  <- setNames(scales::hue_pal(l = 70, c = 60)(length(classes_all)), classes_all)

# big groups (>=10 tips) using merged classes
classes_big     <- tip_metadata %>% count(Class_merge, name = "n") %>% filter(n >= 10) %>% pull(Class_merge) %>% as.character()
big_non_insecta <- setdiff(classes_big, "Insecta")

# build non-insecta groups by merged class
non_insecta_list <- tip_metadata %>%
  filter(Class_merge %in% big_non_insecta) %>%
  select(Class_merge, label) %>%
  mutate(across(everything(), as.character)) %>%
  split(.$Class_merge) %>% lapply(`[[`, "label")

# insecta subgroups unchanged
is_insecta <- str_to_lower(tip_metadata$Class_merge) == "insecta"
ord_lower  <- str_to_lower(tip_metadata$Order)
tips_coleop       <- as.character(tip_metadata$label[ is_insecta & ord_lower %in% c("coleoptera","cleoptera") ])
tips_diptera      <- as.character(tip_metadata$label[ is_insecta & ord_lower == "diptera" ])
tips_insecta_other <- setdiff(as.character(tip_metadata$label[is_insecta]), c(tips_coleop, tips_diptera))

insecta_list <- list(
  `Insecta—Coleoptera` = tips_coleop,
  `Insecta—Diptera`    = tips_diptera,
  `Insecta—other`      = tips_insecta_other
)

# explicit colors (note: Copepoda removed; Maxillopoda added)
explicit_groups <- c("Insecta—Coleoptera","Insecta—Diptera","Insecta—other",
                     "Mammalia","Arachnida","Maxillopoda")
dyn_non_insecta <- setdiff(big_non_insecta, explicit_groups)

# pick a hex for Maxillopoda; reuse your old Copepoda purple to keep consistency
branch_cols <- c(
  setNames(class_cols[dyn_non_insecta], dyn_non_insecta),  # dynamic others
  "Insecta—Coleoptera" = "#37C789",
  "Insecta—Diptera"    = "#27C8F5",
  "Insecta—other"      = "#8EF527",
  "Mammalia"           = "#F54927",
  "Arachnida"          = "#F5B040",
  "Maxillopoda"        = "#D658D4"
)

# rebuild branch_groups with the merged non-insecta lists
branch_groups <- c(non_insecta_list, insecta_list)

# =========================
# 7) Plot
# =========================
radius <- 3.5
tight_tree_theme <- ggtree::theme_tree() +
  theme(
    plot.margin  = margin(2, 2, 2, 2, "pt"),
    panel.grid   = element_blank(),
    panel.border = element_blank(),
    axis.text    = element_blank(),
    axis.title   = element_blank(),
    axis.ticks   = element_blank()
  )

branch_groups     <- c(non_insecta_list, insecta_list)
grouped_tree_col  <- ggtree::groupOTU(grouped_tree, branch_groups, group_name = "branch_group")

p0 <- (ggtree(grouped_tree_col, layout = "circular",
              branch.length = "none", open.angle = 0,
              size = 0, color = NA) %<+% tip_metadata) +
  xlim_tree(radius) + tight_tree_theme +
  scale_x_continuous(expand = expansion(mult = 0, add = 0))

p0$layers <- Filter(function(ly) !inherits(ly$geom, "GeomTree"), p0$layers)

S <- 4
p <- p0 +
  geom_tree(aes(color = branch_group),
            linewidth = 0.01 / S, lineend = "butt", linejoin = "round") +
  scale_color_manual(values = branch_cols, na.value = "lightgrey", name = "Branch group")

p <- open_tree(p, 20); p <- rotate_tree(p, 285)

# --- Ring sizing (thicker) ---
ring_w <- 0.5; gutter <- 0.012
off1 <- 0.10; off2 <- off1 + gutter; off3 <- off2 + gutter

# --- Ring 1 colors from BrBG: medium + dark yellow/tan ---
orangepal <- RColorBrewer::brewer.pal(9, "YlOrBr")
ring1_cols <- c(
  "Achiasmatic" = orangepal[4],  
  "Asynaptic"   = orangepal[6]  
)

# Ring 1 – Segregation Type 
p1 <- p + ggnewscale::new_scale_fill() +
  geom_fruit(
    geom  = geom_tile,
    aes(y = label, fill = type_meiosis_ring),
    pwidth    = ring_w,
    offset    = off1,
    color     = NA,
    linewidth = 0
  ) +
  scale_fill_manual(
    values       = ring1_cols,
    na.value     = "white",
    na.translate = FALSE,
    drop         = FALSE,
    name   = "Segregation Type",
    breaks = c("Achiasmatic", "Asynaptic"),
    labels = c(
      "Achiasmatic Meiosis/Sex Chromosomes",
      "Asynaptic Sex Chromosomes"
    ),
    guide  = guide_legend(order = 1)
  )


# Ring 2 – Sex without Recombination (SECOND)
p2 <- p1 + ggnewscale::new_scale_fill() +
  geom_fruit(
    geom = geom_tile, aes(y = label, fill = sex_no_recomb_ring),
    pwidth = ring_w, offset = off2, color = NA, linewidth = 0
  ) +
  scale_fill_brewer(
    palette      = "YlGn",
    direction    = -1,
    na.value     = "white",
    na.translate = FALSE,
    drop         = FALSE,
    name   = "Sex without Recombination",
    limits = c("Both","Heterogametic","Hermaphroditic"),
    guide  = guide_legend(order = 2)
  )

# =========================
# Ring 3 – Sex Chromosome System (THIRD)
# =========================
tip_metadata <- tip_metadata %>%
  mutate(
    sex_chr_sys = na_if(str_trim(as.character(sex_chr_sys)), ""),
    sex_chr_sys = ifelse(sex_chr_sys == "Complex XY systems", "XY", sex_chr_sys),
    sex_chr_sys = factor(sex_chr_sys, levels = c("XO","XY","ZO","ZW"))
  )

# --- use 4 colors from RdBu, but mapped to match your screenshot ---
sx_cols <- c(
  "XO" = "#2171B5",  
  "XY" = "#6BAED6",  
  "ZO" = "#c55215",  
  "ZW" = "#e28743"   
)
p3 <- p2 + ggnewscale::new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    aes(y = label, fill = sex_chr_sys),
    pwidth = ring_w, offset = off3, color = NA, linewidth = 0
  ) +
  scale_fill_manual(
    values       = sx_cols,
    na.value     = "white",
    na.translate = FALSE,
    drop         = FALSE,
    name         = "Sex Chromosome System",
    guide        = guide_legend(order = 3)
  )

# Outer bar ring
p_final <- p3 +
  geom_fruit(geom = geom_col, aes(y = label, x = diploid_num),
             orientation = "y", pwidth = 0.45, offset = 0.12, linewidth = 0)

# Legends
legend_theme_scaled <- theme(
  legend.position      = "right",
  legend.box           = "vertical",
  legend.justification = "top",
  legend.background    = element_blank(),
  legend.box.margin    = margin(4*S,4*S,4*S,4*S,"pt"),
  legend.margin        = margin(2*S,2*S,2*S,2*S,"pt"),
  legend.key.height    = unit(7*S,  "pt"),
  legend.key.width     = unit(8*S,  "pt"),
  legend.spacing.y     = unit(2*S,  "pt"),
  legend.text          = element_text(size = 5*S),
  legend.title         = element_text(size = 5*S, face = "bold")
)

p_out <- p_final + legend_theme_scaled + guides(color = "none")
p_out

# Save EMF (oversized; shrink in PPT to ~25% since S=4)
ggsave("~/GitHub/achiasmy-synthesis/figures/fulltree.emf",
       p_out,
       device = function(filename, ...) devEMF::emf(file = filename, emfPlus = TRUE, ...),
       width = dev.size("in")[1] * S,
       height = dev.size("in")[2] * S,
       bg = "white")
