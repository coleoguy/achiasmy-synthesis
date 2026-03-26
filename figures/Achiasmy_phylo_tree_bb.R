# ==============================================================================
# ACHIASMY SYNTHESIS: MULTI-DATABASE BACKBONE & NUANCE-SAFE TREE PIPELINE
# ==============================================================================

# ---- 0) Setup ----
options(timeout = 600) 
library(ape); library(ggtree); library(dplyr); library(tidyr)
library(readr); library(stringr); library(ggtreeExtra); library(ggnewscale)
library(scales); library(ggplot2); library(taxize)

# =========================
# 1) Load & Rescuing Names
# =========================
# Set your working directory and load data
setwd("~/GitHub/achiasmy-synthesis/data")
Taxa <- read_csv("Achiasmy_full_meiosis_data.csv", show_col_types = FALSE)

clean_val <- function(x) x |> str_trim() |> str_replace_all("[\\(\\),;]", "_") |> str_squish()

Taxa <- Taxa %>%
  mutate(
    # display_name preserves original nuances like "sp. 1" or "cytotype I"
    display_name = clean_val(Species),
    
    # query_species extracts the canonical name (e.g., "tincta" from "tincta?")
    query_species = display_name %>%
      str_replace_all("(?i)sp\\.|ssp\\.|nr\\.|[0-9]+", "") %>%
      str_replace_all("[\\?,\\(\\)]", "") %>%
      str_trim() %>%
      str_extract("^[A-Za-z]+"),
    
    # Search string for NCBI/GBIF
    query_name = paste(Genus, query_species, sep = " ")
  )

# ========================================================
# 2) Multi-Backbone Lookup (NCBI + GBIF)
# ========================================================
unique_names <- unique(Taxa$query_name)
message(paste("Querying databases for", length(unique_names), "taxa..."))

# NCBI Fetch
ids_ncbi <- get_ids(unique_names, db = "ncbi", ask = FALSE)
hier_ncbi <- classification(ids_ncbi$ncbi, db = "ncbi", pause = 0.5)
bb_ncbi <- rbind(hier_ncbi) %>%
  filter(rank %in% c("kingdom", "class", "order", "family", "genus", "species")) %>%
  pivot_wider(names_from = rank, values_from = name, values_fn = first)

# ========================================================
# 3) Merge, Standardize, and Hierarchy Prep
# ========================================================

Taxa_resolved <- Taxa %>%
  left_join(bb_ncbi, by = c("query_name" = "query")) %>%
  mutate(
    f_kingdom = as.character(coalesce(kingdom, Kingdom)),
    f_class   = as.character(coalesce(class,   Class)),
    f_order   = as.character(coalesce(order,   Order)),
    f_family  = as.character(coalesce(family,  Family)),
    f_genus   = as.character(coalesce(genus,   Genus)),
    f_species = as.character(coalesce(species, display_name))
  ) %>%
  mutate(
    f_class = ifelse(f_class %in% c("Copepoda", "Hexanauplia"), "Maxillopoda", f_class),
    f_order = ifelse(f_order == "Harpaticoida", "Harpacticoida", f_order),
    path_key = paste(f_kingdom, f_class, f_order, f_family, f_genus, display_name, sep="|"),
    Species_key = factor(paste(f_genus, display_name, sep="_"))
  )

# ========================================================
# 4) Build Tree (Hardened Hierarchy)
# ========================================================

# STEP A: Define the builder function
make_phylo_safely <- function(df) {
  edges <- list()
  cols <- colnames(df)
  for (i in 1:(length(cols) - 1)) {
    tmp <- df[, c(cols[i], cols[i+1])]
    colnames(tmp) <- c("parent", "child")
    edges[[i]] <- distinct(tmp)
  }
  edge_table <- bind_rows(edges) %>%
    filter(!is.na(parent) & !is.na(child)) %>% 
    mutate(across(everything(), as.character))
  
  # Crucial: Ensure the graph is connected before as.phylo
  g <- igraph::graph_from_data_frame(edge_table, directed = TRUE)
  return(tidytree::as.phylo(g))
}

# STEP B: Build the HARDENED Hierarchy 
# (This ensures no NAs break the tree chain)
hierarchy <- Taxa_resolved %>%
  mutate(
    f_kingdom = ifelse(is.na(f_kingdom), "Unknown_Kingdom", f_kingdom),
    f_class   = ifelse(is.na(f_class),   "Unknown_Class",   f_class),
    f_order   = ifelse(is.na(f_order),   "Unknown_Order",   f_order),
    f_family  = ifelse(is.na(f_family),  "Unknown_Family",  f_family),
    f_genus   = ifelse(is.na(f_genus),   "Unknown_Genus",   f_genus)
  ) %>%
  mutate(
    root      = "Life",
    v_kingdom = f_kingdom,
    v_class   = paste(v_kingdom, f_class, sep="_"),
    v_order   = paste(v_class, f_order, sep="_"),
    v_family  = paste(v_order, f_family, sep="_"),
    v_genus   = paste(v_family, f_genus, sep="_"),
    v_tip     = as.character(path_key) 
  ) %>%
  select(root, v_kingdom, v_class, v_order, v_family, v_genus, v_tip) %>%
  distinct()

# STEP C: Execute the Build
phylo_tree <- make_phylo_safely(hierarchy)

# STEP D: Prune and Clean Geometry
if (!is.null(phylo_tree)) {
  actual_data_tips <- intersect(phylo_tree$tip.label, Taxa_resolved$path_key)
  if (length(actual_data_tips) > 0) {
    phylo_tree <- ape::keep.tip(phylo_tree, actual_data_tips)
    phylo_tree <- ape::collapse.singles(phylo_tree)
    phylo_tree <- ape::multi2di(phylo_tree)
    message(paste("Success! Tree built with", length(phylo_tree$tip.label), "tips."))
  }
}
grouped_tree <- ladderize(phylo_tree, right = TRUE)

# ========================================================
# 4.5) Calculate Diploid Numbers (MUST RUN BEFORE SEC 5)
# ========================================================
set.seed(1)
parse_diploid_cell <- function(x) {
  if (is.na(x) || !nzchar(x)) return(integer(0))
  s <- str_replace_all(x, "\u2013|\u2014|\u2212", "-") %>% 
    str_replace_all("[^0-9,;\\-]", " ") %>% str_squish()
  nums <- suppressWarnings(as.integer(str_split(s, "[,;\\s]+", simplify = TRUE)))
  unique(nums[is.finite(nums)])
}

diplo_col_name <- names(Taxa_resolved)[grep("Diploid", names(Taxa_resolved), ignore.case = TRUE)][1]

diplo_choice <- Taxa_resolved %>%
  mutate(diploid_num_vals = lapply(.data[[diplo_col_name]], parse_diploid_cell)) %>%
  unnest_longer(diploid_num_vals, values_to = "val", keep_empty = TRUE) %>%
  filter(!is.na(val)) %>% 
  count(Species_key, val) %>%
  group_by(Species_key) %>% 
  filter(n == max(n)) %>% 
  slice_sample(n = 1) %>% 
  ungroup() %>% 
  transmute(Species_key, diploid_num_chosen = as.integer(val))

# ========================================================
# 5) Metadata Alignment & Grouping (Locking the Order)
# ========================================================
# Assign physical sequence: Copepoda -> Mammalia -> Arachnida -> Insects
tip_metadata <- tibble(label = as.character(phylo_tree$tip.label)) %>%
  left_join(
    Taxa_resolved %>% 
      rename(type_meiosis = `Meiosis Type`, sex_no_recomb = `Achiasmatic Sex`) %>%
      select(path_key, type_meiosis, sex_no_recomb, `Simplified SCS`, f_class, f_order, Species_key) %>%
      distinct(path_key, .keep_all = TRUE), 
    by = c("label" = "path_key")
  ) %>%
  left_join(diplo_choice, by = "Species_key") %>%
  mutate(
    clade_group = case_when(
      str_to_lower(f_class) == "maxillopoda"                ~ "Copepoda",
      str_to_lower(f_class) == "mammalia"                   ~ "Mammalia",
      str_to_lower(f_class) == "arachnida"                  ~ "Arachnida",
      str_to_lower(f_order) %in% c("coleoptera", "cleoptera") ~ "Coleoptera",
      str_to_lower(f_order) == "diptera"                    ~ "Diptera",
      str_to_lower(f_class) == "insecta"                    ~ "Insecta-other",
      TRUE                                                  ~ "Other"
    ) %>% factor(levels = c("Copepoda", "Mammalia", "Arachnida", 
                            "Insecta-other", "Coleoptera", "Diptera")),
    
    diploid_num = tidyr::replace_na(as.integer(diploid_num_chosen), 0L),
    meiosis_raw_clean = str_to_lower(str_trim(type_meiosis)),
    type_meiosis_ring = case_when(
      str_detect(meiosis_raw_clean, "achiasmatic meiosis") ~ "Achiasmatic Meiosis",
      str_detect(meiosis_raw_clean, "achiasmatic sex") ~ "Achiasmatic Sex Chromosomes",
      str_detect(meiosis_raw_clean, "asynap|distance") ~ "Asynaptic Sex Chromosomes",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Asynaptic Sex Chromosomes", "Achiasmatic Meiosis", "Achiasmatic Sex Chromosomes")),
    sex_no_recomb_ring = case_when(
      str_to_lower(sex_no_recomb) %in% c("male","m","female","f") ~ "Heterogametic",
      str_detect(str_to_lower(sex_no_recomb), "^both") ~ "Both",
      str_detect(str_to_lower(sex_no_recomb), "hermaphrod") ~ "Hermaphroditic",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Both","Heterogametic","Hermaphroditic")),
    sex_chr_sys = case_when(
      `Simplified SCS` == "XO" ~ "XO",
      `Simplified SCS` %in% c("XY", "Complex XY systems") ~ "XY",
      `Simplified SCS` == "ZO" ~ "ZO",
      `Simplified SCS` == "ZW" ~ "ZW",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("XO","XY","ZO","ZW"))
  )

# ========================================================
# 6) Deterministic Physical Reordering & Branch Coloring
# ========================================================
# 1. Generate master target sequence
target_tip_sequence <- tip_metadata %>% arrange(clade_group) %>% pull(label)

# 2. PHYSICALLY REORDER THE TREE
grouped_tree <- ape::rotateConstr(phylo_tree, target_tip_sequence)

# 3. CRITICAL SYNC: Re-index Metadata to the physical tree order
tip_metadata <- tip_metadata[match(grouped_tree$tip.label, tip_metadata$label), ]

# 4. Initialize Node Data
all_nodes_count <- length(grouped_tree$tip.label) + grouped_tree$Nnode
root_node <- length(grouped_tree$tip.label) + 1
node_data <- data.frame(node = 1:all_nodes_count) %>%
  mutate(branch_color = "gray90", branch_size = 0.05)

# 5. Define group colors (Match the names in clade_group above)
group_colors <- c(
  "Copepoda"      = "#D658D4",
  "Mammalia"      = "#F54927",
  "Arachnida"     = "#F5B040",
  "Insecta-other" = "#8EF527", 
  "Coleoptera"    = "#37C789",
  "Diptera"       = "#27C8F5"
)

# 6. Branch Coloring Loop
final_groups <- split(tip_metadata$label, tip_metadata$clade_group)
for (group_name in names(final_groups)) {
  tips_in_tree <- intersect(final_groups[[group_name]], grouped_tree$tip.label)
  if (length(tips_in_tree) == 0) next
  
  m_node <- if(length(tips_in_tree) > 1) {
    ape::getMRCA(grouped_tree, tips_in_tree)
  } else {
    which(grouped_tree$tip.label == tips_in_tree)
  }
  
  if (!is.null(m_node) && !is.na(m_node)) {
    # Color the 'stem' back to root
    path_to_root <- tryCatch({ nodepath(grouped_tree, from = root_node, to = m_node) }, error = function(e) NULL)
    clade_nodes <- c(m_node, tidytree::offspring(grouped_tree, m_node))
    full_path <- unique(c(path_to_root, clade_nodes))
    node_data$branch_color[node_data$node %in% full_path] <- group_colors[group_name]
    node_data$branch_size[node_data$node %in% full_path]  <- 0.2 
  }
}

# ========================================================
# 7) Plotting 
# ========================================================
# 1. Plot Constants (Recalibrated for a radius of 20)
# 'pwidth' controls the thickness of the ring.
# 'offset' is the distance from the previous layer.
radius <- 20
ring_w_thick <- 8  # Much smaller for vector precision
bar_w_thick  <- 0.5   # Width of the outer diploid bar chart
library(RColorBrewer)

# Ring 1 Colors
meiosis_cols <- c(
  "Asynaptic Sex Chromosomes" = "#FEC44F",
  "Achiasmatic Meiosis"= "#EC7014",
  "Achiasmatic Sex Chromosomes" = "#A22A58"
)

# Ring 3 Colors (Sex Chromosome Systems)
sx_cols <- c(
  "XO" = "#2171B5", 
  "XY" = "#6BAED6", 
  "ZO" = "#C55215", 
  "ZW" = "#E28743"  
)

# 2. Base Tree Construction
p <- ggtree(grouped_tree, layout = "circular", branch.length = "none", color="transparent") %<+% 
  node_data %<+% 
  tip_metadata +
  geom_tree(aes(color = I(branch_color), linewidth = I(branch_size))) +
  xlim_tree(radius + 40) + # Increased xlim to make room for 3 rings + bars
  theme_tree()

# 3. Open and Rotate
p <- open_tree(p, 20) %>% rotate_tree(285)

# --- RING 1: Segregation Type (The Orangepal Ring) ---
p <- p + new_scale_fill() +
  geom_fruit(geom = geom_tile, aes(y = label, fill = type_meiosis_ring), 
             pwidth = ring_w_thick, offset = 0.1, color = NA) +
  scale_fill_manual(
    values = meiosis_cols,
    name = "Alternative Mieosis Type",
    na.translate = FALSE
  )

# --- RING 2: Sex w/o Recomb (RESTORED) ---
p <- p + new_scale_fill() +
  geom_fruit(geom = geom_tile, aes(y = label, fill = sex_no_recomb_ring), 
             pwidth = ring_w_thick, offset = 0.08, color = NA) +
  scale_fill_brewer(palette = "YlGn", direction = -1, 
                    name = "Sex without Recombination", 
                    na.translate = FALSE)

# --- RING 3: Sex Chromosome System (RESTORED) ---
p <- p + new_scale_fill() +
  geom_fruit(geom = geom_tile, aes(y = label, fill = sex_chr_sys), 
             pwidth = ring_w_thick, offset = 0.08, color = NA) +
  scale_fill_manual(values = sx_cols, 
                    name = "Sex Chromosome System", 
                    na.translate = FALSE)

# --- FINAL LAYER: Diploid Number Bars ---
p_final <- p + 
  geom_fruit(geom = geom_col, 
             aes(y = label, x = as.numeric(diploid_num_chosen)), 
             orientation = "y", 
             pwidth = bar_w_thick, 
             offset = 0.08, # Extra gap before the bars
             fill = "grey30")

# Preview
p_final
# ========================================================
# 8) Export - Use Standard EMF (Most stable)
# ========================================================
devEMF::emf(file = "~/GitHub/achiasmy-synthesis/figures/fulltree_meiosis.emf", 
            width = 12, height = 12, 
            emfPlus = TRUE) # FALSE prevents PPT from trying to 'help' with rendering
print(p_final)
dev.off()

# ========================================================
# 9) Export Final Standardized Dataset (Full Hierarchy)
# ========================================================

# 1. Start with the master resolved data
final_data_export <- Taxa_resolved %>%
  # 2. Only keep species that are actually tips in your final tree
  filter(path_key %in% grouped_tree$tip.label) %>%
  # 3. Join the Ring/Bar data from tip_metadata
  left_join(
    tip_metadata %>% 
      select(label, type_meiosis_ring, sex_no_recomb_ring, sex_chr_sys, diploid_num_chosen), 
    by = c("path_key" = "label")
  ) %>%
  # 4. Create the clean hierarchy structure
  transmute(
    Kingdom              = f_kingdom,
    Class                = f_class,
    Order                = f_order,
    Family               = f_family,
    Genus                = f_genus,
    Species              = f_species,
    Full_Path_Key        = path_key,
    
    # Ring Data (Biological Traits)
    Segregation_Type     = type_meiosis_ring,
    Sex_Without_Recomb   = sex_no_recomb_ring,
    Sex_Chromosome_Sys   = sex_chr_sys,
    
    # Bar Data
    Diploid_Number       = diploid_num_chosen
  ) %>%
  # 5. Sort for a professional spreadsheet look
  arrange(Kingdom, Class, Order, Family, Genus, Species) %>%
  distinct()

# 6. Final Write
write_csv(final_data_export, 
          file = "~/GitHub/achiasmy-synthesis/data/Standardized_Tree_Data_Full.csv")

message("Export Success! Check your data folder for Standardized_Tree_Data_Full.csv")