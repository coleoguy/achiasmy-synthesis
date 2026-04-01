# ==============================================================================
# ACHIASMY SYNTHESIS: FULL MULTI-RING PHYLOGENETIC PIPELINE
# ==============================================================================

# ---- 0) Global Setup & Libraries ----
options(timeout = 600) 
library(ape); library(ggtree); library(dplyr); library(tidyr)
library(readr); library(stringr); library(ggtreeExtra); library(ggnewscale)
library(scales); library(ggplot2); library(taxize); library(RColorBrewer)
library(tidytree); library(igraph)

# ---- 1) Load Data & Prep Query ----
setwd("~/GitHub/achiasmy-synthesis/data")
Taxa <- read_csv("Achiasmy_full_meiosis_data.csv", show_col_types = FALSE)

# Standardize search string based on CSV columns
Taxa <- Taxa %>%
  mutate(
    query_name = case_when(
      !is.na(Species) & !is.na(Genus) ~ paste(Genus, Species),
      !is.na(Species) ~ Species,
      !is.na(Genus) ~ Genus,
      TRUE ~ Family
    )
  )

unique_names <- unique(na.omit(Taxa$query_name))
message(paste("Querying NCBI for", length(unique_names), "unique taxa..."))

# ---- 2) NCBI Backbone Lookup (Batch-Safe) ----
ids_ncbi <- get_ids(unique_names, db = "ncbi", ask = FALSE)
found_ids <- ids_ncbi$ncbi[!is.na(ids_ncbi$ncbi)]

# Query in chunks of 50 to avoid connection resets/throttling
hier_list <- list()
chunks <- split(found_ids, ceiling(seq_along(found_ids)/50))

message("Querying NCBI in batches...")
for(i in seq_along(chunks)){
  message(paste("Processing batch", i, "of", length(chunks)))
  hier_list[[i]] <- classification(chunks[[i]], db = "ncbi")
  Sys.sleep(1) # Safety pause to respect NCBI rate limits
}

hier_ncbi <- do.call(rbind, hier_list)

# ---- 2.5) Pivot NCBI Hierarchy into Backbone ----
bb_ncbi <- hier_ncbi %>%
  filter(rank %in% c("kingdom", "class", "order", "family", "genus", "species")) %>%
  pivot_wider(names_from = rank, values_from = name, values_fn = first)

# ---- 3) Standardize & Build Tree (Hardened Hierarchy) ----

# STEP A: Resolve Taxonomy and Create path_key
Taxa_resolved <- Taxa %>%
  left_join(bb_ncbi, by = c("query_name" = "query")) %>%
  mutate(
    f_kingdom = as.character(coalesce(kingdom, Kingdom)),
    f_class   = as.character(coalesce(class,   Class)),
    f_order   = as.character(coalesce(order,   Order)),
    f_family  = as.character(coalesce(family,  Family)),
    f_genus   = as.character(coalesce(genus,   Genus)),
    f_species = as.character(coalesce(species, Species)) 
  ) %>%
  mutate(
    f_class = ifelse(f_class %in% c("Copepoda", "Hexanauplia"), "Maxillopoda", f_class),
    f_order = ifelse(f_order == "Harpaticoida", "Harpacticoida", f_order),
    path_key = paste(f_kingdom, f_class, f_order, f_family, f_genus, f_species, sep="|"),
    Species_key = factor(paste(f_genus, f_species, sep="_"))
  )

# STEP B: Define Tree Builder Function
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
  g <- igraph::graph_from_data_frame(edge_table, directed = TRUE)
  return(tidytree::as.phylo(g))
}

# STEP C: Build Hierarchy Table with Placeholders
hierarchy <- Taxa_resolved %>%
  mutate(across(c(f_kingdom, f_class, f_order, f_family, f_genus), 
                ~ifelse(is.na(.), paste0("Unknown_", cur_column()), .))) %>%
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

# STEP D: Execute Build and Ladderize
phylo_tree <- make_phylo_safely(hierarchy)
if (!is.null(phylo_tree)) {
  actual_data_tips <- intersect(phylo_tree$tip.label, Taxa_resolved$path_key)
  phylo_tree <- ape::keep.tip(phylo_tree, actual_data_tips)
  phylo_tree <- ape::collapse.singles(phylo_tree)
  phylo_tree <- ape::multi2di(phylo_tree)
}
grouped_tree <- ladderize(phylo_tree, right = TRUE)

# ---- 4) Calculate Diploid Numbers ----
set.seed(1)
parse_diploid_cell <- function(x) {
  if (is.na(x) || !nzchar(x)) return(integer(0))
  s <- str_replace_all(x, "\u2013|\u2014|\u2212", "-") %>% 
    str_replace_all("[^0-9,;\\-]", " ") %>% str_squish()
  nums <- suppressWarnings(as.integer(str_split(s, "[,;\\s]+", simplify = TRUE)))
  unique(nums[is.finite(nums)])
}

diplo_col <- names(Taxa_resolved)[grep("Diploid", names(Taxa_resolved), ignore.case = TRUE)][1]

diplo_choice <- Taxa_resolved %>%
  mutate(diploid_num_vals = lapply(.data[[diplo_col]], parse_diploid_cell)) %>%
  unnest_longer(diploid_num_vals, values_to = "val", keep_empty = TRUE) %>%
  filter(!is.na(val)) %>% 
  count(Species_key, val) %>%
  group_by(Species_key) %>% 
  filter(n == max(n)) %>% slice_sample(n = 1) %>% ungroup() %>% 
  transmute(Species_key, diploid_num_chosen = as.integer(val))

# ---- 5) Metadata Alignment & Factor Logic (Fixed for Dots) ----
tip_metadata <- tibble(label = as.character(phylo_tree$tip.label)) %>%
  left_join(Taxa_resolved %>% 
              rename(type_meiosis = Meiosis.Type, 
                     sex_no_recomb = Achiasmatic.Sex, 
                     scs_raw = Simplified.SCS) %>%
              select(path_key, type_meiosis, sex_no_recomb, scs_raw, f_class, f_order, Species_key) %>%
              distinct(path_key, .keep_all = TRUE), by = c("label" = "path_key")) %>%
  left_join(diplo_choice, by = "Species_key") %>%
  mutate(
    meiosis_raw_clean = str_to_lower(str_trim(type_meiosis)),
    
    type_meiosis_ring = case_when(
      str_detect(meiosis_raw_clean, "achiasmatic") & !str_detect(meiosis_raw_clean, "sex") ~ "Achiasmatic Meiosis",
      str_detect(meiosis_raw_clean, "achiasmatic") & str_detect(meiosis_raw_clean, "sex") ~ "Achiasmatic Sex Chromosomes",
      str_detect(meiosis_raw_clean, "asynap|distance") ~ "Asynaptic Sex Chromosomes",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Achiasmatic Meiosis", "Asynaptic Sex Chromosomes", "Achiasmatic Sex Chromosomes")),
    
    sex_no_recomb_ring = case_when(
      str_to_lower(sex_no_recomb) %in% c("male","m","female","f") ~ "Heterogametic",
      str_detect(str_to_lower(sex_no_recomb), "^both") ~ "Both",
      str_detect(str_to_lower(sex_no_recomb), "hermaphrod") ~ "Hermaphroditic",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Both","Heterogametic","Hermaphroditic")),
    
    sex_chr_sys = case_when(
      scs_raw == "XO" ~ "XO",
      scs_raw %in% c("XY", "Complex XY systems") ~ "XY",
      scs_raw == "ZO" ~ "ZO",
      scs_raw == "ZW" ~ "ZW",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("XO","XY","ZO","ZW"))
  )

# ---- 6) Define Colors & Exclusive Branch Logic ----
is_insecta <- str_to_lower(tip_metadata$f_class) == "insecta"
ord_lower  <- str_to_lower(tip_metadata$f_order)

tips_coleop  <- as.character(tip_metadata$label[is_insecta & ord_lower %in% c("coleoptera","cleoptera")])
tips_diptera <- as.character(tip_metadata$label[is_insecta & ord_lower == "diptera"])
tips_ins_oth <- setdiff(as.character(tip_metadata$label[is_insecta]), c(tips_coleop, tips_diptera))
tips_mammal  <- as.character(tip_metadata$label[tip_metadata$f_class == "Mammalia"])
tips_arach   <- as.character(tip_metadata$label[tip_metadata$f_class == "Arachnida"])
tips_maxillo <- as.character(tip_metadata$label[tip_metadata$f_class == "Maxillopoda"])

node_data <- data.frame(node = 1:(length(grouped_tree$tip.label) + grouped_tree$Nnode)) %>%
  mutate(branch_color = "gray90", branch_size = 0.05, assigned = FALSE) 

group_list <- list("Insecta-Coleoptera" = tips_coleop, "Insecta-Diptera" = tips_diptera,
                   "Mammalia" = tips_mammal, "Arachnida" = tips_arach, 
                   "Maxillopoda" = tips_maxillo, "Insecta-other" = tips_ins_oth)

group_colors <- c("Insecta-Coleoptera" = "#37C789", "Insecta-Diptera" = "#27C8F5",
                  "Mammalia" = "#F54927", "Arachnida" = "#F5B040", 
                  "Maxillopoda" = "#D658D4", "Insecta-other" = "#8EF527")

root_node <- length(grouped_tree$tip.label) + 1
for (group_name in names(group_list)) {
  tips_in_tree <- intersect(group_list[[group_name]], grouped_tree$tip.label)
  if (length(tips_in_tree) == 0) next
  m_node <- if(length(tips_in_tree) > 1) ape::getMRCA(grouped_tree, tips_in_tree) else which(grouped_tree$tip.label == tips_in_tree)
  if (!is.na(m_node)) {
    path <- tryCatch({ nodepath(grouped_tree, from = root_node, to = m_node) }, error = function(e) NULL)
    clade <- c(m_node, tidytree::offspring(grouped_tree, m_node))
    full_path <- unique(c(path, clade))
    target_nodes <- if(group_name == "Insecta-other") full_path[!node_data$assigned[node_data$node %in% full_path]] else full_path
    node_data$branch_color[node_data$node %in% target_nodes] <- group_colors[group_name]
    node_data$branch_size[node_data$node %in% target_nodes]  <- 0.2 
    node_data$assigned[node_data$node %in% target_nodes]     <- TRUE
  }
}

# ---- 7) Final Visualization (Outline-Free & Vector-Ready) ----
radius <- 20; ring_w_thick <- 1.2; offset_val <- 0.1
orangepal <- brewer.pal(9, "Oranges")
sx_cols <- c("XO" = "#E41A1C", "XY" = "#377EB8", "ZO" = "#4DAF4A", "ZW" = "#984EA3")

# Set color to transparent to kill the black background layer
p <- ggtree(grouped_tree, layout = "circular", branch.length = "none", color = "transparent") %<+% 
  node_data %<+% tip_metadata +
  geom_tree(aes(color = I(branch_color), linewidth = I(branch_size))) +
  xlim_tree(radius + 150) + theme_tree()

p <- open_tree(p, 20) %>% rotate_tree(285)

# Ring 1: Segregation Type
p <- p + new_scale_fill() +
  geom_fruit(geom = geom_tile, aes(y = label, fill = type_meiosis_ring), 
             pwidth = ring_w_thick, offset = offset_val, color = NA) +
  scale_fill_manual(values = c("Achiasmatic Meiosis"=orangepal[3], 
                               "Asynaptic Sex Chromosomes"=orangepal[8], 
                               "Achiasmatic Sex Chromosomes"=orangepal[5]), name = "Segregation Type") +
  guides(fill = guide_legend(order = 1))

# Ring 2: Sex w/o Recomb
p <- p + new_scale_fill() +
  geom_fruit(geom = geom_tile, aes(y = label, fill = sex_no_recomb_ring), 
             pwidth = ring_w_thick, offset = offset_val, color = NA) +
  scale_fill_brewer(palette = "YlGn", direction = -1, name = "Sex w/o Recomb.") +
  guides(fill = guide_legend(order = 2))

# Ring 3: SCS
p <- p + new_scale_fill() +
  geom_fruit(geom = geom_tile, aes(y = label, fill = sex_chr_sys), 
             pwidth = ring_w_thick, offset = offset_val, color = NA) +
  scale_fill_manual(values = sx_cols, name = "Sex Chromosome System") +
  guides(fill = guide_legend(order = 3))

p_final <- p + 
  geom_fruit(geom = geom_col, aes(y = label, x = as.numeric(diploid_num_chosen)), 
             orientation = "y", pwidth = 0.5, offset = 0.2, fill = "grey30") +
  theme(legend.position = "right", legend.box = "vertical", legend.spacing.y = unit(0.2, "cm"))

p_final

# ---- 8) Fail-Safe Export to CSV ----
final_data_export <- Taxa_resolved %>%
  filter(path_key %in% grouped_tree$tip.label) %>%
  left_join(tip_metadata %>% 
              select(label, type_meiosis_ring, sex_no_recomb_ring, sex_chr_sys, diploid_num_chosen), 
            by = c("path_key" = "label")) %>%
  transmute(Kingdom = f_kingdom, Class = f_class, Order = f_order, Family = f_family, Genus = f_genus, Species = f_species,
            Full_Path_Key = path_key, Segregation_Type = type_meiosis_ring, Sex_Without_Recomb = sex_no_recomb_ring,
            Sex_Chromosome_Sys = sex_chr_sys, Diploid_Number = diploid_num_chosen) %>%
  arrange(Kingdom, Class, Order, Family, Genus, Species) %>% 
  distinct()

write_csv(final_data_export, "Standardized_Tree_Data_Full.csv")
message("Pipeline Complete. Figure generated and CSV exported.")