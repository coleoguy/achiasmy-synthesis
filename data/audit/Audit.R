# ==============================================================================
# ACHIASMY SYNTHESIS: THREE-WAY DATA CENSUS & ACCOUNTABILITY PIPELINE
# ==============================================================================
library(dplyr); library(readr); library(tidyr); library(stringr)

# 1. Load the three snapshots
raw_df      <- read_csv("~/GitHub/achiasmy-synthesis/data/Achiasmy_full_data.csv", show_col_types = FALSE)
backbone_df <- read_csv("~/GitHub/achiasmy-synthesis/data/backbone_tree_tip_metadata.csv", show_col_types = FALSE)
final_df    <- read_csv("~/GitHub/achiasmy-synthesis/data/final_tree_tip_labels.csv", show_col_types = FALSE)

# 2. Standardize Full Path Keys for Comparison
# Raw: Original manual hierarchy
raw_paths <- raw_df %>%
  mutate(full_path = paste(Kingdom, Class, Order, Family, Genus, Species, sep="|")) %>%
  pull(full_path) %>% unique()

# Method A (Backbone): Database-resolved hierarchy (f_ columns)
backbone_paths <- backbone_df %>%
  mutate(full_path = paste(f_kingdom, f_class, f_order, f_family, f_genus, display_name, sep="|")) %>%
  pull(full_path) %>% unique()

# Method B (Final): Current memory-safe tree hierarchy
final_paths <- final_df %>%
  mutate(full_path = paste(Kingdom, Class, Order, Family, Genus, Species, sep="|")) %>%
  pull(full_path) %>% unique()

# 3. intersection Logic: Finding what was lost where
shared_all        <- intersect(intersect(raw_paths, backbone_paths), final_paths)
lost_raw_to_back  <- setdiff(raw_paths, backbone_paths)  # Paths dropped/changed at Method A
lost_back_to_fin  <- setdiff(backbone_paths, final_paths) # Paths in A but not B (the 67)
lost_raw_to_fin   <- setdiff(raw_paths, final_paths)      # Paths in Raw but not B (the 85)

# 4. Generate Census Table
census_summary <- tibble(
  Metric = c(
    "1. Unique Paths in Raw CSV",
    "2. Unique Paths in Backbone (Method A)",
    "3. Unique Paths in Final (Method B)",
    "4. Paths Identical across all 3",
    "5. Paths in Backbone but Lost/Changed in Final",
    "6. Raw Paths Missing from Final Tree",
    "7. Raw Paths Missing from Backbone"
  ),
  Count = c(
    length(raw_paths),
    length(backbone_paths),
    length(final_paths),
    length(shared_all),
    length(lost_back_to_fin),
    length(lost_raw_to_fin),
    length(lost_raw_to_back)
  )
)

# 5. Extract the Specific 85 Missing Raw Paths
# These are paths in raw data that didn't appear in Method A OR Method B
the_missing_85 <- raw_df %>%
  mutate(full_path = paste(Kingdom, Class, Order, Family, Genus, Species, sep="|")) %>%
  filter(!full_path %in% final_paths) %>%
  filter(!full_path %in% backbone_paths) %>%
  select(Kingdom, Class, Order, Family, Genus, Species, full_path) %>%
  distinct()

# 6. Extract the 67 Paths Lost between Backbone and Final
the_lost_67 <- backbone_df %>%
  mutate(full_path = paste(f_kingdom, f_class, f_order, f_family, f_genus, display_name, sep="|")) %>%
  filter(!full_path %in% final_paths) %>%
  select(f_kingdom, f_class, f_order, f_family, f_genus, display_name, full_path) %>%
  distinct()

# ==============================================================================
# REPORTING
# ==============================================================================
print(as.data.frame(census_summary))

message("\nDetailed Audit:")
message(paste("Total Raw paths accounted for in Final Tree:", length(intersect(raw_paths, final_paths))))
message(paste("Raw paths that were RENAMED by Database (Path changed but species kept):", 
              length(intersect(raw_df$Species, final_df$Species)) - length(intersect(raw_paths, final_paths))))

# Save the lists for your records
write_csv(the_missing_85, "~/GitHub/achiasmy-synthesis/data/audit_missing_85_raw.csv")
write_csv(the_lost_67, "~/GitHub/achiasmy-synthesis/data/audit_lost_67_backbone.csv")

# ==============================================================================
# AUDIT: EXPLAINING THE GAP (5609 Rows vs 5271 Unique Paths)
# ==============================================================================

# 1. Generate the full path for every row in the original data
raw_with_paths <- raw_df %>%
  mutate(full_path = paste(Kingdom, Class, Order, Family, Genus, Species, sep="|"))

# 2. Identify which paths appear more than once (The Duplicates)
duplicated_paths <- raw_with_paths %>%
  group_by(full_path) %>%
  filter(n() > 1) %>%
  summarise(
    row_count = n(),
    original_species_names = paste(unique(Species), collapse = "; "),
    meiosis_types = paste(unique(`Meiosis Type`), collapse = "; "),
    diploid_nums  = paste(unique(`Diploid Number`), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(row_count))

# 3. Summary Report
message("--- Duplication Audit Results ---")
message(paste("Total Rows in Raw Data:   ", nrow(raw_df)))
message(paste("Unique Paths Generated:   ", length(unique(raw_with_paths$full_path))))
message(paste("Rows lost to Duplication: ", nrow(raw_df) - length(unique(raw_with_paths$full_path))))
message(paste("Number of Path Groups containing duplicates: ", nrow(duplicated_paths)))

# 4. View the top duplicates
print("Top duplicated taxa (Rows that collapsed into 1 path):")
print(head(duplicated_paths, 15))

# 5. Save the list of duplicates for manual review
write_csv(duplicated_paths, "~/GitHub/achiasmy-synthesis/data/audit_duplicates_collapsed.csv")

# ========================================================
# Audit: Identifying the 8 Initial Backbone Failures
# ========================================================

# 1. Reconstruct the full paths from the raw data
raw_paths_full <- raw_df %>%
  mutate(full_path = paste(Kingdom, Class, Order, Family, Genus, Species, sep="|"))

# 2. Identify paths that are in the raw data (5271 unique) 
#    but NOT in the backbone metadata (5263 unique)
missing_from_backbone <- raw_paths_full %>%
  filter(!full_path %in% backbone_df$tree_label) %>%
  select(Kingdom, Class, Order, Family, Genus, Species, full_path) %>%
  distinct()

# 3. Print the results
message(paste("Displaying the", nrow(missing_from_backbone), "taxa that failed the Backbone lookup:"))
print(as.data.frame(missing_from_backbone))

# 4. Save to CSV for your records
write_csv(missing_from_backbone, "~/GitHub/achiasmy-synthesis/data/audit_8_missing_from_backbone.csv")