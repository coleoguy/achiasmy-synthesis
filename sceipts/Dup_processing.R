# Sean Chien 
library(writexl)
library(readxl)
# Checking data 
# dat <- read.csv("../data/Achiasmy_full_meiosis_data.csv")

# remove space before and after genus and species
cols <- 1:6
char_cols <- sapply(dat[, cols], is.character)
dat[, cols][, char_cols] <- lapply(dat[, cols][, char_cols], trimws)
dat <- dat[,-13]
write_xlsx(dat, "../data/Achiasmy_full_meiosis_data.xlsx")

dat <- read_xlsx('../data/Achiasmy_full_meiosis_data.xlsx')
dat$Sci_name <- paste(dat$Genus,dat$Species,sep = '_')
all.sp <- unique(paste(dat$Genus,dat$Species,sep = '_'))
#### checking taxonomic level ####
# if data entry that have same scientific name have same Kindon Class Order...
diff_results <- data.frame()
for (i in 1:length(all.sp)){
  if (length(which(dat$Sci_name == all.sp[i]))>1){
    
    # Check if all rows are identical
    subset_rows <- dat[dat$Sci_name == all.sp[i], ][,c(1:6)]
    all_same <- all(apply(subset_rows, 2, function(col) length(unique(col)) == 1))

    if (!all_same) {
      print("===============================================================")
      print(subset_rows)
      print("===============================================================")
      diff_results <- rbind(diff_results, subset_rows)
    }
  }
}

# write.csv(diff_results, "../data/different_rows.csv", row.names = FALSE)

#### Duplication ####
dup <- 0
for (i in 1:length(all.sp)){
  idx <- which(dat$Sci_name == all.sp[i])
  if (length(idx) > 1) {
    subdat <- dat[idx, ]
    dup_rows <- duplicated(subdat) | duplicated(subdat, fromLast = TRUE)
    if (any(dup_rows)) {
      print("===============================================================")
      dup <- dup + 1
      print(i)
      print(all.sp[i])
      print("Identical rows:")
      print(subdat[dup_rows, ])
      print("===============================================================")
    }
  }
}
dat_clean <- dat[!duplicated(dat), ]
dat <- dat_clean[,-13]
# write_xlsx(dat, "../data/Achiasmy_full_meiosis_data.xlsx")


#### Type 1 #### 
# combining data 
# multiple citation contribute to
dat <- read_xlsx('../data/Achiasmy_full_meiosis_data.xlsx')
dat$Sci_name <- paste(dat$Genus,dat$Species,sep = '_')
all.sp <- unique(paste(dat$Genus,dat$Species,sep = '_'))
dat[, 7:11] <- apply(dat[, 7:11], 2, trimws)
need.check.count <- 0
for (i in 1:length(all.sp)) {
  idx <- which(dat$Sci_name == all.sp[i])
  
  if (length(idx) > 1) {
    
    subdat <- dat[idx, ]
    
    # check columns 1–6 are identical
    cols_1_6_same <- all(apply(subdat[, c(1:6,12)], 2, function(x) length(unique(x)) == 1))
    
    # check columns 7–11: some rows have data, some are empty
    cols_7_11_issue <- any(apply(subdat[, 7:11], 2, function(x) {
      has_data <- x != "" & !is.na(x)
      any(has_data) & any(!has_data)
    }))
    
    if (cols_1_6_same & cols_7_11_issue) {
      print("===============================================================")
      print(i)
      need.check.count <- need.check.count + 1
      print(subdat)
      print("===============================================================")
    }
  }
}

merged_list <- list()

for (i in 1:length(all.sp)) {
  idx <- which(dat$Sci_name == all.sp[i])
  subdat <- dat[idx, ]
  
  if (nrow(subdat) == 1) {
    merged_list[[i]] <- subdat
  } else {
    
    # check columns 1–6 identical
    cols_1_6_same <- all(apply(subdat[, 1:6], 2, function(x) length(unique(x)) == 1))
    
    if (!cols_1_6_same) {
      # do not merge if base columns differ
      merged_list[[i]] <- subdat
      next
    }
    
    new_row <- subdat[1, ]
    conflict_flag <- FALSE
    
    for (col in 7:10) {
      vals <- subdat[, col]
      vals_clean <- vals[vals != "" & !is.na(vals)]
      unique_vals <- unique(vals_clean)
      
      if (length(unique_vals) == 1) {
        # safe → fill value
        new_row[, col] <- unique_vals[1]
        
      } else if (length(unique_vals) == 0) {
        # all NA/empty
        new_row[, col] <- ""
        
      } else {
        # TRUE conflict (e.g., "Xyp" vs "Xyyp")
        conflict_flag <- TRUE
        break
      }
    }
    
    if (conflict_flag) {
      # NOT merge → keep all rows
      merged_list[[i]] <- subdat
    } else {
      # safe → merged row
      merged_list[[i]] <- new_row
    }
  }
}

dat_merged <- do.call(rbind, merged_list)
rownames(dat_merged) <- NULL
dat <- dat_merged[,-13]
# write_xlsx(dat, "../data/Achiasmy_full_meiosis_data.xlsx")


#### Type 2 ####
# data without conflict 
# either with data or NA or the data are the same 
# multiple citations
dat <- read_xlsx('../data/Achiasmy_full_meiosis_data.xlsx')
dat$Sci_name <- paste(dat$Genus,dat$Species,sep = '_')
all.sp <- unique(paste(dat$Genus,dat$Species,sep = '_'))

need.check.count <- 0
for (i in 1:length(all.sp)) {
  idx <- which(dat$Sci_name == all.sp[i])
  
  if (length(idx) > 1) {
    subdat <- dat[idx, ]
    
    # columns 1–6 identical
    cols_1_6_same <- all(apply(subdat[, 1:6], 2, function(x) length(unique(x)) == 1))
    
    # columns 7–10: allow same value OR NA/"" mixed in, but no conflicts
    cols_7_10_ok <- all(apply(subdat[, 7:10], 2, function(x) {
      vals <- x[x != "" & !is.na(x)]  # keep only real data
      length(unique(vals)) <= 1       # at most one unique value
    }))
    
    # column 11 must be different
    col11_diff <- length(unique(subdat[, 11])) > 1
    
    if (cols_1_6_same & cols_7_10_ok & col11_diff) {
      print("===============================================================")
      print(i)
      need.check.count <- need.check.count + 1
      print(subdat)
      print("===============================================================")
    }
  }
}


#### Type 3 #### 
# data with conflict 
need.check.count <- 0
conflict_list <- list()

for (i in 1:length(all.sp)) {
  idx <- which(dat$Sci_name == all.sp[i])
  
  if (length(idx) > 1) {
    subdat <- dat[idx, ]
    
    # columns 1–6 identical
    cols_1_6_same <- all(apply(subdat[, 1:6], 2, function(x) length(unique(x)) == 1))
    
    # columns 7–10: detect conflict
    cols_7_10_conflict <- any(apply(subdat[, 7:10], 2, function(x) {
      vals <- x[x != "" & !is.na(x)]
      length(unique(vals)) > 1
    }))
    
    if (cols_1_6_same & cols_7_10_conflict) {
      print("===============================================================")
      print(i)
      need.check.count <- need.check.count + 1
      print(subdat)
      print("===============================================================")
      
      # 👉 ADD THIS: save to list
      subdat$group_id <- need.check.count
      conflict_list[[need.check.count]] <- subdat
    }
  }
}

# combine after loop
conflict_dat <- do.call(rbind, conflict_list)
rownames(conflict_dat) <- NULL

write.xlsx(conflict_dat, "Type3_conflicts.xlsx")