# Sean Chien 
# Checking data 
dat <- read.csv("Achiasmy_full_meiosis_data.csv")

# remove space before and after genus and species
cols <- 1:11
char_cols <- sapply(dat[, cols], is.character)
dat[, cols][, char_cols] <- lapply(dat[, cols][, char_cols], trimws)

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
<<<<<<< HEAD
# write.csv(diff_results, "../data/different_rows.csv", row.names = FALSE)
=======
write.csv(diff_results, "different_rows.csv", row.names = FALSE)
>>>>>>> parent of 4b40860 (script)

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
# dat_clean <- dat[!duplicated(dat), ]
# write.csv(dat_clean, "../data/Achiasmy_full_meiosis_data.csv", row.names = FALSE)


#### Type 1 ####
need.check.count <- 0
for (i in 1:length(all.sp)) {
  idx <- which(dat$Sci_name == all.sp[i])
  if (length(idx) > 1) {
    subdat <- dat[idx, ]
    # check columns 1–6 are identical
    cols_1_6_same <- all(apply(subdat[, c(1:6,11)], 2, function(x) length(unique(x)) == 1))
    # check columns 7–11: some rows have data, some are empty
    cols_7_11_issue <- any(apply(subdat[, 7:10], 2, function(x) {
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

# merge #
# the data with same taxonomic level and citation 
# but with some empty entry
# fill with those NAs 
merged_list <- list()
for (i in 1:length(all.sp)) {
  idx <- which(dat$Sci_name == all.sp[i])
  subdat <- dat[idx, ]
  if (nrow(subdat) == 1) {
    merged_list[[i]] <- subdat
  } else {
    new_row <- subdat[1, ]  # start with first row
    
    # columns 7–10: fill with non-empty values
    for (col in 7:10) {
      vals <- subdat[, col]
      
      # remove NA and empty
      vals_clean <- vals[vals != "" & !is.na(vals)]
      
      if (length(vals_clean) > 0) {
        new_row[, col] <- vals_clean[1]  # take first available value
      } else {
        new_row[, col] <- ""  # keep empty if nothing exists
      }
    }
    
    merged_list[[i]] <- new_row
  }
}

# combine back into a data frame
dat_merged <- do.call(rbind, merged_list)
# reset rownames
rownames(dat_merged) <- NULL
dat <- dat_merged[,-13]
# wrtie new file
# write.csv(dat, "../data/Achiasmy_full_meiosis_data.csv", row.names = FALSE)



### Type 2 ###
# data without conflict
# either with data or NA or the data are the same
# multiple citations

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


### Type 3 ###
# data with conflict
need.check.count <- 0

for (i in 1:length(all.sp)) {
  idx <- which(dat$Sci_name == all.sp[i])
  
  if (length(idx) > 1) {
    subdat <- dat[idx, ]
    
    # columns 1–6 identical
    cols_1_6_same <- all(apply(subdat[, 1:6], 2, function(x) length(unique(x)) == 1))
    
    # columns 7–10: detect conflict (more than one real value)
    cols_7_10_conflict <- any(apply(subdat[, 7:10], 2, function(x) {
      vals <- x[x != "" & !is.na(x)]  # keep only real data
      length(unique(vals)) > 1        # TRUE if conflict exists
    }))
    
    if (cols_1_6_same & cols_7_10_conflict) {
      print("===============================================================")
      print(i)
      need.check.count <- need.check.count + 1
      print(subdat)
      print("===============================================================")
    }
  }
}