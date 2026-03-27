# Sean Chien 
# Checking data 
dat <- read.csv("../data/Achiasmy_full_meiosis_data.csv")

# remove space before and after genus and species
cols <- 1:6
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
write.csv(diff_results, "../data/different_rows.csv", row.names = FALSE)

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



# combining data
# multiple citation contribute to 
need.check.count <- 0
for (i in 1:length(all.sp)){
  if (length(which(dat$Sci_name == all.sp[i]))>1){
    need.check.count <- need.check.count + 1
    print("===============================================================")
    print(i)
    print(all.sp[i])
    
    # type 1
    # combine data 
    
    
    print(dat[which(dat$Sci_name == all.sp[i]),])
    print("===============================================================")
  }
}




#### Type 1####


# combining data
# multiple citation contribute to 
need.check.count <- 0
for (i in 1:length(all.sp)){
  if (length(which(dat$Sci_name == all.sp[i]))>1){
    need.check.count <- need.check.count + 1
    print("===============================================================")
    print(i)
    print(all.sp[i])
    
    # type 1
    # combine data 
    
    
    print(dat[which(dat$Sci_name == all.sp[i]),])
    print("===============================================================")
  }
}



