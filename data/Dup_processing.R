dat <- read.csv("Achiasmy_full_meiosis_data.csv")
dat$Sci_name <- paste(dat$Genus,dat$Species,sep = '_')
all.sp <- unique(paste(dat$Genus,dat$Species,sep = '_'))

need.check.count <- 0
for (i in 1:length(all.sp)){
  if (length(which(dat$Sci_name == all.sp[i]))>1){
    need.check.count <- need.check.count + 1
    print("===============================================================")
    print(all.sp[i])
    print(dat[which(dat$Sci_name == all.sp[i]),])
    print("===============================================================")
  }
}
