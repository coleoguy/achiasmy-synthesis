# Terrence Sylvester
# pradakshanas@gmail.com
# 14 February 2022

# read in data
data.chroms <- read.csv("../data/mammal-TS-edit.csv", as.is = T)
data.achiasmy <- read.csv("../data/achiasmy-TS-ed.csv", as.is = T)

# get genus names for achiasmy data
data.achiasmy$Genus <- NA
for(i in 1:nrow(data.achiasmy)){
  data.achiasmy$Genus[i] <- unlist(strsplit(data.achiasmy$Species[i], split = " "))[1]
}

# if there is a range present then get a single chromosome number
for(i in 1:nrow(data.chroms)){
  if(data.chroms$female2n.range[i] != ""){
    hit <- as.numeric(unlist(strsplit(data.chroms$female2n.range[i], split = "-")))
    hit.seq <- seq(from = hit[1], to = hit[2], by = 2)
    data.chroms$female2n[i] <- sample(hit.seq,1)
  }
}

for(i in 1:nrow(data.achiasmy)){
  if(data.achiasmy$N.haploid.range[i] != ""){
    hit <- as.numeric(unlist(strsplit(data.achiasmy$N.haploid.range[i], split = "-")))
    hit.seq <- seq(from = hit[1], to = hit[2], by = 1)
    data.achiasmy$N..haploid.autosome.[i] <- sample(hit.seq,1)
  }
}

# masupial orders
masupials <- c("Didelphimorphia",
               "Paucituberculata",
               "Microbiotheria",
               "Dasyuromorphia",
               "Peramelemorphia",
               "Notoryctemorphia",
               "Diprotodontia")

# seperate chromosome dataset to eutherians and masupials
data.chroms.masupials <- data.chroms[(data.chroms$Order %in% masupials),]
data.chroms.eutherians <- data.chroms[!(data.chroms$Order %in% masupials),]

# get achiasmatic chordatas
data.achiasmy.chordats <- data.achiasmy[data.achiasmy$Phyla == "Chordata",]
data.chroms.achiasmy <- data.chroms.eutherians[data.chroms.eutherians$Genus %in% data.achiasmy.chordats$Genus,]

# get non achiasmatic eutherians
data.chiasmy.eutheria <- data.chroms.eutherians[!(data.chroms.eutherians$Genus %in% data.achiasmy.chordats$Genus),]
data.chiasmy.rodents <- data.chiasmy.eutheria[data.chiasmy.eutheria$Order == "Rodentia",]
data.chiasmy.Cricetidae <- data.chiasmy.rodents[data.chiasmy.rodents$Family == "Cricetidae",]

# plot
# chiasmatic eutherians vs achiasmatic chordates

# plot
plot(xlim = c(0,110),
     ylim = c(0,0.04),
     x = NULL,
     y = NULL,
     xlab = "Chromosome number",
     ylab = "Density",
     main = "Chiasmatic chordates vs achiasmatic chordates")
polygon(density(data.chiasmy.eutheria$female2n, na.rm = T),
        col = rgb(1,0,0,0.3),
        border = rgb(1,0,0,1),
        lwd = 2)
polygon(density(c(data.achiasmy.chordats$N..haploid.autosome.*2,data.chroms.achiasmy$female2n), na.rm = T),
        col = rgb(0,0,1,0.3),
        border = rgb(0,0,1,1),
        lwd = 2)
legend("topright",
       legend = c("Chiasmatic chordates", "Achiasmatic chordates"),
       inset = 0.01,
       fill = c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)),
       border = c(rgb(1,0,0,1),rgb(0,0,1,1)))

# plot
# chiasmatic rodents vs achiasmatic chordates

# plot
plot(xlim = c(0,110),
     ylim = c(0,0.04),
     x = NULL,
     y = NULL,
     xlab = "Chromosome number",
     ylab = "Density",
     main = "Chiasmatic Rodents vs achiasmatic chordates")
polygon(density(data.chiasmy.Cricetidae$female2n, na.rm = T),
        col = rgb(1,0,0,0.3),
        border = rgb(1,0,0,1),
        lwd = 2)
polygon(density(c(data.achiasmy.chordats$N..haploid.autosome.*2,data.chroms.achiasmy$female2n), na.rm = T),
        col = rgb(0,0,1,0.3),
        border = rgb(0,0,1,1),
        lwd = 2)
legend("topright",
       legend = c("Chiasmatic chordates", "Achiasmatic chordates"),
       inset = 0.01,
       fill = c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)),
       border = c(rgb(1,0,0,1),rgb(0,0,1,1)))


# plot
plot(xlim = c(0,110),
     ylim = c(0,0.04),
     x = NULL,
     y = NULL,
     xlab = "Chromosome number",
     ylab = "Density",
     main = "Chiasmatic Rodents vs achiasmatic chordates")
polygon(density(data.chiasmy.rodents$female2n, na.rm = T),
        col = rgb(1,0,0,0.3),
        border = rgb(1,0,0,1),
        lwd = 2)
polygon(density(c(data.achiasmy.chordats$N..haploid.autosome.[data.achiasmy.chordats$Genus == "Microtus"]*2,
                  data.chroms.achiasmy$female2n[data.chroms.achiasmy$Genus == "Microtus"]), na.rm = T),
        col = rgb(0,0,1,0.3),
        border = rgb(0,0,1,1),
        lwd = 2)
legend("topright",
       legend = c("Chiasmatic chordates", "Achiasmatic chordates"),
       inset = 0.01,
       fill = c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)),
       border = c(rgb(1,0,0,1),rgb(0,0,1,1)))














ach.microtus <- dat.ach[dat.ach$Genus == "Microtus",]
dat.ach.chords <- dat.ach[dat.ach$Phyla == "Chordata",]




# subset data
dat.microtus <- dat[dat$Genus == "Microtus",]
dat.rodents <- dat[c(dat$Genus != "Microtus" & dat$Order == "Rodentia"),]
dat.masupials <- dat[dat$Order %in% c("Didelphimorphia",
                                      "Paucituberculata",
                                      "Microbiotheria",
                                      "Dasyuromorphia",
                                      "Peramelemorphia",
                                      "Notoryctemorphia",
                                      "Diprotodontia"),]
dat.eutherian <- dat[!(dat$Order %in% c("Didelphimorphia",
                                        "Paucituberculata",
                                        "Microbiotheria",
                                        "Dasyuromorphia",
                                        "Peramelemorphia",
                                        "Notoryctemorphia",
                                        "Diprotodontia")),]
dat.eutherian.non.voles <- dat.eutherian[!(dat.eutherian$Genus %in% dat.ach.chords$Genus),]


microtus.chroms <- c(dat.microtus$female2n, ach.microtus$N..haploid.autosome.*2)
rodents.chroms <- dat.rodents$female2n
eutherians.chroms <- dat.eutherian$female2n
eutherian.non.voles.chroms <- dat.eutherian.non.voles$female2n
masupials.chroms <- dat.masupials$female2n

{
  par(mfcol = c(1,3))
  
  # microtus vs eutherians
  # plot
  plot(xlim = c(0,110),
       ylim = c(0,0.06),
       x = NULL,
       y = NULL,
       xlab = "Chromosome number",
       ylab = "Density",
       main = "Voles vs Eutherians")
  polygon(density(dat.ach.chords$N..haploid.autosome.*2, na.rm = T),
          col = rgb(0,0,1,0.3),
          border = rgb(0,0,1,1),
          lwd = 2)
  polygon(density(eutherian.non.voles.chroms, na.rm = T),
          col = rgb(1,0,0,0.3),
          border = rgb(1,0,0,1),
          lwd = 2)
  legend("topright",
         legend = c("Eutherians", "Voles"),
         inset = 0.01,
         fill = c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)),
         border = c(rgb(1,0,0,1),rgb(0,0,1,1)))
  
  # Voles vs Other rodents
  # plot
  plot(xlim = c(0,110),
       ylim = c(0,0.06),
       x = NULL,
       y = NULL,
       xlab = "Chromosome number",
       ylab = "density",
       main = "Voles vs Other rodents")
  polygon(density(microtus.chroms, na.rm = T),
          col = rgb(0,0,1,0.3),
          border = rgb(0,0,1,1),
          lwd = 2)
  polygon(density(rodents.chroms, na.rm = T),
          col = rgb(1,0,0,0.3),
          border = rgb(1,0,0,1),
          lwd = 2)
  legend("topright",
         legend = c("Rodents", "Voles"),
         inset = 0.01,
         fill = c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)),
         border = c(rgb(1,0,0,1),rgb(0,0,1,1)))
  
  # Voles vs Masupials
  # plot
  plot(xlim = c(0,70),
       ylim = c(0,0.15),
       x = NULL,
       y = NULL,
       xlab = "Chromosome number",
       ylab = "density",
       main = "Voles vs Masupials")
  polygon(density(microtus.chroms, na.rm = T),
          col = rgb(0,0,1,0.3),
          border = rgb(0,0,1,1),
          lwd = 2)
  polygon(density(masupials.chroms, na.rm = T),
          col = rgb(1,0,0,0.3),
          border = rgb(1,0,0,1),
          lwd = 2)
  legend("topright",
         legend = c("Masupials","Voles"),
         inset = 0.01,
         fill = c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)),
         border = c(rgb(1,0,0,1),rgb(0,0,1,1)))
  par(mfcol = c(1,1))
}

