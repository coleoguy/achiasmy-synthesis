# Terrence Sylvester
# pradakshanas@gmail.com
# 14 February 2022

# read in helper functions
source("helper.functions.R")

# read data
dat <- read.csv("../data/achiasmy.csv", as.is = T)

# get mammalia data only
dat <- dat[dat$Class == "Mammalia",]

# sample a single chromosome number for species that have a range of values
dat <- chromSampler(dat, info.column = 10,data.column = 10)
dat$Female.Diploid.. <- as.numeric(dat$Female.Diploid..)

# remove entries that do not have chromosome number for female diploid value
dat <- dat[!is.na(dat$Female.Diploid..),]

# get marsupial orders
marsupials <- c("Didelphimorphia",
               "Paucituberculata",
               "Microbiotheria",
               "Dasyuromorphia",
               "Peramelemorphia",
               "Notoryctemorphia",
               "Diprotodontia")
dat.marsupials <- dat[dat$Order %in% marsupials,]

# get eutherians
dat.eutherians <- dat[!(dat$Order %in% marsupials),] 
# remove monotremes from eutherian subset
dat.eutherians <- dat.eutherians[dat.eutherians$Order != "Monotremata",]
# get rodents
dat.rodents <- dat[dat$Order == "Rodentia",]
# get cricetidae
dat.cricetidae <- dat[dat$Family == "Cricetidae",]
# get microtus
dat$genera <- NA
for(i in 1:nrow(dat)){
  dat$genera[i] <- unlist(strsplit(dat$Species[i], split = " "))[1]
}
dat.microtus <- dat[dat$genera == "Microtus",]

par(mfrow = c(2,3))

#### plot ####
#### Eutherians (chiasmatic) \n vs \n Marsupials ####
plot(xlim = c(0,110),
     ylim = c(0,0.25),
     x = NULL,
     y = NULL,
     xlab = "Diploid chromosome number",
     ylab = "Density",
     main = "Eutherians (chiasmatic) \n vs \n Marsupials (achiasmatic)")
polygon(density(c(dat.marsupials$Female.Diploid..), na.rm = T),
        col = rgb(0,0,1,0.3),
        border = rgb(0,0,1,1),
        lwd = 2)
polygon(density(c(dat.eutherians$Female.Diploid..[dat.eutherians$Meiosis.Type == "chiasmatic"]), na.rm = T),
        col = rgb(1,0,0,0.3),
        border = rgb(1,0,0,1),
        lwd = 2)
length(dat.eutherians$Female.Diploid..[dat.eutherians$Meiosis.Type == "chiasmatic"])
legend("topright",
       legend = c("Eutherians (n=1397)", "Marsupials (n=38)"),
       inset = 0.01,
       fill = c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)),
       border = c(rgb(1,0,0,1),rgb(0,0,1,1)))

#### chiasmatic eutherians \n vs \n achiasmatic chordates ####
plot(xlim = c(0,110),
     ylim = c(0,0.06),
     x = NULL,
     y = NULL,
     xlab = "Diploid chromosome number",
     ylab = "Density",
     main = "Eutherians \n chiasmatic vs achiasmatic")
polygon(density(c(dat.eutherians$Female.Diploid..[dat.eutherians$Meiosis.Type != "chiasmatic"])),
        col = rgb(0,0,1,0.3),
        border = rgb(0,0,1,1),
        lwd = 2)
polygon(density(dat.eutherians$Female.Diploid..[dat.eutherians$Meiosis.Type == "chiasmatic"]),
        col = rgb(1,0,0,0.3),
        border = rgb(1,0,0,1),
        lwd = 2)
length(dat.eutherians$Female.Diploid..[dat.eutherians$Meiosis.Type == "chiasmatic"])
length(dat.eutherians$Female.Diploid..[dat.eutherians$Meiosis.Type != "chiasmatic"])
legend("topright",
       legend = c("Chiasmatic (n=1397)", "Achiasmatic (n=23)"),
       inset = 0.01,
       fill = c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)),
       border = c(rgb(1,0,0,1),rgb(0,0,1,1)))


#### chiasmatic rodents \n vs \n achiasmatic rodents ####
plot(xlim = c(0,110),
     ylim = c(0,0.06),
     x = NULL,
     y = NULL,
     xlab = "Diploid chromosome number",
     ylab = "Density",
     main = "Rodents \n chiasmatic vs achiasmatic")
polygon(density(c(dat.rodents$Female.Diploid..[dat.rodents$Meiosis.Type != "chiasmatic"])),
        col = rgb(0,0,1,0.3),
        border = rgb(0,0,1,1),
        lwd = 2)
polygon(density(dat.rodents$Female.Diploid..[dat.rodents$Meiosis.Type == "chiasmatic"]),
        col = rgb(1,0,0,0.3),
        border = rgb(1,0,0,1),
        lwd = 2)
length(dat.rodents$Female.Diploid..[dat.rodents$Meiosis.Type == "chiasmatic"])
length(dat.rodents$Female.Diploid..[dat.rodents$Meiosis.Type != "chiasmatic"])
legend("topright",
       legend = c("Chiasmatic (n=453)", "Achiasmatic (n=23)"),
       inset = 0.01,
       fill = c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)),
       border = c(rgb(1,0,0,1),rgb(0,0,1,1)))


#### Cricetidae (Chiasmatic) \n vs \n Cricetidae (aChiasmatic) ####
plot(xlim = c(0,110),
     ylim = c(0,0.1),
     x = NULL,
     y = NULL,
     xlab = "Diploid chromosome number",
     ylab = "Density",
     main = "Cricetidae \n chiasmatic vs achiasmatic")
polygon(density(c(dat.cricetidae$Female.Diploid..[dat.cricetidae$Meiosis.Type != "chiasmatic"])),
        col = rgb(0,0,1,0.3),
        border = rgb(0,0,1,1),
        lwd = 2)
polygon(density(dat.cricetidae$Female.Diploid..[dat.cricetidae$Meiosis.Type == "chiasmatic"]),
        col = rgb(1,0,0,0.3),
        border = rgb(1,0,0,1),
        lwd = 2)
length(dat.cricetidae$Female.Diploid..[dat.cricetidae$Meiosis.Type == "chiasmatic"])
length(dat.cricetidae$Female.Diploid..[dat.cricetidae$Meiosis.Type != "chiasmatic"])
legend("topright",
       legend = c("Chiasmatic (n=193)", "Achiasmatic (n=17)"),
       inset = 0.01,
       fill = c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)),
       border = c(rgb(1,0,0,1),rgb(0,0,1,1)))



#### Microtus (Chiasmatic) \n vs \n Microtus (aChiasmatic) ####
plot(xlim = c(0,110),
     ylim = c(0,0.1),
     x = NULL,
     y = NULL,
     xlab = "Diploid chromosome number",
     ylab = "Density",
     main = "Microtus \n chiasmatic vs achiasmatic")
polygon(density(c(dat.microtus$Female.Diploid..[dat.microtus$Meiosis.Type != "chiasmatic"])),
        col = rgb(0,0,1,0.3),
        border = rgb(0,0,1,1),
        lwd = 2)
polygon(density(dat.microtus$Female.Diploid..[dat.microtus$Meiosis.Type == "chiasmatic"]),
        col = rgb(1,0,0,0.3),
        border = rgb(1,0,0,1),
        lwd = 2)
length(dat.microtus$Female.Diploid..[dat.microtus$Meiosis.Type == "chiasmatic"])
length(dat.microtus$Female.Diploid..[dat.microtus$Meiosis.Type != "chiasmatic"])
legend("topright",
       legend = c("Chiasmatic (n=36)", "Achiasmatic (n=16)"),
       inset = 0.01,
       fill = c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)),
       border = c(rgb(1,0,0,1),rgb(0,0,1,1)))