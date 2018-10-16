library(ggplot2)
library(ggmap)
library(leaflet)
library(sp)
numeric(length = 15)

setwd("/Volumes/PCCR_Backup/Ayu /Research/Density_Project/data_aug_reject_smpl/flow_cyt/crime_app/data")

crime  <- read.csv("Crimes_-_2001_to_present_-_Map.csv", header = TRUE) 
murder <- subset(crime, Primary.Type == "HOMICIDE")
murder <- subset(murder, Latitude > 40)
save(murder, file = "chi_murder.RData")

murder_loc <- murder$Location 
murder_loc <- chartr('()', '  ', murder_loc)
murder_loc <- gsub(" ", "", murder_loc)
murder_loc <- strsplit(murder_loc, ",")
murder_loc <- as.double(unlist(murder_loc))

qmplot(Longitude, Latitude, data = murder, colour = I('red'), 
       size = I(3), darken = 0.3)





