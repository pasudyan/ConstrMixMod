rm(list = ls())

library(sp)
library(leaflet)
library(magrittr)
library(ggplot2)
library(ggmap)

# setwd("/Volumes/PCCR_Backup/Ayu /Research/Density_Project/data_aug_reject_smpl/flow_cyt/crime_app/data")
setwd("~/Documents/Research/Density_Project/data_aug_reject_smpl/flow_cyt/crime_app/data")

options(digits = 15)

# --------------- Read in data for neighborhood coordinates --------------- #
neigh_coord_df <- read.csv("CommAreas.csv", header = T)
neigh_coord    <- gsub("MULTIPOLYGON", "", neigh_coord_df$the_geom)
neigh_coord    <- chartr('()', '  ', neigh_coord)
num_neigh      <- 77

douglas <- neigh_coord[1]
a <- strsplit(douglas, ",")[[1]]
a <- gsub("  ", "", a, fixed = T)
aa <- gsub(" -", "-" , a, fixed = T)
aa <- strsplit(aa, " ")

matrix(unlist(lapply(aa, function(x){ as.numeric(x)})), ncol = 2, byrow = T)

coord_list <- rep(list(list()), nrow(neigh_coord_df))
for (i in 1:nrow(neigh_coord_df)){
  print(i)
  neigh_coord_df[i, 1] %>% gsub("MULTIPOLYGON", "", ., fixed = T) %>% 
    chartr('()', '  ', .) -> temp
    strsplit(temp, ",")[[1]] %>%
    gsub("   ", "", ., fixed = T) %>%
    gsub(" -", "-" , ., fixed = T) -> temp 
    a <- strsplit(temp, " ")
    coord_list[[i]] <- matrix(unlist(lapply(a, function(x){as.numeric(x)})), ncol = 2, byrow = T)
}

# -------------- Turn neighborhood into Polygon --------------- #
# Summary of the coordinates
summ_coord <- lapply(coord_list, function(x) summary(x))

fix_list   <- coord_list[[36]]
fix_list_a <- fix_list[1329:nrow(fix_list), ]
fix_list_b <- fix_list[1:1328, ]

df_fix <- data.frame(fix_list_a[-380, 2])
df_fix <- cbind(df_fix, fix_list_a[-1, 1])

colnames(df_fix) <- c("Longitude", "Latitude")
colnames(fix_list_b) <- c("Longitude", "Latitude")
fixed_mat <- rbind(fix_list_b, df_fix)

coord_list[[36]] <- as.matrix(fixed_mat)

fix_list   <- coord_list[[75]]
which(is.na(fix_list[,1]) == TRUE)
which(is.na(fix_list[,2]) == TRUE)

fix_list_a <- fix_list[2228:2317, ]

df_fix <- data.frame(fix_list_a[-90, 2])
df_fix <- cbind(df_fix, fix_list_a[-1, 1])

colnames(df_fix) <- c("Longitude", "Latitude")

fix_list <- fix_list[-c(2228, 2317), ]

fix_list_b <- fix_list[1:2227, ]
fix_list_c <- fix_list[2316:nrow(fix_list), ]

colnames(fix_list_b) <- c("Longitude", "Latitude")
colnames(fix_list_c) <- c("Longitude", "Latitude")

fixed_mat  <- rbind(fix_list_b, df_fix)
colnames(fixed_mat) <- c("Longitude", "Latitude")

fixed_mat  <- rbind(fixed_mat, fix_list_c)

coord_list[[75]] <- as.matrix(fixed_mat)

coord_list <- lapply(1:num_neigh, function(x) coord_list[[x]][, c(2,1)])

min_coord <- matrix(unlist(
  lapply(1:num_neigh, function(x) c(min(coord_list[[x]][,1]), min(coord_list[[x]][,2])))),
  ncol = 2, byrow = T)
min_all_coord   <- apply(min_coord, 2, min)

max_coord <- matrix(unlist(
  lapply(1:num_neigh, function(x) c(max(coord_list[[x]][,1]), max(coord_list[[x]][,2])))),
  ncol = 2, byrow = T)
max_all_coord   <- apply(max_coord, 2, max)

save(coord_list, min_all_coord, max_all_coord, file = "Neigh_coord_list.RData")

# =============== Attempt to scale coordinates ================ #

coord_list_sc <- lapply(1:num_neigh, function(x)
  t(t(coord_list[[x]]) - min_all_coord))
lapply(1:num_neigh, function(x) summary(coord_list[[x]]))

min_coord_sc <- matrix(unlist(
  lapply(1:num_neigh, function(x) c(min(coord_list_sc[[x]][,1]), min(coord_list_sc[[x]][,2])))),
  ncol = 2, byrow = T)
min_all_coord_sc <- apply(min_coord_sc, 2, min)

max_coord_sc <- matrix(unlist(
  lapply(1:num_neigh, function(x) c(max(coord_list_sc[[x]][,1]), max(coord_list_sc[[x]][,2])))),
  ncol = 2, byrow = T)
max_all_coord_sc <- apply(max_coord_sc, 2, max)

coord_list_sc_t <- lapply(1:num_neigh, function(x) coord_list_sc[[x]][, c(2,1)])

# ================ Turning list into Polygons - original ================= #
num_neigh <- length(coord_list_sc_t)
l_poly <- lapply(1:num_neigh, function(x) {
  temp <- Polygon(coord_list_sc_t[[x]])
  Polygons(list(temp), neigh_coord_df$COMMUNITY[x])
})
SpPol <- SpatialPolygons(l_poly)
plot(SpPol, pbg = "red")
abline(v = 0, col = "red")
abline(v = 0.42, col = "red")
abline(h = 0, col = "red")
abline(h = 0.38, col = "red")

tt = matrix(0, nrow = 10000, ncol = 2)
tt[,1] <- seq(0, .42, length = 100)
tt[,2] s<- rep(seq(0, .38, length = 100), each = 100)

a <- constr(t(tt))
df_a <- data.frame(x = tt[,1], y = tt[,2], z = a)
ggplot(df_a) + geom_point(aes(x = y, y = x, col = z))

save(SpPol, l_poly, max_all_coord_sc, file = "chicago_poly.RData")

# ================ Turning list into Polygons - unit square ================= #
# Scaling polygon into a unit square distance instead of a half square
coord_list_sc_t_2 <- lapply(1:num_neigh, function(x) coord_list_sc_t[[x]]*2)

num_neigh <- length(coord_list_sc_t_2)
l_poly <- lapply(1:num_neigh, function(x) {
  temp <- Polygon(coord_list_sc_t_2[[x]])
  Polygons(list(temp), neigh_coord_df$COMMUNITY[x])
})
SpPol <- SpatialPolygons(l_poly)
plot(SpPol, pbg = "red")
abline(v = 0, col = "red")
abline(v = 0.42*2, col = "red")
abline(h = 0, col = "red")
abline(h = 0.38*2, col = "red")

max_all_coord_sc <- max_all_coord_sc*2

save(SpPol, l_poly, max_all_coord_sc, file = "chicago_poly_sc.RData")

# --------------- Read in data for crime events --------------- #
crime  <- read.csv("Crimes_-_2001_to_present_-_Map.csv") 
murder <- subset(crime, Primary.Type == "HOMICIDE")
murder <- subset(murder, Latitude > 40)
murder_set <- subset(murder, Year >= 2012 & Year <= 2017)
qmplot(Longitude, Latitude, data = murder_set, map_type = "toner-lite", colour = I('red'), 
       size = I(3), darken = 0.3)

murder_coord <- data.frame(Longitude = murder_set$Longitude,
                           Latitude = murder_set$Latitude)

murder_loc <- murder_set$Location 
murder_loc <- chartr('()', '  ', murder_loc)
murder_loc <- gsub(" ", "", murder_loc)
murder_loc <- strsplit(murder_loc, ",")

murder_loc <- matrix(unlist(lapply(murder_loc, function(x) as.numeric(x))), ncol = 2, byrow = T)

colnames(murder_loc) <- c("latitude", "longitude")

# =================== Format Murder Data ================== # 
murder_loc_sc <- t(t(murder_loc) - min_all_coord)
colnames(murder_loc_sc) <- c("latitude", "longitude")

save(murder_loc_sc, murder_loc, min_all_coord, file = "chi_murder_2012_2017.RData")

# ============ Check if Murder is in any neighborhood ============ #
# Naive way of doing things
in_obs <- lapply(1:num_neigh, function(x){
  point.in.polygon(point.x = murder_loc[, "latitude"], 
                   point.y = murder_loc[, "longitude"], 
                   pol.x = l_poly[[x]]@Polygons[[1]]@coords[, 2],
                   pol.y = l_poly[[x]]@Polygons[[1]]@coords[, 1])
})

in_obs_neigh <- unlist(lapply(in_obs, function(x) sum(x)))
max_obs_neigh <- which.max(in_obs_neigh)
neigh_coord_df[max_obs_neigh, ]
# Most homicide happened in neighborhood Austin 

all_in_obs <- Reduce("+", in_obs)
sum_in_obs <- sum(all_in_obs)
print(sum_in_obs)

# Function to determine if point is in polygon
cnstr_chicago  <- function(XX, l_poly){
  num_neigh <- 77
  in_obs <- lapply(1:num_neigh, function(x){
    point.in.polygon(point.x = XX$Latitude, 
                     point.y = XX$Longitude, 
                     pol.x = l_poly[[x]]@Polygons[[1]]@coords[, 2],
                     pol.y = l_poly[[x]]@Polygons[[1]]@coords[, 1])
  })
  in_obs <- Reduce("+", in_obs)
  sum(in_obs)
}

