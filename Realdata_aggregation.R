rm(list=ls())
set.seed(1)

library(dplyr)
library(ggplot2)
library(viridis)

####### Load air temperature data
load("tmp_monthly_2015-2100.RData") # Original temperature data
load(file="residuals_given_month_standardized_arma.RData") # Residuals

locations <- mydata$locs # Location information

lons <- unique(locations[,1]) # longitude
lats <- unique(locations[,2]) # latitude

p=2 # Aggregation level (we use squares of size p x p, see section 4.4)
length(lats) / p
length(lons) / p

temp <- c()
K <- seq(1, length(lons), by=p)
J <- seq(1, length(lats), by=p)
for(k in K){
  for(j in J){
    for(i in k:(k+p-1)){
      temp <- c(temp, ((i-1)*length(lats)+j):(((i-1)*length(lats)+j) + (p-1)))
    }
  }
}

batch <- (p*p)
B <- dim(locations)[1] / batch
locations.agg <- matrix(NA, B, 2)
data.agg <- matrix(NA, dim(resi.given.month)[1], B)

for(b in 1:B){
  locations.agg[b,] <- colMeans(locations[temp[((b-1)*batch+1):(b*batch)],])   
  data.agg[,b] <- rowMeans(resi.given.month[,temp[((b-1)*batch+1):(b*batch)]])
}

locs <- locations.agg
N <- dim(locs)[1] # Number of locations for aggregated data
data.for.test <- data.agg

## Heat map
tmp <- data.for.test[1,] # January in 2015
dt <- data.frame(lon=locs[,1], lat=locs[,2], tmp=tmp)

gg <- ggplot(dt, aes(x=lon, y=lat, fill=tmp))
gg <- gg + geom_tile(color="white", linewidth=0.1, show.legend = F)
gg <- gg + scale_fill_viridis(name="Residuals",option="inferno",
                              limits = c(-3.5, 3.5), breaks = seq(-2,2, by=1))
gg <- gg + labs(x="Longitude", y="Latitude")

gg <- gg + theme_bw() +
  coord_cartesian(ylim = c(-90,90)) +
  scale_y_continuous(breaks = seq(-80,80, by=20)) +
  # scale_y_continuous(breaks = seq(-80,80, by=20), position = "right") +
  scale_x_continuous(breaks = seq(0,360, by=50))

gg <- gg + theme(plot.title=element_text(hjust=0))
gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + theme(axis.text=element_text(size=10))

My_Theme <- theme(
  axis.title.x = element_text(size = 25),
  axis.title.y = element_text(size = 25),
  axis.text.x = element_text(size = 22),
  axis.text.y = element_text(size = 22),
  legend.text = element_text(size=15),
  legend.title = element_text(size=22))

gg + My_Theme # Heat map for aggregated residuals

lon <- locs[,1]
lat <- locs[,2]
lon <- (lon*pi) / 180
lat <- (lat*pi) / 180

R <- 6371 # radius of earth 

x = R*cos(lat)*cos(lon)
y = R*cos(lat)*sin(lon)
z = R*sin(lat)

xyz <- data.frame(x=x, y=y, z=z)

dist.mat <- as.matrix(dist(xyz)) ## Chordal distance matrix

