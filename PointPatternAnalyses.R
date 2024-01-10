###Spatial point pattern analyses 
##UNM herbarium spatial patterns and gaps project
##code written by eml on March21 2023
#https://mgimond.github.io/Spatial/chp11_0.html

#load libraries
library(sf)
library(stars)
library(tidyr)
library(dplyr)
library(spatstat)
library(AOI)

load(file="script4_envData_complete.RData")
load(file="Data_both-allandgeo.RData")


try<-st_as_sf(x=env.v,
              coords = c("decimalLongitude", "decimalLatitude"),
              crs= "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

dat.sf = st_transform(try, 32113)
#dats = dat.sf %>% dplyr::select(rowID, Accepted_family, geometry)
#dat.ppp = as.ppp(env.v[, c("decimalLongitude", "decimalLatitude")], W=owin(aoi.ppp)) #these are the useful coordinates with the species as marks



#Define area of interest
aoi=aoi_get(state="NM") %>% st_transform(crs=4326)
class(aoi)
aoi = st_transform(aoi, 32113)
aoi.ppp =st_geometry(aoi)
aoi.owin = as.owin(aoi.ppp) #this is the useful polygon

# nm_owin = maptools::as.ppp.SpatialPointsDataFrame(as_Spatial(aoi))
# nm_owin


#check for aoi and coordinate alignment
plot(st_geometry(aoi))
plot(st_geometry(dat.sf), add = TRUE, col = "lightblue") #looks good. things line up



###START HERE####

# #####created marked point process for collections
#dc <- read.table("dat_cl_ppp.txt", header=TRUE)

dat.sf<- st_as_sf(x=env.v,
                 coords = c("decimalLongitude", "decimalLatitude"),
                 crs= "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

dat.sf = st_transform(dat.sf, 32113)
dats = st_geometry(dat.sf)

dat.ppp = as.ppp(dats) #these are the useful coordinates


#ppp <- as.ppp(dat.sf, window=aoi.owin)




###check your spatstat objects
class(dat.ppp)
class(aoi.owin)

#link the window and point pattern
ppp = dat.ppp[aoi.owin] #there are some records being dropped out of the analysis here. 

plot(ppp, pch='.') #good! `ppp` is the important bit now



###quadrat count
q = quadratcount(ppp, nx = 10, ny = 10)
plot(q)

quadrat.test(ppp, nx = 10, ny = 10) #data do not conform to CSR using quadrat counts


###check kernel density 
den = density(ppp)
den

plot(den, main = "")
#plot(ppp, pch = '.', cex = 0.05, alpha=.3, add = TRUE)
contour(den, add = TRUE)


###First order spatial analyses: density based analyses of collections in NM

####Second order spatial analyses: distance based analyses of collections in NM

#calculate mean distance to the first nearest neighbor for each collection point

nn.all <- nndist(ppp, k=1)

##Add columns to env.v with ANN for nearest five neighbors (ANN) and the nearest neighbor (NN[,1])
env.v <- env.v %>%
  mutate(NN = nndist(dat.ppp, k=1)) 

env.v <- env.v %>%
  mutate(ANN = rowSums(nndist(dat.ppp, k=1:5))/5)

##save rda
save(env.v, file="env.vPointPattern.RData")


#check nearest neighbor calculations
nn.mean <-mean(nndist(ppp, k=1)) #89.02779 meters? check units

mean(env.v$ANN) #average distance to the nearest five neighbors = 0.2283501 meters (is this right??)
median(env.v$ANN)

mean(env.v$NN[,1]) #this is about right. ~99 meters to the nearest neighbor
mean(nn.all)
median(env.v$NN[,1]) #median is 0!!


##Write the environment and data env.v with point processes in columns
#save new CSV with point pattern in env.v
write.csv(env.v, "/Users/elizabethlombardi/Desktop/Research/UNM/Herbarium collections patterns and gaps project/env.vFull_ppp.csv")


##save rda
save(env.v, file=" ")


  

#ANN
#https://mgimond.github.io/Spatial/point-pattern-analysis-in-r.html#distance-based-analysis-1
ANN <- apply(nndist(ppp, k=1:100),2,FUN=mean) #calculate average nn for first 100 closest neighbors
plot(ANN ~ eval(1:100), type="b", main=NULL, las=1)


#Kest function
K <- Kest(ppp, nlarge=1000) #bump up N to match dataframe size eventually
plot(K, main=NULL, las=1, legendargs=list(cex=.5, xpd=TRUE, inset=c(1.01, 0) ))



####Hypothesis testing
#H0:all collections from NM have been made at random locations throughout the state
#H1: collections made from NM are more clustered than CSR

nn.mean <- mean(nndist(ppp, k=1))

#Null model! Generate homogenous point process using Monte Carlo methods 
n     <- 1000L             # Number of simulations
nn.rand <- vector(length = n) # Create an empty object to be used to store simulated ANN values
for (i in 1:n){  
  rand.p   <- rpoint(n=ppp$n, win=aoi.owin)
  nn.rand[i] <- mean(nndist(rand.p, k=1))  # Tally the ANN values
}
plot(rand.p, pch=16, cex=.1, main=NULL, cols=rgb(0,0,0,0.5)) # randomly distributed points. Equal number to real data.


##Histogram of expected (CSR) versus observed distribution
hist(nn.rand, main=NULL, las=1, breaks=1000, col="#999933", xlim=range(nn.mean, 1000))
hist(nn.rand, main=NULL, las=1, breaks=10, col="#999933", xlim=range(0, 700), ylim=range(0, 50))
abline(v=nn.mean, col="#AA4499") #the observed nearest neighbor distance is WAY smaller than expected under CSR

mean(nn.rand) #mean NN if random
nn.mean #observed NN


##Pseudo-pvalue for MC distributions
N.greater <- sum(nn.rand > nn.mean) #get the percent of simulated nn values that are greater than observed NN (all of them...)

p <- min(N.greater + 1, n + 1 - N.greater) / (n +1)
p #p-value 0.0099 REPORT THIS



###
###
####NEXT: add elevation or other rasters as covariates that might explain the clustered point patterns
https://geobgu.xyz/r/point-pattern-analysis.html


####Determine colors to match the selected TOL muted palette
## Show the colour palette
library(khroma)
plot_scheme(colour("muted")(9), colours=TRUE)

####Reproject the dat.sf to match the raster projectsion
dat.sf <- st_transform(dat.sf, crs = "+proj=tmerc +lat_0=31 +lon_0=-106.25 +k=0.9999 +x_0=500000 +y_0=0 +datum=NAD83 +units=m +no_defs")


#Prepare covariate data
elevation <- get_elev_raster(aoi, z = 9)
elev <- read_stars("/Users/elizabethlombardi/Desktop/Research/UNM/Herbarium collections patterns and gaps project/Environmental Data/nm_relief_color_tif/nm_relief_color.tif")


elev = elevation

crop_elev <- crop(elev, aoi)
mask_elev <- mask(crop_elev, mask=crop_elev)

par(mar=c(1,1,1,1))
plot(crop_elev, main = "NM Elevation Map")

quantile(mask_elev)
elev_zones <- reclassify(crop_elev,
                         c(0, 865, 1,
                           865, 1580, 2,
                           1580, 2298, 3,
                           2298, 3985, 4))
elev_zones <- ratify(elev_zones)
plot(elev_zones, main = "Elevation Zones")

#Convert to a Spatstat-compatible object
elev_zones <- as.im.RasterLayer(elev_zones)

#Tesselate the image
tes_elev <- tess(image = elev_zones)
plot(tes_elev, main = "Tesselated Elevation Zones")
plot(st_geometry(dat.sf), add = TRUE, cex=0.1, col = "white")
plot(ppp, add=T, main = "Riley County Quadrat Well Count", cex = 1, pch = "+", cols = "black", legend = FALSE, use.marks = FALSE)




######Do the same thing but tesselate and bin precipitation and temperature data (rasterstack?)
climate_data <- prism_archive_ls() %>%  
  pd_stack(.)  
clim <-projectRaster(climate_data, crs="+proj=tmerc +lat_0=31 +lon_0=-106.25 +k=0.9999 +x_0=500000 +y_0=0 +datum=NAD83 +units=m +no_defs")


crop_clim <- crop(clim, aoi)
mask_clim <- mask(crop_clim, mask=crop_clim)

par(mar=c(1,1,1,1))
plot(crop_clim, main = "NM climate map")

quantile(mask_clim)

#define climate zones/bins
temp_zones <- reclassify(crop_clim$PRISM_tmean_30yr_normal_4kmM4_annual_bil,
                         c(0, 0.784, 1,
                           0.784, 10.18, 2,
                           10.18, 12.99, 3,
                           12.99, 15.48, 4,
                           15.48, 18.91, 5))

temp_zones <- ratify(temp_zones)
plot(temp_zones, main = "Temperature Zones")
temp_zones <- as.im.RasterLayer(temp_zones) #Convert to a Spatstat-compatible object



#Also zone the precipitation data from PRISM normals
ppt_zones <- reclassify(crop_clim$PRISM_ppt_30yr_normal_4kmM4_annual_bil,
                      c(0, 153.41, 1,
                        153.1, 283.7, 2,
                        283.7, 345.5, 3,
                        345.5, 409, 4,
                        409, 1058.61, 5))
ppt_zones <- ratify(ppt_zones)
plot(ppt_zones, main = "Precipitation Zones")
ppt_zones <- as.im.RasterLayer(ppt_zones) #Convert to a Spatstat-compatible object


#Tesselate the images of climate data
tes_ppt <- tess(image = ppt_zones)
tes_temp <- tess(image = temp_zones)


plot(tes_ppt, main = "Tesselated Precipitation Zones")
plot(tes_temp, main = "Tesselated Temperature Zones")



####
####
####Model fitting with the tesselations for elevation and climate
#Try model fitting with elevation
fit1 = ppm(ppp, ~tes_elev)
fit2 = ppm(ppp, ~elev)
summary(fit1)
summary(fit2)

#model fitting with precipitation
fit3 = ppm(ppp, ~tes_ppt)
fit4 = ppm(ppp, ~ppt_zones)

#model fitting with temperature
fit5 = ppm(ppp, ~tes_temp)
fit6 = ppm(ppp, ~temp_zones)





####EXTRACT tesselated climate and elevation data and add to main dataframe (env.v)
# Extract the  values from the raster stack for those sites 
tessElev<-st_extract(ppt_zones, st_geometry(dat.sf))



#######

#######

#######
###TRY Chatgpt for writing code to characterize ppp across different landscape grouping/categories.
# Load required packages
library(spatstat)
dat<-st_sf(geometry = dat.sf)

# Convert sf object to ppp object
dat2 <- as.ppp(st_as_sfc(dat))

# Define the function
characterize_patterns <- function(data, category_col, pattern_col) {
  
  # Extract categories
  categories <- levels(data[[category_col]])
  
  # Initialize empty list to store results
  results <- list()
  
  # Loop through categories
  for (i in seq_along(categories)) {
    
    # Subset the data by category
    subset <- data[data[[category_col]] == categories[i], ]
    
    # Create a point pattern from the subset
    pp <- as.ppp(subset[[pattern_col]])
    
    # Calculate summary statistics for the point pattern
    summary_stats <- list(category = categories[i],
                          intensity = density(pp),
                          csr_pvalue = csr.test(pp)$p.value,
                          jitter_pvalue = jitter.test(pp)$p.value)
    
    # Add summary statistics to results list
    results[[i]] <- summary_stats
  }
  
  # Combine results into a data frame
  results_df <- do.call(rbind, results)
  
  # Return results data frame
  return(results_df)
}

# Example usage
#data(dat.sf)
dat.sf$ecoreg3 <- factor(dat.sf$ecoreg3)
results <- characterize_patterns(dat.sf, "ecoreg3", "geometry")
print(results)