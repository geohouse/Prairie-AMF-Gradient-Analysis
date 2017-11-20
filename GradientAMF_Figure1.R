# This creates the map of sampled sites with points colored by site precipitation
# and with different shapes for remnant or disturbed sites.

# This is working correctly 111917

library(maps)
library(mapdata)
library(raster)
library(ggplot2)
library(RColorBrewer)

# import sampling sites
gradientSites <- read.table("~/Prairie-AMF-Gradient-Analysis-master/GradientAMF_siteLocations_forFigure1.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

uniquePrecipValues <- unique(gradientSites$NOAA_Precip[order(gradientSites$NOAA_Precip)])

# length 11 from Color Brewer
colors <- c("#8E0152","#C51B7D","#DE77AE","#F1B6DA","#FDE0EF","#F7F7F7","#E6F5D0","#B8E186","#7FBC41","#4D9221","#276419")

uniquePrecipValues <- cbind(uniquePrecipValues, colors)

matchPrecipToColors <- function(row){
    precip <- as.numeric(row[6])
    colorMatch <- which(uniquePrecipValues[,1] == precip)
    color <- uniquePrecipValues[colorMatch,2]
    return(color)
}

precipColors <- apply(gradientSites,MARGIN = 1,FUN = matchPrecipToColors)

gradientSites$pointShape <- ifelse(gradientSites$Disturbed_Remnant %in% c("Disturbed", "Disturbed_Brome", "Disturbed_OWB"),24,21)

stateMapData <- map_data('state')

selectedStates <- stateMapData$region %in%  c("oklahoma", "kansas","missouri",
                                              "illinois","iowa", "minnesota",
                                              "indiana", "nebraska")

selectedStateMapData <- stateMapData[selectedStates,]

remnantData <- gradientSites[which(gradientSites$Disturbed_Remnant == "Remnant"),]

disturbedData <- gradientSites[which(gradientSites$Disturbed_Remnant %in% c("Disturbed", "Disturbed_Brome", "Disturbed_OWB")),]

gradientMap <- ggplot() +
    geom_polygon(data = stateMapData, aes(x = long, y = lat, group = group), fill = "gray90", linetype = 1, size = 0.2, color = "gray70") + 
    
    geom_point(data  = disturbedData, aes(x = Avg_longitude, y = Avg_latitude, fill = NOAA_Precip),
             pch = disturbedData$pointShape,
             cex = 6, stroke = 0.4) +

    geom_point(data  = remnantData, aes(x = Avg_longitude, y = Avg_latitude, fill = NOAA_Precip), 
             pch = remnantData$pointShape,
             cex = 6, stroke = 0.4) + 
    scale_fill_gradientn(name = "Avg. annual \n precip. (mm)", colours=brewer.pal(9,"BrBG"), limits = c(500,1200), breaks = seq(500,1200,100)) + 
    coord_map(xlim = c(-101,-87), ylim = c(34,42.5)) + 
    scale_y_continuous(name="Latitude") + 
    scale_x_continuous(name="Longitude") + 
    annotate("text", x = c(-97,-97,-92,-89), y = c(35, 38, 38, 39), label = c("OK", "KS", "MO", "IL"), size = 7) + 
    geom_vline(xintercept = -96, color = "gray30", linetype = "dashed", size = 0.2) + 
  theme(legend.key.height = unit(2,"cm"), legend.key.width = unit(1,"cm"),
        legend.text = element_text(size = 18), legend.title = element_text(size = 18), 
        legend.position = "right", axis.title.x=element_text(size = 18), 
        axis.text.x = element_text(size = 18), axis.title.y=element_text(size = 18),
        axis.text.y = element_text(size = 18))

print(gradientMap)

