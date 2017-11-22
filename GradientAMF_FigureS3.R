# This is for canonical analysis of PCoA (CAP, or capscale in vegan) using only the samples with 
# full soil var panel measured (including C/N) to see how diff soil vars may predict AM fungal community 
# composition in a constrained ordination method that allows use of Morisita dissimilarity metric.
# Using partial capscale to partial out the geographic pseudo-replication within sites, and then
# estimating effects due to soil vars. This is Figure S3.

# This is working correctly 112117.

library("vegan")

# Import the OTU table
OTUTableGradient <- read.table("~/Prairie-AMF-Gradient-Analysis-master/GradientAMF_OTUTable_forR.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Import the AM fungal OTU numbers (identified from the phylogeny with reference sequences, Appendix 1, Figure S1).
AMFOTUNumbers <- scan("~/Prairie-AMF-Gradient-Analysis-master/GradientAMF_AMF_OTUNumbers_forR.txt")

# Remove the labels and other information columns to keep only the OTU table counts
OTUframeGradient <- data.frame(OTUTableGradient[,5:length(OTUTableGradient[1,])])

# This works correctly to subset only the AMF OTU columns by moving them by name (eg X1,X12) 
gradientSamplesOnlyGradientFrame <- OTUframeGradient[,paste("X",AMFOTUNumbers,sep="")]

# Import environmental and soil chemical values for each seq barcode
gradientEnvAndSoilChemValues <- read.table("~/Prairie-AMF-Gradient-Analysis-master/GradientAMF_environ_soilNutrient_values_forR.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Subset the environmental values to only the samples that contain C/N measures, and exclude those from
# Klemme (so no samples from Oklahoma are considered)
soilNutrValuesForCapscale <- gradientEnvAndSoilChemValues[which(!is.na(gradientEnvAndSoilChemValues$C_N)),]

soilNutrValuesForCapscale <- soilNutrValuesForCapscale[which(soilNutrValuesForCapscale$Site != "Klemme"),]

# Subset the OTU table keeping only the sample numbers that are still in the environmental data table
gradientSamplesOnlyGradientFrame <- gradientSamplesOnlyGradientFrame[which(OTUTableGradient$Coordinated_sampleNumber %in% soilNutrValuesForCapscale$CoordinatedSampleNumber),]

geogrDistGroup <- c(rep("KS_Dist.",7),rep("MO_Dist.",4),
                    rep("KS_Remn.",16),rep("KS_Remn.",5),rep("MO_Remn.",10),rep("IL_Dist.",3),
                    rep("KS_Remn.",1),rep("KS_Dist.",1),rep("KS_Remn.",1),
                    rep("KS_Dist.",3),rep("KS_Remn.",1),rep("KS_Dist.",2),rep("KS_Remn.",1),
                    rep("KS_Dist.",1),rep("IL_Dist.",3),rep("KS_Remn.",1),rep("KS_Dist.",1),rep("KS_Remn.",1),
                    rep("KS_Dist.",1),rep("KS_Remn.",1),rep("KS_Dist.",1),rep("IL_Dist.",1),rep("IL_Remn.",6))

#------------------------------


# Constrained partial PCoA with constraint based on site location (combination of latitude and longitude)
constrPCoAResults <- capscale(gradientSamplesOnlyGradientFrame ~ perc_OM  + log(Bray1P_ppm) + log(Bray2P_ppm) + log(bicarbP_ppm) + soil_pH + log(CEC_meqPer100g) + C_N + log(K_ppm) + 
                                log(Ca_ppm) + log(Mg_ppm) + Condition(Avg_longitude + Avg_latitude),
                              data = soilNutrValuesForCapscale, dist = 'morisita')


# plot(constrPCoAResults, scaling = 'species')

# -----------------
## This is the plotting.

# Need to change this to get the plots for the different soil variables.
# For OM, use:
# environVarToFit <- soilNutrValuesForCapscale$perc_OM
# environVarToFit <- soilNutrValuesForCapscale$soil_pH
# environVarToFit <- log10(soilNutrValuesForCapscale$Bray2P_ppm)
# environVarToFit <- log10(soilNutrValuesForCapscale$Mg_ppm)

environVarToFit <- soilNutrValuesForCapscale$perc_OM
plotTitle <- "Percent Organic Matter"

# Can use choices=c(1,3) to plot PCoA1 and PCoA3 instead (or any other pair). C/N is a better fit on PCoA 1 and 3 (it's terrible
# with PCoA 1 and 2).
morisita.ordisurf <- ordisurf(constrPCoAResults ~ environVarToFit, choices = c(1,2), scaling = "sites", cex = 0, labcex = 1.5, lwd.cl = 2, col = "black", xlab = "PCoA 1", ylab = "PCoA 2", main = plotTitle)

distPointColor <- "#377eb8"
# Orange
remnPointColor <- "#ff7f00"

# Calculates and plots "spiders" connecting points from each geogr/dist.

meanKSDistLoc <- c(mean(scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "KS_Dist.",1]), mean(scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "KS_Dist.",2]))
numKSDistLoc <- length(geogrDistGroup[geogrDistGroup == "KS_Dist."])
meanKSDistLocVec <- matrix(rep(meanKSDistLoc), numKSDistLoc, ncol = 2, byrow = TRUE)

meanKSRemnLoc <- c(mean(scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "KS_Remn.",1]), mean(scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "KS_Remn.",2]))
numKSRemnLoc <- length(geogrDistGroup[geogrDistGroup == "KS_Remn."])
meanKSRemnLocVec <- matrix(rep(meanKSRemnLoc), numKSRemnLoc, ncol = 2, byrow = TRUE)

meanMODistLoc <- c(mean(scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "MO_Dist.",1]), mean(scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "MO_Dist.",2]))
numMODistLoc <- length(geogrDistGroup[geogrDistGroup == "MO_Dist."])
meanMODistLocVec <- matrix(rep(meanMODistLoc), numMODistLoc, ncol = 2, byrow = TRUE)

meanMORemnLoc <- c(mean(scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "MO_Remn.",1]), mean(scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "MO_Remn.",2]))
numMORemnLoc <- length(geogrDistGroup[geogrDistGroup == "MO_Remn."])
meanMORemnLocVec <- matrix(rep(meanMORemnLoc), numMORemnLoc, ncol = 2, byrow = TRUE)


meanILDistLoc <- c(mean(scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "IL_Dist.",1]), mean(scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "IL_Dist.",2]))
numILDistLoc <- length(geogrDistGroup[geogrDistGroup == "IL_Dist."])
meanILDistLocVec <- matrix(rep(meanILDistLoc), numILDistLoc, ncol = 2, byrow = TRUE)

meanILRemnLoc <- c(mean(scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "IL_Remn.",1]), mean(scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "IL_Remn.",2]))
numILRemnLoc <- length(geogrDistGroup[geogrDistGroup == "IL_Remn."])
meanILRemnLocVec <- matrix(rep(meanILRemnLoc), numILRemnLoc, ncol = 2, byrow = TRUE)

segments(meanKSDistLocVec[,1], meanKSDistLocVec[,2], scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "KS_Dist.",1], scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "KS_Dist.",2], col = rgb(0.6,0.6,0.6))
segments(meanKSRemnLocVec[,1], meanKSRemnLocVec[,2], scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "KS_Remn.",1], scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "KS_Remn.",2], col = rgb(0.6,0.6,0.6))

segments(meanMODistLocVec[,1], meanMODistLocVec[,2], scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "MO_Dist.",1], scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "MO_Dist.",2], col = rgb(0.6,0.6,0.6))
segments(meanMORemnLocVec[,1], meanMORemnLocVec[,2], scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "MO_Remn.",1], scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "MO_Remn.",2], col = rgb(0.6,0.6,0.6))

segments(meanILDistLocVec[,1], meanILDistLocVec[,2], scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "IL_Dist.",1], scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "IL_Dist.",2], col = rgb(0.6,0.6,0.6))
segments(meanILRemnLocVec[,1], meanILRemnLocVec[,2], scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "IL_Remn.",1], scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "IL_Remn.",2], col = rgb(0.6,0.6,0.6))

points(scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "KS_Dist.",1],scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "KS_Dist.",2], col = "black", bg = distPointColor, pch = 24, cex = 1)
points(scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "KS_Remn.",1],scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "KS_Remn.",2], col = "black", bg = remnPointColor,pch = 24, cex = 1)
points(scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "MO_Dist.",1],scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "MO_Dist.",2], col = "black", bg = distPointColor,pch = 22, cex = 1)
points(scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "MO_Remn.",1],scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "MO_Remn.",2], col = "black", bg = remnPointColor,pch = 22, cex = 1)
points(scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "IL_Dist.",1],scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "IL_Dist.",2], col = "black", bg = distPointColor, pch = 23, cex = 1)
points(scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "IL_Remn.",1],scores(constrPCoAResults, scaling = "sites")$sites[geogrDistGroup == "IL_Remn.",2], col = "black", bg = remnPointColor,pch = 23, cex = 1)

# Make the points for the mean location text
points(meanKSDistLoc[1],meanKSDistLoc[2], col = "black", bg = distPointColor, pch = 24, cex = 5)
points(meanKSRemnLoc[1],meanKSRemnLoc[2], col = "black", bg = remnPointColor, pch = 24, cex = 5)
points(meanMODistLoc[1],meanMODistLoc[2], col = "black", bg = distPointColor, pch = 22, cex = 5)
points(meanMORemnLoc[1],meanMORemnLoc[2], col = "black", bg = remnPointColor, pch = 22, cex = 5)
points(meanILDistLoc[1],meanILDistLoc[2], col = "black", bg = distPointColor, pch = 23, cex = 5)
points(meanILRemnLoc[1],meanILRemnLoc[2], col = "black", bg = remnPointColor, pch = 23, cex = 5)

# The text labels for the mean location
text(meanKSDistLoc[1], meanKSDistLoc[2], labels = "KS", col = "white")
text(meanKSRemnLoc[1], meanKSRemnLoc[2], labels = "KS")

text(meanMODistLoc[1], meanMODistLoc[2], labels = "MO", col = "white")
text(meanMORemnLoc[1], meanMORemnLoc[2], labels = "MO")

text(meanILDistLoc[1], meanILDistLoc[2], labels = "IL", col = "white")
text(meanILRemnLoc[1], meanILRemnLoc[2], labels = "IL")

print(summary(morisita.ordisurf))
