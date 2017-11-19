# Unconstrained PCoA for visualizing AM fungal community differences with site disturbance for all samples,
# and differences with side of the precip. gradient (West/East) for samples from only remnant or disturbed sites.
# This makes the plots in Figure 2 (all three plots are made automatically; in RStudio can scroll through them in the Plots window).

# Working correctly 080117

library("vegan")

# Import the OTU table
OTUTableGradient <- read.table("~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/GradientCoordinated_OTUTable_forR.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Import the AM fungal OTU numbers (identified from the phylogeny with reference sequences, Appendix 1, Figure S1).
AMFOTUNumbers <- scan("~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/GradientCoordinated_AMF_OTUNumbers_forR.txt")

# Remove the labels and other information columns to keep only the OTU table counts
OTUframeGradient <- data.frame(OTUTableGradient[,5:length(OTUTableGradient[1,])])

# This works correctly to subset only the AMF OTU columns by moving them by name (eg X1,X12) 
gradientSamplesOnlyGradientFrame <- OTUframeGradient[,paste("X",AMFOTUNumbers,sep="")]

# Add the sample names back to the otu table data frame as rownames
gradientSamplesOnlyGradientFrame_names <- OTUTableGradient[,4]
rownames(gradientSamplesOnlyGradientFrame) <- gradientSamplesOnlyGradientFrame_names

# Import environmental and soil chemical values for each seq barcode; This is only used here for the remnant/disturbed
# site labels for each sample.
gradientEnvAndSoilChemValues <- read.table("~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/GradientCoordinated_environ_soilNutrient_values_forR.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

gradientEnvAndSoilChemValues <- gradientEnvAndSoilChemValues[,-1:-3]
rownames(gradientEnvAndSoilChemValues) <- gradientSamplesOnlyGradientFrame_names

# East/West geographic groups.
# West: Klemme, Stillwater, Hays, Ft. Riley (2013), Konza (2015)
# East: Lawrence (Rockefeller and Welda), SW MO, Morris, IL (2013), IL (2015)
geogrGroupWE <- c(rep("W", 7), rep("E", 2), rep("W", 4), rep("E", 8), rep("W", 17), rep("E", 13),
                  rep("W", 4), rep("E", 17), rep("W", 20), rep("E", 3), rep("W", 6), rep("E", 7))

# make a new synthetic grouping factor by merging the sampling site and dist hist factors
EWSiteHistFactor <- factor(x = paste(geogrGroupWE, as.factor(gradientEnvAndSoilChemValues$Disturbed_Remnant)), levels = c("W Disturbed", "E Disturbed", "W Remnant", "E Remnant"))

# Run the PCoA using vegan::capscale using Morisita's dissimilarity index from 
# the raw (untransformed) counts in the OTU table.
# capscale requires a formula; this is how to get it to run an unconstrained PCoA (following the manual)
morisita.PCoA <- vegan::capscale(gradientSamplesOnlyGradientFrame ~ 1, dist = "morisita", binary = FALSE)

# This is hardcoded for the model as of 020817; check the prop. of var using summary(morisita.PCoA)
xlabel <- "PCoA 1 (15% of variance explained)"
ylabel <- "PCoA 2 (13% of variance explained)"

# All samples with points in the plot colored by disturbance history of the site (orange = remnant, blue = disturbed)
# Fig 2A
# ==================================

assignPointColors <- function(distLabel){
  colorForPoint <- ifelse(distLabel == "Remnant", "#FF7F00", "#377EB8")
  return(colorForPoint)
}

# Need to assign a color for each point (sample) based on its disturbance history.
colorHolder <- sapply(gradientEnvAndSoilChemValues$Disturbed_Remnant, FUN = assignPointColors)

# Start by making an empty plot.
initialPlot <- plot(morisita.PCoA, type = "n", scaling = 'sites', xlab = xlabel, ylab = ylabel, ylim = c(-1, 1), xlim = c(-1, 1))
# Then put the points on the plot made in the previous step, with the correct colors based on the disturbance history
# for the site.
points(morisita.PCoA, col = colorHolder, scaling = 'sites', pch = 19, cex = 2)
#===============================


# Only samples from remnant sites with points in the plot colored by the location of the site
# on the precipitation gradient (purple = West, green = East)
# Fig 2B
# ==================================

# Associate each point from the overall PCoA a color - gray for points representing samples from
# disturbed sites, purple = West remnants, green = East remnants
assignPointColors <- function(siteLabel){
  # Get the index of the site for the current point out of the siteGroupFactor
  factorLevelForPoint <- which(levels(EWSiteHistFactor) == siteLabel)
  print(siteLabel)
  if(siteLabel == "W Disturbed" | siteLabel == "E Disturbed"){
    colorForPoint <- "gray90"
  } else{
    colorForPoint <- ifelse(siteLabel == "W Remnant", "#984EA3", "#4DAF4A")
  }
  return(colorForPoint)
}

# Need to assign a color for each point based on its siteGroupFactor.
colorHolder_remn <- sapply(EWSiteHistFactor, FUN = assignPointColors)

# Start by making an empty plot.
initialPlot <- plot(morisita.PCoA, type = "n", scaling = 'sites', xlab = xlabel, ylab = ylabel, ylim = c(-1, 1), xlim = c(-1, 1))
# Then put the points on the plot made in the previous step, with the correct colors based on the site's location on 
# the precipitation gradient (West/East)
points(morisita.PCoA, col = colorHolder_remn, scaling = 'sites', pch = 19, cex = 2)
#===============================

# Only samples from disturbed sites with points in the plot colored by the location of the site
# on the precipitation gradient (purple = West, green = East)
# Fig 2C
# ==================================

# Associate each point from the overall PCoA a color - gray for points representing samples from
# remnant sites, purple = West disturbed sites, green = East disturbed sites
assignPointColors <- function(siteLabel){
    # Get the index of the site for the current point out of the siteGroupFactor
    factorLevelForPoint <- which(levels(EWSiteHistFactor) == siteLabel)
    print(siteLabel)
    if(siteLabel == "W Remnant" | siteLabel == "E Remnant"){
        colorForPoint <- "gray90"
    } else{
        colorForPoint <- ifelse(siteLabel == "W Disturbed", "#984EA3", "#4DAF4A")
    }
    return(colorForPoint)
}

# Need to assign a color for each point based on its siteGroupFactor.
colorHolder_dist <- sapply(EWSiteHistFactor, FUN = assignPointColors)

# Start by making an empty plot.
initialPlot <- plot(morisita.PCoA, type = "n", scaling = 'sites', xlab = xlabel, ylab = ylabel, ylim = c(-1, 1), xlim = c(-1, 1))
# Then put the points on the plot made in the previous step, with the correct colors based on the site's location on 
# the precipitation gradient (West/East)
points(morisita.PCoA, col = colorHolder_dist, scaling = 'sites', pch = 19, cex = 2)
#===============================


