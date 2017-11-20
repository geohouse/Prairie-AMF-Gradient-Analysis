# Permanova tests with differences in site disturbance history and side of the precipiation gradient
# to generate results in Table 1.

# Working correctly 111917

library("vegan")

# Import the OTU table
OTUTableGradient <- read.table("~/Prairie-AMF-Gradient-Analysis-master/GradientAMF_OTUTable_forR.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Import the AM fungal OTU numbers (identified from the phylogeny with reference sequences, Appendix 1, Figure S1).
AMFOTUNumbers <- scan("~/Prairie-AMF-Gradient-Analysis-master/GradientAMF_AMF_OTUNumbers_forR.txt")

# Remove the labels and other information columns to keep only the OTU table counts
OTUframeGradient <- data.frame(OTUTableGradient[,5:length(OTUTableGradient[1,])])

# This works correctly to subset only the AMF OTU columns by moving them by name (eg X1,X12) 
gradientSamplesOnlyGradientFrame <- OTUframeGradient[,paste("X",AMFOTUNumbers,sep="")]

# Add the sample names back to the otu table data frame as rownames
gradientSamplesOnlyGradientFrame_names <- OTUTableGradient[,4]
rownames(gradientSamplesOnlyGradientFrame) <- gradientSamplesOnlyGradientFrame_names

# Import environmental and soil chemical values for each seq barcode
gradientEnvAndSoilChemValues <- read.table("~/Prairie-AMF-Gradient-Analysis-master/GradientAMF_environ_soilNutrient_values_forR.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

gradientEnvAndSoilChemValues <- gradientEnvAndSoilChemValues[,-1:-3]
rownames(gradientEnvAndSoilChemValues) <- gradientSamplesOnlyGradientFrame_names

# East/West geographic groups.
# West: Klemme, Stillwater, Hays, Ft. Riley (2013), Konza (2015)
# East: Lawrence (Rockefeller and Welda), SW MO, Morris, IL (2013), IL (2015)
geogrGroupWE <- c(rep("W", 7), rep("E", 2), rep("W", 4), rep("E", 8), rep("W", 17), rep("E", 13),
                  rep("W", 4), rep("E", 17), rep("W", 20), rep("E", 3), rep("W", 6), rep("E", 7))

gradientEnvAndSoilChemValues$locEW <- geogrGroupWE
# --------------------------

numSeqs <- rowSums(gradientSamplesOnlyGradientFrame)

gradientEnvAndSoilChemValues$numSeqs <- numSeqs

# Subset for just remnant or disturbed sites
# For Remnant sites
gradientEnvAndSoilChemValues_remn <- gradientEnvAndSoilChemValues[gradientEnvAndSoilChemValues$Disturbed_Remnant == 'Remnant',]
gradientSamplesOnlyGradientFrame_remn <- gradientSamplesOnlyGradientFrame[gradientEnvAndSoilChemValues$Disturbed_Remnant == 'Remnant',]

# For Disturbed sites
gradientEnvAndSoilChemValues_dist <- gradientEnvAndSoilChemValues[gradientEnvAndSoilChemValues$Disturbed_Remnant == 'Disturbed',]
gradientSamplesOnlyGradientFrame_dist <- gradientSamplesOnlyGradientFrame[gradientEnvAndSoilChemValues$Disturbed_Remnant == 'Disturbed',]


# ----------
# This is for Table 1 all samples. 072517.
set.seed(618)
siteHist_then_eastWest_interactionTest_logSeqCovariate <- vegan::adonis2(formula = gradientSamplesOnlyGradientFrame ~ Disturbed_Remnant + locEW + Disturbed_Remnant * locEW + log10(numSeqs), 
                                                      data = gradientEnvAndSoilChemValues, permutations = 5000, method = "morisita", binary = FALSE)

print(siteHist_then_eastWest_interactionTest_logSeqCovariate)
# --------------

# --------------
# This is for Table 1 only remnant samples 072517

set.seed(522)
geogrLocControlMorisitaPermanova_remn <- vegan::adonis2(formula = gradientSamplesOnlyGradientFrame_remn ~  locEW + log10(numSeqs),
                                                 data = gradientEnvAndSoilChemValues_remn, permutations = 5000, by = 'margin', method = "morisita", binary = FALSE)
print(geogrLocControlMorisitaPermanova_remn)

# ------------

# ------------
# This is for Table 1 only disturbed samples 110717
set.seed(258)
geogrLocControlMorisitaPermanova_dist <- vegan::adonis2(formula = gradientSamplesOnlyGradientFrame_dist ~  locEW + log10(numSeqs) + Known_History_Of_Tilling_Or_Construction,
                                                 data = gradientEnvAndSoilChemValues_dist, permutations = 5000, by = 'margin', method = "morisita", binary = FALSE)
print(geogrLocControlMorisitaPermanova_dist)

# -------------
