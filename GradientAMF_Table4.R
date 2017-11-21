# This is a permanova test for the capscale constrained PCoA with environmental and soil factors. 
# This tests for significant differences in AM fungal communities among the 
# samples that have full soil test results (including C/N) in response to disturbance history, each soil var + precip, and all
# interactions between disturbance history and soil vars/precip.

# This is working correctly 112017.

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

distGroup <- c(rep("Disturbed",11),rep("Remnant",31),rep("Disturbed",3),rep("Remnant",1),rep("Disturbed",1),rep("Remnant",1),rep("Disturbed",3),
               rep("Remnant",1),rep("Disturbed",2),rep("Remnant",1),
               rep("Disturbed",4),rep("Remnant",1),rep("Disturbed",1),rep("Remnant",1),rep("Disturbed",1),
               rep("Remnant",1),rep("Disturbed",2),rep("Remnant",6))


# add the distGroup and the locationAbbrvs to the Env and Soil chem values
soilNutrValuesForCapscale$distGroup <- distGroup

numSeqs <- rowSums(gradientSamplesOnlyGradientFrame)

soilNutrValuesForCapscale$numSeqs <- numSeqs

# ========================
# Including or excluding Avg_longitude in the model makes no diff to the R2 for the main effects and only minimal diff for interaction for the basic model, Bray2 and Ca (spot tested)

permanovaOutputHolder <- numeric(length = 7)
names(permanovaOutputHolder) <- c("changing_var", "precip", "disturbHist", "avgLongitude", "log10NumSeqs", "precip_x_distubHist", "residuals")

set.seed(605)
soilVarPermanova_basicInteract <- adonis(formula = gradientSamplesOnlyGradientFrame ~ NOAA_AnnualPrecip_mm * distGroup + Avg_longitude + log10(numSeqs),
                                            data = soilNutrValuesForCapscale, permutations = 5000, method = "morisita", binary = FALSE)

holder <- cbind(c(NA, soilVarPermanova_basicInteract$aov.tab$R2[1:6]), c(NA, soilVarPermanova_basicInteract$aov.tab$`Pr(>F)`[1:6]))

colnames(holder) <- c("basicInteract_R2", "basicInteract_Pval")

permanovaOutputHolder <- holder
rownames(permanovaOutputHolder) <- c("changing_var", "precip", "disturbHist", "avgLongitude", "log10NumSeqs", "precip_x_distubHist", "residuals")

print("soilVarPermanova_basicInteract:")
print(soilVarPermanova_basicInteract)


set.seed(605)
soilVarPermanova_interact_w_pH <- adonis(formula = gradientSamplesOnlyGradientFrame ~ soil_pH + NOAA_AnnualPrecip_mm * distGroup + Avg_longitude + log10(numSeqs),
                                         data = soilNutrValuesForCapscale, permutations = 5000, method = "morisita", binary = FALSE)

holder <- cbind(soilVarPermanova_interact_w_pH$aov.tab$R2[1:7], soilVarPermanova_interact_w_pH$aov.tab$`Pr(>F)`[1:7])
colnames(holder) <- c("pH_R2", "pH_Pval")
permanovaOutputHolder <- cbind(permanovaOutputHolder, holder)

print("soilVarPermanova_interact_w_pH")
print(soilVarPermanova_interact_w_pH)


set.seed(605)
soilVarPermanova_interact_w_Bray1P <- adonis(formula = gradientSamplesOnlyGradientFrame ~ log10(Bray1P_ppm) + NOAA_AnnualPrecip_mm * distGroup + Avg_longitude + log10(numSeqs),
                                         data = soilNutrValuesForCapscale, permutations = 5000, method = "morisita", binary = FALSE)

holder <- cbind(soilVarPermanova_interact_w_Bray1P$aov.tab$R2[1:7], soilVarPermanova_interact_w_Bray1P$aov.tab$`Pr(>F)`[1:7])
colnames(holder) <- c("log10(Bray1P)_R2", "log10(Bray1P)_Pval")
permanovaOutputHolder <- cbind(permanovaOutputHolder, holder)

print("soilVarPermanova_interact_w_Bray1P")
print(soilVarPermanova_interact_w_Bray1P)


set.seed(605)
soilVarPermanova_interact_w_Bray2P <- adonis(formula = gradientSamplesOnlyGradientFrame ~ log10(Bray2P_ppm) + NOAA_AnnualPrecip_mm * distGroup + Avg_longitude + log10(numSeqs),
                                             data = soilNutrValuesForCapscale, permutations = 5000, method = "morisita", binary = FALSE)

holder <- cbind(soilVarPermanova_interact_w_Bray2P$aov.tab$R2[1:7], soilVarPermanova_interact_w_Bray2P$aov.tab$`Pr(>F)`[1:7])
colnames(holder) <- c("log10(Bray2P)_R2", "log10(Bray2P)_Pval")
permanovaOutputHolder <- cbind(permanovaOutputHolder, holder)

print("soilVarPermanova_interact_w_Bray2P")
print(soilVarPermanova_interact_w_Bray2P)


set.seed(605)
soilVarPermanova_interact_w_bicarbP <- adonis(formula = gradientSamplesOnlyGradientFrame ~ log10(bicarbP_ppm) + NOAA_AnnualPrecip_mm * distGroup + Avg_longitude + log10(numSeqs),
                                             data = soilNutrValuesForCapscale, permutations = 5000, method = "morisita", binary = FALSE)

holder <- cbind(soilVarPermanova_interact_w_bicarbP$aov.tab$R2[1:7], soilVarPermanova_interact_w_bicarbP$aov.tab$`Pr(>F)`[1:7])
colnames(holder) <- c("log10(bicarbP)_R2", "log10(bicarbP)_Pval")
permanovaOutputHolder <- cbind(permanovaOutputHolder, holder)

print("soilVarPermanova_interact_w_bicarbP")
print(soilVarPermanova_interact_w_bicarbP)


set.seed(605)
soilVarPermanova_interact_w_CEC <- adonis(formula = gradientSamplesOnlyGradientFrame ~ log10(CEC_meqPer100g) + NOAA_AnnualPrecip_mm * distGroup + Avg_longitude + log10(numSeqs),
                                              data = soilNutrValuesForCapscale, permutations = 5000, method = "morisita", binary = FALSE)

holder <- cbind(soilVarPermanova_interact_w_CEC$aov.tab$R2[1:7], soilVarPermanova_interact_w_CEC$aov.tab$`Pr(>F)`[1:7])
colnames(holder) <- c("log10(CEC)_R2", "log10(CEC)_Pval")
permanovaOutputHolder <- cbind(permanovaOutputHolder, holder)


print("soilVarPermanova_interact_w_CEC")
print(soilVarPermanova_interact_w_CEC)


set.seed(605)
soilVarPermanova_interact_w_CN <- adonis(formula = gradientSamplesOnlyGradientFrame ~ C_N + NOAA_AnnualPrecip_mm * distGroup + Avg_longitude + log10(numSeqs),
                                          data = soilNutrValuesForCapscale, permutations = 5000, method = "morisita", binary = FALSE)

holder <- cbind(soilVarPermanova_interact_w_CN$aov.tab$R2[1:7], soilVarPermanova_interact_w_CN$aov.tab$`Pr(>F)`[1:7])
colnames(holder) <- c("CN_R2", "CN_Pval")
permanovaOutputHolder <- cbind(permanovaOutputHolder, holder)

print("soilVarPermanova_interact_w_CN")
print(soilVarPermanova_interact_w_CN)


set.seed(605)
soilVarPermanova_interact_w_OM <- adonis(formula = gradientSamplesOnlyGradientFrame ~ perc_OM + NOAA_AnnualPrecip_mm * distGroup + Avg_longitude + log10(numSeqs),
                                         data = soilNutrValuesForCapscale, permutations = 5000, method = "morisita", binary = FALSE)

holder <- cbind(soilVarPermanova_interact_w_OM$aov.tab$R2[1:7], soilVarPermanova_interact_w_OM$aov.tab$`Pr(>F)`[1:7])
colnames(holder) <- c("OM_R2", "OM_Pval")
permanovaOutputHolder <- cbind(permanovaOutputHolder, holder)


print("soilVarPermanova_interact_w_OM")
print(soilVarPermanova_interact_w_OM)


set.seed(605)
soilVarPermanova_interact_w_Mg <- adonis(formula = gradientSamplesOnlyGradientFrame ~ log10(Mg_ppm) + NOAA_AnnualPrecip_mm * distGroup + Avg_longitude + log10(numSeqs),
                                          data = soilNutrValuesForCapscale, permutations = 5000, method = "morisita", binary = FALSE)

holder <- cbind(soilVarPermanova_interact_w_Mg$aov.tab$R2[1:7], soilVarPermanova_interact_w_Mg$aov.tab$`Pr(>F)`[1:7])
colnames(holder) <- c("log10(Mg)_R2", "log10(Mg)_Pval")
permanovaOutputHolder <- cbind(permanovaOutputHolder, holder)


print("soilVarPermanova_interact_w_Mg")
print(soilVarPermanova_interact_w_Mg)


set.seed(605)
soilVarPermanova_interact_w_Ca <- adonis(formula = gradientSamplesOnlyGradientFrame ~ log10(Ca_ppm) + NOAA_AnnualPrecip_mm * distGroup + Avg_longitude + log10(numSeqs),
                                         data = soilNutrValuesForCapscale, permutations = 5000, method = "morisita", binary = FALSE)

holder <- cbind(soilVarPermanova_interact_w_Ca$aov.tab$R2[1:7], soilVarPermanova_interact_w_Ca$aov.tab$`Pr(>F)`[1:7])
colnames(holder) <- c("log10(Ca)_R2", "log10(Ca)_Pval")
permanovaOutputHolder <- cbind(permanovaOutputHolder, holder)


print("soilVarPermanova_interact_w_Ca")
print(soilVarPermanova_interact_w_Ca)


set.seed(605)
soilVarPermanova_interact_w_K <- adonis(formula = gradientSamplesOnlyGradientFrame ~ log10(K_ppm) + NOAA_AnnualPrecip_mm * distGroup + Avg_longitude + log10(numSeqs),
                                         data = soilNutrValuesForCapscale, permutations = 5000, method = "morisita", binary = FALSE)

holder <- cbind(soilVarPermanova_interact_w_K$aov.tab$R2[1:7], soilVarPermanova_interact_w_K$aov.tab$`Pr(>F)`[1:7])
colnames(holder) <- c("log10(K)_R2", "log10(K)_Pval")
permanovaOutputHolder <- cbind(permanovaOutputHolder, holder)


print("soilVarPermanova_interact_w_K")
print(soilVarPermanova_interact_w_K)


# Now determine which soil factor had the greatest 'pull' of assigned variation away from siteDistHist, aridity, and their interaction compared to the 
# basic model with no soil factors. These are the values in Table 4.

# Extract just the R2 values for the different variables for the diff permanova runs that included a soil var
# The 'changing_var' is the contribution of each soil variable (listed in the column names)
rsquareHolder <- permanovaOutputHolder[,c(3,5,7,9,11,13,15,17,19,21)]

# Now calc the diff between the proportion var explained by aridity and site history in the basic model (with no soil factors) with the proportion of var explained with 
# each of the soil vars.

# This is precip
precipVarExplainedDiff <- ((rsquareHolder[2,] - permanovaOutputHolder[2,1])/permanovaOutputHolder[2,1]) * 100
# This is site history
siteHistoryVarExplainedDiff <- ((rsquareHolder[3,] - permanovaOutputHolder[3,1])/permanovaOutputHolder[3,1]) * 100
# This is the interaction
precipSiteHistInteractVarExplainedDiff <- ((rsquareHolder[6,] - permanovaOutputHolder[6,1])/permanovaOutputHolder[6,1]) * 100

# The min value in each of these variance explained diff holders is the soil var that 'pulled' the most amount of var otherwise explained by either aridity or site hist,
# and therefore is most correlated with them:

soilVarMostPulledPrecipVar <- precipVarExplainedDiff[which(precipVarExplainedDiff == min(precipVarExplainedDiff))]
soilVarMostPulledSiteHistoryVar <- siteHistoryVarExplainedDiff[which(siteHistoryVarExplainedDiff == min(siteHistoryVarExplainedDiff))]
