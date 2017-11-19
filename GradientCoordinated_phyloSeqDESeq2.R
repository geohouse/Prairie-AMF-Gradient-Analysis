# Finding OTUs with differential abundance either with site history (for all samples)
# or with side of the precip. gradient (West/East) for samples just from either remnant or disturbed
# sites. This uses DESEQ2 as an alternative to rarefaction using all OTU counts and models
# the variance from diff sample sizes.

# This gives the plot output shown in Figure 3 (after rotating, reflecting the bars for the West/East 
# comparisons across the axis, changing the order of genera (esp. Rhizoglomus) and other cosmetic modifications).

# These DESEQ2 results are required for the scripts that re-create Figures 4 and 5.

# The DESEQ2 and phyloseq packages are available on Bioconductor.

# Working correctly 080117

library('phyloseq')
library('DESeq2')

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

# Import environmental and soil chemical values for each seq barcode
gradientEnvAndSoilChemValues <- read.table("~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/GradientCoordinated_environ_soilNutrient_values_forR.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

subsetGradientEnvAndSoilChemValues_abbrv <- gradientEnvAndSoilChemValues[,-1:-3]
rownames(subsetGradientEnvAndSoilChemValues_abbrv) <- gradientSamplesOnlyGradientFrame_names

# Import the phylogeny-based genus attributions to each of the OTUs (from Figure S1). These are listed in 
# phylogenetic order e.g. phylogenetically closely related genera are listed immediately before/after one another)
genusAttributions <- read.table("~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/GradientCoordinated_OTU_GenusAttributions_forR.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

OTUNamesMatchOTUTable <- paste("X", genusAttributions$OTU_num, sep = "")

# Match is a great function to find the indices of the OTUs in the gradientFrame that match the 
# order of the phylogenetically-sorted OTU names in the OTUNamesMatch. This will be used for sorting the 
# DESeq2 output in phylogenetic order to see if there are consistent differences across the aridity gradient.
OTUTableIndicesInPhylogOrderForLabels <- match(OTUNamesMatchOTUTable, colnames(gradientSamplesOnlyGradientFrame))

distGroup <- subsetGradientEnvAndSoilChemValues_abbrv$Disturbed_Remnant

# For the revised analysis 120416
# East/West geographic groups.
# West: Klemme, Stillwater, Hays, Ft. Riley (2013), Konza (2015)
# East: Lawrence (Rockefeller and Welda), SW MO, Morris, IL (2013), IL (2015)
geogrGroupWE <- c(rep("W", 7), rep("E", 2), rep("W", 4), rep("E", 8), rep("W", 17), rep("E", 13),
                  rep("W", 4), rep("E", 17), rep("W", 20), rep("E", 3), rep("W", 6), rep("E", 7))


# Control for geographic location (need to lump remn/dist site at each geogr loc to 
# get a matrix model that is full rank)
siteGroup <- c(rep("Klemme", 3), rep("Hays", 4), rep("Lawrence", 2), rep("Konza", 2),
               rep("Klemme",1), rep("Konza", 1), rep("SW_MO", 4), rep("Morris",2),
               rep("Lawrence", 2), rep("Klemme", 8), rep("Hays", 5), rep("Konza", 4),
               rep("Lawrence", 7), rep("IL_2015", 6), rep("Konza", 4), rep("Lawrence", 1),
               rep("SW_MO", 8), rep("Morris",2), rep("IL_2015", 6), rep("Stillwater", 7),
               rep("FtRiley", 13), rep("IL_2013", 3), rep("FtRiley", 6), rep("IL_2013", 7))


subsetGradientEnvAndSoilChemValues_abbrv$siteGroup <- as.factor(siteGroup)

subsetGradientEnvAndSoilChemValues_abbrv$distGroup <- as.factor(distGroup)

subsetGradientEnvAndSoilChemValues_abbrv$geogrGroupWE <- as.factor(geogrGroupWE)

# Add pseudo count - this is REQUIRED  to be able to run the DESeq analysis.
gradientSamplesOnlyGradientFrame <- gradientSamplesOnlyGradientFrame + 1


#-----------------
# This is testing for disturbed/remn differences across both W and E sites, controlling for
# geographic location the best I can with DESeq2

pathToSaveDESeq2Output <- "~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/autoGeneratedScriptOutput/"
outputFileName <- "GradientCoordinated_DESEQ2_results_allSamples.r"
outputFile <- paste(pathToSaveDESeq2Output, "/", outputFileName, sep = "")

# Convert OTU table to phyloseq OTU table
phyloseq_OTUTable <- otu_table(gradientSamplesOnlyGradientFrame, taxa_are_rows = FALSE)

# This seems to work, but View() fails on it....
phyloseq_sampleData <- sample_data(subsetGradientEnvAndSoilChemValues_abbrv)

phyloseq_combined <- phyloseq(phyloseq_OTUTable, phyloseq_sampleData)

# Design '~ 1' is valid for not specifying a design in phyloseq_to_deseq2

# siteGroup and distGroup are factors in phyloseq_combined (from the data frame subsetGradientEnvAndSoilChemValues_abbrv used in the sample_data() call above)
# Listing siteGroup first in the model lets DESEQ2 take general geographic location into account, but by default it will
# test for differences in abundance for the last factor in the model (distGroup), like we want.
allSample_combinedDataForDeseq2 <- phyloseq_to_deseq2(phyloseq_combined, ~ siteGroup + distGroup)

# This tests the log2 fold change of the number of seqs in remnant compared to disturbed sites for each OTU (+ val = higher
# in remnant; - val = higher in disturbed sites), along with p-values for the log2 fold change.
allSample_DEseqOutput <- DESeq(allSample_combinedDataForDeseq2, fitType = 'local')

# access the raw results
DESeq_distHistCompare_rawResults <- results(allSample_DEseqOutput)

# sort results by p-value 
DESeq_distHistCompare_pValueSorted <- DESeq_distHistCompare_rawResults[order(DESeq_distHistCompare_rawResults$padj),]

# sort results by effect size. Commented out on purpose - not needed here.
# DESeq_distHistCompare_effectSizeSorted <- DESeq_distHistCompare_rawResults[order(abs(DESeq_distHistCompare_rawResults$log2FoldChange), decreasing = TRUE),]

# sort the results in phylogenetic order (using the order of entries of the genus attributions because those
# are listed in phylogenetic order e.g. phylogenetically closely related genera are listed immediately before/after one another)
DESeq_distHistCompare_phyloSorted <- DESeq_distHistCompare_rawResults[OTUTableIndicesInPhylogOrderForLabels,]

save(list = c("allSample_DEseqOutput","DESeq_distHistCompare_rawResults",
              "DESeq_distHistCompare_pValueSorted",
              "DESeq_distHistCompare_phyloSorted"), file = outputFile)

# Make the bar graph (Figure 3A)
# Ordered by phylogeny
# Leave width at 0.835 for the segments to fit well.

assignPointColors <- function(valueRow){
    print(valueRow[6])
    
    # adj. p value
    if(!is.na(valueRow[6]) & valueRow[6] < 0.05){
        # Log-2 fold change
        colorForBar <- ifelse(valueRow[2] > 0, "#FF7F00", "#377EB8")
    } else{
        colorForBar <- "white"    
    }
    return(colorForBar)
} 

# Need to assign a color for each point (sample) based on its disturbance history.
colorHolder <- apply(DESeq_distHistCompare_phyloSorted, MARGIN = 1, FUN = assignPointColors)


# phyloAdjPValColor <- ifelse(DESeq_distHistCompare_pValueSorted$padj < 0.05, 'green', 'white')
barplot(DESeq_distHistCompare_phyloSorted$log2FoldChange, col = colorHolder,
        ylab = 'Log 2 fold diff Remn (+) /Dist (-)', main = 'All Sites\n ordered by phylogeny',
        ylim = c(-6,6), width = 0.835)

#-----------------

#-----------------
# This is testing for East/West differences FOR ONLY Remnant sites; Based on how DESeq2 requires balanced cross designs
# for when have 2 diff factor vectors, I can't control for site location in a reasonable way within East/West; site location
# is only controlled as a group within East sites and as a group within West sites.

pathToSaveDESeq2Output <- "~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/autoGeneratedScriptOutput/"
outputFileName <- "GradientCoordinated_DESEQ2_results_remnantSites.r"
outputFile <- paste(pathToSaveDESeq2Output, "/", outputFileName, sep = "")

gradientSamplesOnlyGradientFrame_Remnant <- gradientSamplesOnlyGradientFrame[distGroup == 'Remnant',]
subsetGradientEnvAndSoilChemValues_abbrv_Remnant <- subsetGradientEnvAndSoilChemValues_abbrv[distGroup == 'Remnant',]

# Convert OTU table to phyloseq OTU table
phyloseq_OTUTable <- otu_table(gradientSamplesOnlyGradientFrame_Remnant, taxa_are_rows = FALSE)

# This seems to work, but View() fails on it....
phyloseq_sampleData <- sample_data(subsetGradientEnvAndSoilChemValues_abbrv_Remnant)

phyloseq_combined <- phyloseq(phyloseq_OTUTable, phyloseq_sampleData)

# Design '~ 1' is valid for not specifying a design in phyloseq_to_deseq2

remnantSample_combinedDataForDeseq2 <- phyloseq_to_deseq2(phyloseq_combined, ~ geogrGroupWE)

# This tests the log2 fold change of the number of seqs in remnant compared to disturbed sites for each OTU (+ val = higher
# in west; - val = higher in east sites), along with p-values for the log2 fold change.
remnantSample_DEseqOutput <- DESeq(remnantSample_combinedDataForDeseq2, fitType = 'local')

# access the raw results
DESeq_remnant_WE_Compare_rawResults <- results(remnantSample_DEseqOutput)

# sort results by p-value 
DESeq_remnant_WE_Compare_pValueSorted <- DESeq_remnant_WE_Compare_rawResults[order(DESeq_remnant_WE_Compare_rawResults$padj),]

# sort results by effect size. Commented out on purpose - not needed here.
# DESeq_remnant_WE_Compare_effectSizeSorted <- DESeq_remnant_WE_Compare_rawResults[order(abs(DESeq_remnant_WE_Compare_rawResults$log2FoldChange), decreasing = TRUE),]

# sort the results in phylogenetic order (using the order of entries of the genus attributions because those
# are listed in phylogenetic order e.g. phylogenetically closely related genera are listed immediately before/after one another)
DESeq_remnant_WE_Compare_phyloSorted <- DESeq_remnant_WE_Compare_rawResults[OTUTableIndicesInPhylogOrderForLabels,]

save(list = c("remnantSample_DEseqOutput","DESeq_remnant_WE_Compare_rawResults",
              "DESeq_remnant_WE_Compare_pValueSorted",
              "DESeq_remnant_WE_Compare_phyloSorted"), file = outputFile)

# Make the bar graph (Figure 3B)
# Ordered by phylogeny
# Leave width at 0.835 for the segments to fit well.

assignPointColors <- function(valueRow){
    print(valueRow[6])
    
    # adj. p value
    if(!is.na(valueRow[6]) & valueRow[6] < 0.05){
        # Log-2 fold change
        colorForBar <- ifelse(valueRow[2] > 0, "#984EA3", "#4DAF4A")
    } else{
        colorForBar <- "white"    
    }
    return(colorForBar)
} 

# Need to assign a color for each point (sample) based on its disturbance history.
colorHolder <- apply(DESeq_remnant_WE_Compare_phyloSorted, MARGIN = 1, FUN = assignPointColors)

# Ordered by phylogeny
# Leave width at 0.835 for the segments to fit well.
barplot(DESeq_remnant_WE_Compare_phyloSorted$log2FoldChange, col = colorHolder,
        ylab = 'Log 2 fold diff West (+) / East (-)', main = 'Remnant Sites Only\n ordered by phylogeny',
        ylim = c(-6,6), width = 0.835)



#-----------------


#-----------------
# This is testing for East/West differences FOR ONLY Disturbed sites; Based on how DESeq2 requires balanced cross designs
# for when have 2 diff factor vectors, I can't control for site location in a reasonable way within East/West; site location
# is only controlled as a group within East sites and as a group within West sites.

pathToSaveDESeq2Output <- "~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/autoGeneratedScriptOutput/"
outputFileName <- "GradientCoordinated_DESEQ2_results_disturbedSites.r"
outputFile <- paste(pathToSaveDESeq2Output, "/", outputFileName, sep = "")

gradientSamplesOnlyGradientFrame_Disturbed <- gradientSamplesOnlyGradientFrame[distGroup == 'Disturbed',]
subsetGradientEnvAndSoilChemValues_abbrv_Disturbed <- subsetGradientEnvAndSoilChemValues_abbrv[distGroup == 'Disturbed',]

# Convert OTU table to phyloseq OTU table
phyloseq_OTUTable <- otu_table(gradientSamplesOnlyGradientFrame_Disturbed, taxa_are_rows = FALSE)

# This seems to work, but View() fails on it....
phyloseq_sampleData <- sample_data(subsetGradientEnvAndSoilChemValues_abbrv_Disturbed)

phyloseq_combined <- phyloseq(phyloseq_OTUTable, phyloseq_sampleData)

# Design '~ 1' is valid for not specifying a design in phyloseq_to_deseq2

disturbedSample_combinedDataForDeseq2 <- phyloseq_to_deseq2(phyloseq_combined, ~ geogrGroupWE)

# This tests the log2 fold change of the number of seqs in remnant compared to disturbed sites for each OTU (+ val = higher
# in west; - val = higher in east sites), along with p-values for the log2 fold change.
disturbedSample_DEseqOutput <- DESeq(disturbedSample_combinedDataForDeseq2, fitType = 'local')

# access the raw results
DESeq_disturbed_WE_Compare_rawResults <- results(disturbedSample_DEseqOutput)

# sort results by p-value 
DESeq_disturbed_WE_Compare_pValueSorted <- DESeq_disturbed_WE_Compare_rawResults[order(DESeq_disturbed_WE_Compare_rawResults$padj),]

# sort results by effect size. Commented out on purpose - not needed here.
# DESeq_disturbed_WE_Compare_effectSizeSorted <- DESeq_disturbed_WE_Compare_rawResults[order(abs(DESeq_disturbed_WE_Compare_rawResults$log2FoldChange), decreasing = TRUE),]

# sort the results in phylogenetic order (using the order of entries of the genus attributions because those
# are listed in phylogenetic order e.g. phylogenetically closely related genera are listed immediately before/after one another)
DESeq_disturbed_WE_Compare_phyloSorted <- DESeq_disturbed_WE_Compare_rawResults[OTUTableIndicesInPhylogOrderForLabels,]

save(list = c("disturbedSample_DEseqOutput","DESeq_disturbed_WE_Compare_rawResults",
              "DESeq_disturbed_WE_Compare_pValueSorted",
              "DESeq_disturbed_WE_Compare_phyloSorted"), file = outputFile)

# Make the bar graph (Figure 3B)
# Ordered by phylogeny
# Leave width at 0.835 for the segments to fit well.

assignPointColors <- function(valueRow){
    print(valueRow[6])
    
    # adj. p value
    if(!is.na(valueRow[6]) & valueRow[6] < 0.05){
        # Log-2 fold change
        colorForBar <- ifelse(valueRow[2] > 0, "#984EA3", "#4DAF4A")
    } else{
        colorForBar <- "white"    
    }
    return(colorForBar)
} 

# Need to assign a color for each point (sample) based on its disturbance history.
colorHolder <- apply(DESeq_disturbed_WE_Compare_phyloSorted, MARGIN = 1, FUN = assignPointColors)

# Ordered by phylogeny
# Leave width at 0.835 for the segments to fit well.
barplot(DESeq_disturbed_WE_Compare_phyloSorted$log2FoldChange, col = colorHolder,
        ylab = 'Log 2 fold diff West (+) / East (-)', main = 'Disturbed Sites Only\n ordered by phylogeny',
        ylim = c(-6,6), width = 0.835)



#-----------------








# #-----------------
# # This is testing for East/West differences FOR ONLY Disturbed sites; Based on how DESeq2 requires balanced cross designs
# # for when have 2 diff factor vectors, I can't control for site location in a reasonable way within East/West; site location
# # is only controlled as a group within East sites and as a group within West sites.
# 
# pathToSaveDESeq2Output <- "~/Box Sync/Gradient_Coordianted_2013_2015/Graphs/DESeq2_aridityGradientOTUChangesGraphs/"
# outputFileName <- "GradientCoordinated_Module4B_DESeq2_WestEastEffect_onlyDisturbedSites_revised_072317.r"
# 
# gradientSamplesOnlyGradientFrame_noMissing_Disturbed <- gradientSamplesOnlyGradientFrame_noMissing[distGroup == 'Disturbed',]
# subsetGradientEnvAndSoilChemValues_abbrv_Disturbed <- subsetGradientEnvAndSoilChemValues_abbrv[distGroup == 'Disturbed',]
# 
# # Convert OTU table to phyloseq OTU table
# phyloseq_OTUTable <- otu_table(gradientSamplesOnlyGradientFrame_noMissing_Disturbed, taxa_are_rows = FALSE)
# 
# # This seems to work, but View() fails on it....
# phyloseq_sampleData <- sample_data(subsetGradientEnvAndSoilChemValues_abbrv_Disturbed)
# 
# phyloseq_combined <- phyloseq(phyloseq_OTUTable, phyloseq_sampleData)
# 
# # Design '~ 1' is valid for not specifying a design in phyloseq_to_deseq2
# 
# OTU_sampleDataForDESeqForCompareWithRF <- phyloseq_to_deseq2(phyloseq_combined, ~ geogrGroupWE)
# 
# # This tests the log2 fold change of the number of seqs in remnant compared to disturbed sites for each OTU (+ val = higher
# # in west; - val = higher in east sites), along with p-values for the log2 fold change.
# DESeqForCompareWithRFOutput <- DESeq(OTU_sampleDataForDESeqForCompareWithRF, fitType = 'local')
# 
# DESeqForCompareWithRFOutput_results_DisturbedWECompare <- results(DESeqForCompareWithRFOutput)
# DESeqForCompareWithRFOutput_resultsPValSorted_DisturbedWECompare <- DESeqForCompareWithRFOutput_results_DisturbedWECompare[order(DESeqForCompareWithRFOutput_results_DisturbedWECompare$padj),]
# DESeqForCompareWithRFOutput_resultsEffectSizeSorted_DisturbedWECompare <- DESeqForCompareWithRFOutput_results_DisturbedWECompare[order(abs(DESeqForCompareWithRFOutput_results_DisturbedWECompare$log2FoldChange), decreasing = TRUE),]
# DESeqForCompareWithRFOutput_resultsLabelSorted_DisturbedWECompare <- DESeqForCompareWithRFOutput_results_DisturbedWECompare[OTUTableIndicesInPhylogOrderForLabels,]
# 
# # This is importing the Random Forest results (taken from the DESeq2 script)
# # Import the RF binary OTU results for Remnant Sites only to sort the DESeq2 results in the order of RF most importance
# RFFileName <- "~/Box Sync/Gradient_Coordianted_2013_2015/Graphs/RandomForestGraphs/GradientCoordinated_RF_allSamples/GradientCoordinated_Module_4B_revisedDataset_072017/ GradientCoordinated_Module4B_BinaryPresAbs_5FoldCV_SiteEWCompare_DisturbedSamplesOnly_RandomForestGraphs_genusLabelled_072017_revised_ 2000Iter_RFData.r"
# # loads sortedAllCV_meanDecreaseAcc_byOTU, sortedAllCV_meanDecreaseGini_byOTU
# load(RFFileName)
# 
# # Match is a great function to find the indices of the OTUs in the gradientFrame that match the
# # order of the RF binary OTU output for OTUs with the greatest effect on decreased accuracy.
# OTUTableIndicesInRFBinaryDecAccOutputOrderForLabels <- match(names(sortedAllCV_meanDecreaseAcc_byOTU), colnames(gradientSamplesOnlyGradientFrame_noMissing))
# 
# DESeqForCompareWithRFOutput_resultsRFBinaryDecAccOTUSorted_DisturbedWECompare <- DESeqForCompareWithRFOutput_results_DisturbedWECompare[OTUTableIndicesInRFBinaryDecAccOutputOrderForLabels,]
# 
# outputFile <- paste(pathToSaveDESeq2Output, "/", outputFileName, sep = "")
# 
# save(list = c("DESeqForCompareWithRFOutput",
#               "DESeqForCompareWithRFOutput_results_DisturbedWECompare",
#               "DESeqForCompareWithRFOutput_resultsPValSorted_DisturbedWECompare",
#               "DESeqForCompareWithRFOutput_resultsEffectSizeSorted_DisturbedWECompare",
#               "DESeqForCompareWithRFOutput_resultsLabelSorted_DisturbedWECompare",
#               "DESeqForCompareWithRFOutput_resultsRFBinaryDecAccOTUSorted_DisturbedWECompare"),
#      file = outputFile)
# 
# if(plotPhyloOnly == 1){
# 
#   # Ordered by phylogeny
#   # Leave width at 0.835 for the segments to fit well.
#   phyloAdjPValColor <- ifelse(DESeqForCompareWithRFOutput_resultsLabelSorted_DisturbedWECompare$padj < 0.05, 'green', 'white')
#   barplot(DESeqForCompareWithRFOutput_resultsLabelSorted_DisturbedWECompare$log2FoldChange, col = phyloAdjPValColor,
#           ylab = 'Log 2 fold diff West (+) /East (-)', main = 'Disturbed Sites Only\n ordered by phylogeny',
#           ylim = c(-6,6), width = 0.835)
# 
#   # The start of the bars is at x = 0.5, and each bar is 1 segment unit wide
#   # when their width is 0.835.
# 
# 
# } else{
# 
#   par(mfrow = c(3,1))
# 
# 
#   # Ordered by phylogeny
#   phyloAdjPValColor <- ifelse(DESeqForCompareWithRFOutput_resultsLabelSorted_DisturbedWECompare$padj < 0.05, 'green', 'white')
#   barplot(DESeqForCompareWithRFOutput_resultsLabelSorted_DisturbedWECompare$log2FoldChange, col = phyloAdjPValColor,
#           ylab = 'Log 2 fold diff West (+) /East (-)', main = 'Disturbed Sites Only\n ordered by phylogeny',
#           ylim = c(-6,6))
# 
#   adjPValColor <- ifelse(DESeqForCompareWithRFOutput_resultsPValSorted_DisturbedWECompare$padj < 0.05, 'red', 'white')
#   # Ordered by p-Value
#   barplot(DESeqForCompareWithRFOutput_resultsPValSorted_DisturbedWECompare$log2FoldChange, col = adjPValColor,
#           ylab = 'Log 2 fold diff West (+) /East (-)', main = 'Disturbed Sites Only\n ordered by adj. p-value',
#           ylim = c(-6,6))
# 
#   # Ordered by RF binary OTUs with highest to lowest decrease in prediction accuracy.
#   RFadjPValColor <- ifelse(DESeqForCompareWithRFOutput_resultsRFBinaryDecAccOTUSorted_DisturbedWECompare$padj < 0.05, 'blue', 'white')
#   barplot(DESeqForCompareWithRFOutput_resultsRFBinaryDecAccOTUSorted_DisturbedWECompare$log2FoldChange, col = RFadjPValColor,
#           ylab = 'Log 2 fold diff West (+) /East (-)', main = 'Disturbed Sites Only\n ordered by Random Forest Binary Decrease Accuracy OTUs',
#           ylim = c(-6,6))
# 
#   # make a dashed vertical line on the sorted by RF decrease accuracy where the decrease in accuracy goes from
#   # positive (ie helpful OTUs in the model) to zero or negative (ie OTUs that either do not matter in the model
#   # or are actively hindering the prediction accuracy of the model)
# 
#   posRFAccVals <- which(sortedAllCV_meanDecreaseAcc_byOTU > 0)
# 
#   abline(v = posRFAccVals[length(posRFAccVals)], lty = 2, col = "gray")
# }
# #-----------------


# barplot(DESeqForCompareWithRFOutput_resultsEffectSizeSorted_DisturbedWECompare$log2FoldChange, names.arg = rownames(DESeqForCompareWithRFOutput_resultsEffectSizeSorted_DisturbedWECompare)
# color <- ifelse(DESeqForCompareWithRFOutput_resultsPValSorted_DisturbedWECompare$padj < 0.05, 'red', 'white')
# barplot(DESeqForCompareWithRFOutput_resultsPValSorted_DisturbedWECompare$log2FoldChange, col = color)


#=====================
#=====================

# # -----------------------
# # This is to run rLog to transform the OTU table for Random Forests using the West/Central/East site factor as
# # the constraint on the tranformation model
# 
# # Convert OTU table to phyloseq OTU table
# phyloseq_OTUTable_forRLog <- otu_table(gradientSamplesOnlyGradientFrame_noMissing, taxa_are_rows = FALSE)
# 
# # This seems to work, but View() fails on it....
# phyloseq_sampleData_forRLog <- sample_data(subsetGradientEnvAndSoilChemValues_abbrv)
# 
# phyloseq_combined_forRLog <- phyloseq(phyloseq_OTUTable_forRLog, phyloseq_sampleData_forRLog)
# 
# 
# # Design '~ 1' is valid for not specifying a design in phyloseq_to_deseq2 
# 
# OTU_sampleData_forRLog <- phyloseq_to_deseq2(phyloseq_combined_forRLog, ~ aridityGroup)
# 
# rlogTransformation <- rlog(OTU_sampleData_forRLog, fitType = 'local', blind = FALSE)
# 
# #transformedRLog_OTUTable <- assay(rlogTest)
# transformedRLog_OTUTable_forRF <- t(assay(rlogTransformation))
# 
# save(transformedRLog_OTUTable_forRF, file = "~/Box Sync/Gradient_Coordianted_2013_2015/GradientCoordinated_DESeq2_rLogTransformedOTUTable_using3AridityGroups_module4B_noKZDist3_ForRandomForest_112816.r")
# # --------------------

# #------------
# # This is for OTU-specific diagnostics to determine +/- in remanant or disturbed sites for the random forest run
# 
# # Convert OTU table to phyloseq OTU table
# phyloseq_OTUTable <- otu_table(gradientSamplesOnlyGradientFrame_noMissing, taxa_are_rows = FALSE)
# 
# # This seems to work, but View() fails on it....
# phyloseq_sampleData <- sample_data(subsetGradientEnvAndSoilChemValues_abbrv)
# 
# phyloseq_combined <- phyloseq(phyloseq_OTUTable, phyloseq_sampleData)
# 
# 
# # Design '~ 1' is valid for not specifying a design in phyloseq_to_deseq2 
# 
# OTU_sampleDataForDESeqForCompareWithRF <- phyloseq_to_deseq2(phyloseq_combined, ~ aridityGroup + Disturbed_Remnant)
# 
# # This tests the log2 fold change of the number of seqs in remnant compared to disturbed sites for each OTU (+ val = higher
# # in remnant; - val = higher in disturbed sites), along with p-values for the log2 fold change. The highest scoring OTU in the random
# # forests based on pres/abs is also the highest here (OTU110), and the direction (remn/dist) of each OTU found in the random forest
# # pres/abs is almost always the same with the log2 fold change direction found here. 
# DESeqForCompareWithRFOutput <- DESeq(OTU_sampleDataForDESeqForCompareWithRF, fitType = 'local')
# 
# DESeqForCompareWithRFOutput_results <- results(DESeqForCompareWithRFOutput)
# DESeqForCompareWithRFOutput_resultsPValSorted <- DESeqForCompareWithRFOutput_results[order(DESeqForCompareWithRFOutput_results$padj),]
# DESeqForCompareWithRFOutput_resultsEffectSizeSorted <- DESeqForCompareWithRFOutput_results[order(abs(DESeqForCompareWithRFOutput_results$log2FoldChange), decreasing = TRUE),]
# DESeqForCompareWithRFOutput_resultslabelSorted <- DESeqForCompareWithRFOutput_results[OTUTableIndicesInPhylogOrderForLabels,]
# 
# pValForGraph_CompareWithRFOutput <- DESeqForCompareWithRFOutput_resultslabelSorted$padj
# 
# colForGraph_CompareWithRFOutput <- ifelse(pValForGraph_CompareWithRFOutput <= 0.05, rgb(1,0,0), rgb(1,1,1))
# 
# barplot(DESeqForCompareWithRFOutput_resultslabelSorted$log2FoldChange, col= colForGraph_CompareWithRFOutput, 
#         ylab = "Log fold change Remnant/Disturbed", main = "Relative abundance of OTUs \nin remnant compared to disturbed sites", ylim = c(-6, 6))
# 
# 
# # -----------------

# # -----------------
# # To look at consistent differences in OTU abundance across the aridity gradient (as represented by the aridityGroup of Western (Klemme, Stillwater, Hays), 
# # Central (Konza, Ft. Riley, Rockefeller, Welda), and Eastern (all MO and IL) sites) in REMNANT sites only
# 
# CWColor <- "#E41A1C"
# ECColor <- "#377EB8"
# EWColor <- "#4DAF4A"
# 
# gradientSamplesOnlyGradientFrame_noMissing_remnantOnly <- gradientSamplesOnlyGradientFrame_noMissing[distGroup == "Remnant",]
# subsetGradientEnvAndSoilChemValues_abbrv_remnantOnly <- subsetGradientEnvAndSoilChemValues_abbrv[distGroup == "Remnant",]
# 
# # Convert OTU table to phyloseq OTU table
# phyloseq_OTUTable_remnantOnly <- otu_table(gradientSamplesOnlyGradientFrame_noMissing_remnantOnly, taxa_are_rows = FALSE)
# 
# # This seems to work, but View() fails on it....
# phyloseq_sampleData_remnantOnly <- sample_data(subsetGradientEnvAndSoilChemValues_abbrv_remnantOnly)
# 
# phyloseq_combined_remnantOnly <- phyloseq(phyloseq_OTUTable_remnantOnly, phyloseq_sampleData_remnantOnly)
# 
# 
# # Design '~ 1' is valid for not specifying a design in phyloseq_to_deseq2 
# 
# OTU_sampleDataForDESeqForRemnAridChanges <- phyloseq_to_deseq2(phyloseq_combined_remnantOnly, ~ aridityGroup)
# 
# # This tests the log2 fold change of the number of seqs 
# DESeqForRemnAridChanges <- DESeq(OTU_sampleDataForDESeqForRemnAridChanges, fitType = 'local')
# 
# # The numerator for these tests is always the east-most site of the pair, and the denominator is the west-most
# # Testing the Central/West comparison (use print() to verify the contrasts are correct)
# DESeqForRemnAridChanges_results_CW <- results(DESeqForRemnAridChanges, altHypothesis = 'greaterAbs', contrast = c("aridityGroup","C","W"))
# DESeqForRemnAridChanges_results_CW_PValSorted <- DESeqForRemnAridChanges_results_CW[order(DESeqForRemnAridChanges_results_CW$padj),]
# DESeqForRemnAridChanges_results_CW_labelSorted <- DESeqForRemnAridChanges_results_CW[OTUTableIndicesInPhylogOrderForLabels,]
# 
# pValForGraph_CW <- DESeqForRemnAridChanges_results_CW_labelSorted$padj
# colForGraph_CW <- ifelse(pValForGraph_CW <= 0.05, CWColor, rgb(1,1,1))
# barplot(DESeqForRemnAridChanges_results_CW_labelSorted$log2FoldChange, col= colForGraph_CW, 
#         ylab = "Log fold change Central/West", main = "Remnant sites", ylim = c(-6, 6))
# 
# 
# # Testing the East/Central comparison (use print() to verify the contrasts are correct)
# DESeqForRemnAridChanges_results_EC <- results(DESeqForRemnAridChanges, altHypothesis = 'greaterAbs', contrast = c("aridityGroup","E","C"))
# DESeqForRemnAridChanges_results_EC_PValSorted <- DESeqForRemnAridChanges_results_EC[order(DESeqForRemnAridChanges_results_EC$padj),]
# DESeqForRemnAridChanges_results_EC_labelSorted <- DESeqForRemnAridChanges_results_EC[OTUTableIndicesInPhylogOrderForLabels,]
# 
# pValForGraph_EC <- DESeqForRemnAridChanges_results_EC_labelSorted$padj
# colForGraph_EC <- ifelse(pValForGraph_EC <= 0.05, ECColor, rgb(1,1,1))
# barplot(DESeqForRemnAridChanges_results_EC_labelSorted$log2FoldChange, col= colForGraph_EC, 
#         ylab = "Log fold change East/Central", main = "Remnant sites", ylim = c(-6, 6))
# 
# 
# # Testing the Central/West comparison (use print() to verify the contrasts are correct)
# DESeqForRemnAridChanges_results_EW <- results(DESeqForRemnAridChanges, altHypothesis = 'greaterAbs', contrast = c("aridityGroup","E","W"))
# DESeqForRemnAridChanges_results_EW_PValSorted <- DESeqForRemnAridChanges_results_EW[order(DESeqForRemnAridChanges_results_EW$padj),]
# DESeqForRemnAridChanges_results_EW_labelSorted <- DESeqForRemnAridChanges_results_EW[OTUTableIndicesInPhylogOrderForLabels,]
# 
# pValForGraph_EW <- DESeqForRemnAridChanges_results_EW_labelSorted$padj
# colForGraph_EW <- ifelse(pValForGraph_EW <= 0.05, EWColor, rgb(1,1,1))
# barplot(DESeqForRemnAridChanges_results_EW_labelSorted$log2FoldChange, col= colForGraph_EW, 
#         ylab = "Log fold change East/West", main = "Remnant sites", ylim = c(-6, 6))
# #axis(1, at = seq(0, length(DESeqForRemnAridChanges_results_EW_labelSorted$padj),1))
# #rug(jitter(seq(0,length(DESeqForRemnAridChanges_results_EW_labelSorted$padj)),1))
# 
# # ---------------------

# # -----------------
# # To look at consistent differences in OTU abundance across the aridity gradient (as represented by the aridityGroup of Western (Klemme, Stillwater, Hays), 
# # Central (Konza, Ft. Riley, Rockefeller, Welda), and Eastern (all MO and IL) sites) in DISTURBED sites only
# 
# CWColor <- "#E41A1C"
# ECColor <- "#377EB8"
# EWColor <- "#4DAF4A"
# 
# gradientSamplesOnlyGradientFrame_noMissing_disturbedOnly <- gradientSamplesOnlyGradientFrame_noMissing[distGroup == "Disturbed",]
# subsetGradientEnvAndSoilChemValues_abbrv_disturbedOnly <- subsetGradientEnvAndSoilChemValues_abbrv[distGroup == "Disturbed",]
# 
# # Convert OTU table to phyloseq OTU table
# phyloseq_OTUTable_disturbedOnly <- otu_table(gradientSamplesOnlyGradientFrame_noMissing_disturbedOnly, taxa_are_rows = FALSE)
# 
# # This seems to work, but View() fails on it....
# phyloseq_sampleData_disturbedOnly <- sample_data(subsetGradientEnvAndSoilChemValues_abbrv_disturbedOnly)
# 
# phyloseq_combined_disturbedOnly <- phyloseq(phyloseq_OTUTable_disturbedOnly, phyloseq_sampleData_disturbedOnly)
# 
# 
# # Design '~ 1' is valid for not specifying a design in phyloseq_to_deseq2 
# 
# OTU_sampleDataForDESeqForDistAridChanges <- phyloseq_to_deseq2(phyloseq_combined_disturbedOnly, ~ aridityGroup)
# 
# # This tests the log2 fold change of the number of seqs 
# DESeqForDistAridChanges <- DESeq(OTU_sampleDataForDESeqForDistAridChanges, fitType = 'local')
# 
# # The numerator for these tests is always the east-most site of the pair, and the denominator is the west-most
# # Testing the Central/West comparison (use print() to verify the contrasts are correct)
# DESeqForDistAridChanges_results_CW <- results(DESeqForDistAridChanges, altHypothesis = 'greaterAbs', contrast = c("aridityGroup","C","W"))
# DESeqForDistAridChanges_results_CW_PValSorted <- DESeqForDistAridChanges_results_CW[order(DESeqForDistAridChanges_results_CW$padj),]
# DESeqForDistAridChanges_results_CW_labelSorted <- DESeqForDistAridChanges_results_CW[OTUTableIndicesInPhylogOrderForLabels,]
# 
# pValForGraph_CW <- DESeqForDistAridChanges_results_CW_labelSorted$padj
# colForGraph_CW <- ifelse(pValForGraph_CW <= 0.05, CWColor, rgb(1,1,1))
# barplot(DESeqForDistAridChanges_results_CW_labelSorted$log2FoldChange, col= colForGraph_CW, 
#         ylab = "Log fold change Central/West", main = "Disturbed sites", ylim = c(-6,6))
# 
# 
# # Testing the East/Central comparison (use print() to verify the contrasts are correct)
# DESeqForDistAridChanges_results_EC <- results(DESeqForDistAridChanges, altHypothesis = 'greaterAbs', contrast = c("aridityGroup","E","C"))
# DESeqForDistAridChanges_results_EC_PValSorted <- DESeqForDistAridChanges_results_EC[order(DESeqForDistAridChanges_results_EC$padj),]
# DESeqForDistAridChanges_results_EC_labelSorted <- DESeqForDistAridChanges_results_EC[OTUTableIndicesInPhylogOrderForLabels,]
# 
# pValForGraph_EC <- DESeqForDistAridChanges_results_EC_labelSorted$padj
# colForGraph_EC <- ifelse(pValForGraph_EC <= 0.05, ECColor, rgb(1,1,1))
# barplot(DESeqForDistAridChanges_results_EC_labelSorted$log2FoldChange, col= colForGraph_EC, 
#         ylab = "Log fold change East/Central", main = "Disturbed sites", ylim = c(-6,6))
# 
# 
# # Testing the Central/West comparison (use print() to verify the contrasts are correct)
# DESeqForDistAridChanges_results_EW <- results(DESeqForDistAridChanges, altHypothesis = 'greaterAbs', contrast = c("aridityGroup","E","W"))
# DESeqForDistAridChanges_results_EW_PValSorted <- DESeqForDistAridChanges_results_EW[order(DESeqForDistAridChanges_results_EW$padj),]
# DESeqForDistAridChanges_results_EW_labelSorted <- DESeqForDistAridChanges_results_EW[OTUTableIndicesInPhylogOrderForLabels,]
# 
# pValForGraph_EW <- DESeqForDistAridChanges_results_EW_labelSorted$padj
# colForGraph_EW <- ifelse(pValForGraph_EW <= 0.05, EWColor, rgb(1,1,1))
# barplot(DESeqForDistAridChanges_results_EW_labelSorted$log2FoldChange, col= colForGraph_EW, 
#         ylab = "Log fold change East/West", main = "Disturbed sites", ylim = c(-6,6))
# #axis(1, at = seq(0, length(DESeqForDistAridChanges_results_EW_labelSorted$padj),1))
# #rug(jitter(seq(0,length(DESeqForDistAridChanges_results_EW_labelSorted$padj)),1))
# 
# # ---------------------



# ---------------
# This is to add the genus/demarking color line and dashed lines to the ordered by
# phylo bars (regardless of the data source) when they are plotted by themselves. This
# includes flexible tuning parameters for getting the widths right.

# This is now set up for the revised OTU attributions with Mor elo outgroup (072117)

lineThickness <- 3
abLineColor <- "gray"
abLineThickness <- 0.75
segmentBarWidth <- 1
# Keep at 0.3 and save as 16"W by 4" H
barOffset <- 0.3

contrastBarColor1 <- "black"
contrastBarColor2 <- "gray"

# For testing
#yLocation <- 0

# For real
yLocation = -5.9

# For Div. (11)
segments(x0 = barOffset, 
         y0 = yLocation, 
         x1 = barOffset + 11 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor1, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 11 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Aca. (15)
segments(x0 = barOffset + 11 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 26 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor2, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 26 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Scu. (2)
segments(x0 = barOffset + 26 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 28 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor1, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 28 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Cet. (1)
segments(x0 = barOffset + 28 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 29 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor2, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 29 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Gig. (2)
segments(x0 = barOffset + 29 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 31 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor1, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 31 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Cla. (14)
segments(x0 = barOffset + 31 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 45 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor2, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 45 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Glomeraceae. (57)
segments(x0 = barOffset + 45 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 102 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor1, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 102 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Glomus. (21)
segments(x0 = barOffset + 102 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 123 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor2, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 123 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Fun. (3)
segments(x0 = barOffset + 123 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 126 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor1, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 126 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Sep. (12)
segments(x0 = barOffset + 126 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 138 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor2, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 138 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Rhi. (55)
segments(x0 = barOffset + 138 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 193 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor1, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 193 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Amb. (2)
segments(x0 = barOffset + 193 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 195 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor2, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 195 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Arc. (1)
segments(x0 = barOffset + 195 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 196 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor1, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 196 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Par. (3)
segments(x0 = barOffset + 196 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 199 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor2, lwd = lineThickness, lend = "butt")
 


# ------------------------
