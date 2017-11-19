# Combines DESEQ2 differential abundance results with Random Forest classifier results to determine
# which OTUs are promising indicators of remnant prairie sites, or indicators of remnant prairies from
# different sides of the precipitation gradient (West/East). This produces all 6 subplots used in Figure 5. Note - the
# plots for the Random Forest and the DESEQ2 results are produced separately for each data subset (need to scroll through 
# them). The plot
# from the DESEQ2 results includes the full genus attribution and the OTU number, but these are only visible
# after saving as a pdf and then opening the pdf in a program like Illustrator for re-sizing, as was done to make
# Figure 5.

# Working correctly 080717.

library("DESeq2")

# Specify the number of rows to subset out of the RF and the DESEQ2 output (this will be the number of OTUs displayed in the graph)
# Use 199 to show all.
numSubsetRows <- 20

# Import OTU names from phylo attribution using the phylogeny of only AM fungal OTU seqs with refs
# This is revised with Mor elo outgroup 072117
OTUNames_fromPhylo <- read.table("~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/GradientCoordinated_GenusAttribOfOTUs_fromUSEARCHRefSeqs_WITHPhyloSuppl_Mor_elo_outgroup_july_2017_forR.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)


# ========================

producePlots <- function(coloring = NA){
    
    OTUNames_fromRF <- names(sortedAllCV_meanDecreaseAcc_byOTU)
    
    OTUNames_fromDESeq2 <- row.names(DESeqResultsMatrix)
    
    # Match is a great function to find the indices of the OTUs in the DESeq2 results that match the order of the OTUs in the RF results (the RF results take
    # precedence here in terms of keeping the order of the OTUs wrt the greatest effect in the RFs.)
    OTU_indicesInDESeq2Results_matchedToRFResults <- match(OTUNames_fromRF, OTUNames_fromDESeq2)
    
    # Convert the OTU names from eg "X110" to 110 in preparation to compare to the OTU numbers with the associated genus names.
    OTUNames_fromRF_asNumeric <- sapply(OTUNames_fromRF, function(x){as.numeric(substr(x,2,nchar(x)))})
    
    OTU_indicesInGenusNames_matchedToRFResults <- match(OTUNames_fromRF_asNumeric, OTUNames_fromPhylo$OTU_num)
    
    # Now re-order the rows in the phylogeny-based genera attributions to match the OTU name order of the RF results.
    matchedToRFOTUs_Phylo_OTUs <- OTUNames_fromPhylo[OTU_indicesInGenusNames_matchedToRFResults,]
    
    # Now re-order the rows in the DESeq2 results to match the OTU name order of the RF results.
    matchedToRFOTUs_DESeq2Results <- DESeqResultsMatrix[OTU_indicesInDESeq2Results_matchedToRFResults,]
    
    # Keep only the top ## rows (for ease of use with plotting)
    matchedToRFOTUs_DESeq2Results_highestSubset <- matchedToRFOTUs_DESeq2Results[1:numSubsetRows,]
    
    sortedAllCV_meanDecreaseAcc_byOTU_highestSubset <- sortedAllCV_meanDecreaseAcc_byOTU[1:numSubsetRows]
    
    # Set the proper colors to use.
    if(coloring == "siteHist"){
        matchedToRFOTUs_DESeq2Results_colVect <- ifelse(matchedToRFOTUs_DESeq2Results_highestSubset$log2FoldChange > 0, "#FF7F00", "#377EB8")
    } else if(coloring == "Precip"){
        matchedToRFOTUs_DESeq2Results_colVect <- ifelse(matchedToRFOTUs_DESeq2Results_highestSubset$log2FoldChange > 0, "#4DAF4A", "#984EA3")
    }
    
    # If the adjusted p value for an OTU is NA, then change the color in the colVect to gray regardless of the differential representation in the log2 fold column.
    matchedToRFOTUs_DESeq2Results_colVect[is.na(matchedToRFOTUs_DESeq2Results_highestSubset$padj)] <- "gray50"
    
    # Set the shape (which is also the fill) for solid if the padj is < 0.05, and hollow if not.
    matchedToRFOTUs_DESeq2Results_pchVect <- ifelse(matchedToRFOTUs_DESeq2Results_highestSubset$padj < 0.05, 17, 24) 
    
    # clean up for if there were padj values of NA, then add entries to those indices of the pch vector
    matchedToRFOTUs_DESeq2Results_pchVect[is.na(matchedToRFOTUs_DESeq2Results_pchVect)] <- 21
    
    logFoldChangeVector <- matchedToRFOTUs_DESeq2Results_highestSubset$log2FoldChange
    
    # Add the OTU# (without the leading X) to the genus names
    
    RawOTUNames <- row.names(matchedToRFOTUs_DESeq2Results_highestSubset)
    ProcessedOTUNames <- sapply(RawOTUNames, function(x){processed <- substr(x,2,nchar(x))})
    
    OTUNumAddedOTUNames <- paste(matchedToRFOTUs_Phylo_OTUs$Genus_Attribution[1:numSubsetRows], "OTU", ProcessedOTUNames)
    
    logFoldChangeVector <- setNames(logFoldChangeVector, OTUNumAddedOTUNames)

    
    # -------------------
    # DESeq2 Dotchart
    
    plot(logFoldChangeVector, seq(length(logFoldChangeVector),1,-1), 
         col = "white", 
         yaxt ="n", xlab = "Log 2 change in OTU \nrepresentation",
         main = plotTitle, ylab = "", cex.lab = 1.2, xlim = c(-6,6))
    
    mtext(text = OTUNumAddedOTUNames, side = 2, 
          at = seq(length(logFoldChangeVector), 1, -1), las = 2, line = 0.3, cex = 1.3)
    
    x <- sapply(seq(1,length((logFoldChangeVector))), function(x){abline(h = x, col = "gray", lty = 3)})
    
    abline(v = 0, lty = 2, col = "black")
    
    points(logFoldChangeVector, seq(length(logFoldChangeVector),1,-1), 
           col = matchedToRFOTUs_DESeq2Results_colVect, 
           pch = matchedToRFOTUs_DESeq2Results_pchVect, cex = 2.5)
    
    # -------------------
    
    # -------------------
    # This is the dotchart from the RF
    
    plot(sortedAllCV_meanDecreaseAcc_byOTU_highestSubset, seq(length(logFoldChangeVector),1,-1), 
         col = "white", 
         yaxt ="n", xlab = "Mean decrease in accuracy\n across all CV folds", ylab = "", cex.lab = 1.2, xlim = c(2,18), xaxp = c(2,18,4))
    
    x <- sapply(seq(1,length((logFoldChangeVector))), function(x){abline(h = x, col = "gray", lty = 3)})
    
    abline(v = 0, lty = 2, col = "black")
    
    points(sortedAllCV_meanDecreaseAcc_byOTU_highestSubset, seq(length(logFoldChangeVector),1,-1), 
           col = "black", 
           pch = 19, cex = 2.5)
    
}

# ========================

# Select which scenario to data import

#-------------------------------------
# This is for all sites remn/dist

# This is importing the Random Forest results (taken from the DESeq2 script)
# Import the RF binary OTU results for allSites to sort the DESeq2 results in the order of RF most importance

pathToRFOutput <- "~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/autoGeneratedScriptOutput/GradientCoordinated_allSamples_10CV_RandomForest_2000Iter_RFData.r"

# loads sortedAllCV_meanDecreaseAcc_byOTU, sortedAllCV_meanDecreaseGini_byOTU
load(pathToRFOutput)

# loads the DESeq2 results
pathToSaveDESeq2Output <- "~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/autoGeneratedScriptOutput/GradientCoordinated_DESEQ2_results_allSamples.r"

load(pathToSaveDESeq2Output)

# This is still a structure in class DESeq2, but when I subset it into vectors using dataframe indexing, those
# vectors are class numeric as expected
DESeqResultsMatrix <- DESeqResults(DESeq_distHistCompare_pValueSorted)

plotTitle <- "All sites\nL: Dist; R: Remn"

# Call the producePlots function to make the 2 plots.
producePlots(coloring = "siteHist")

#------------------------------------

#-------------------------------------
# This is for East/West differences FOR ONLY Remnant sites

# This is importing the Random Forest results (taken from the DESeq2 script)
# Import the RF binary OTU results for Remnant Sites only to sort the DESeq2 results in the order of RF most importance
# RFFileName <- "~/Box Sync/Gradient_Coordianted_2013_2015/Graphs/RandomForestGraphs/GradientCoordinated_RF_allSamples/GradientCoordinated_RF_usingModule4B_BinaryPresAbs_OTUTable_120616/GradientCoordinated_Module4B_RFBinaryOTU_RemnSites_testEastWest/ GradientCoordinated_Module4B_BinaryPresAbs_5FoldCV_SiteEWCompare_RemnantSamplesOnly_RandomForestGraphs_120416_ 2000Iter_RFData.r"

RFFileName <- "~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/autoGeneratedScriptOutput/GradientCoordinated_RemnantSamples_5CV_RandomForest_2000Iter_RFData.r"

# loads sortedAllCV_meanDecreaseAcc_byOTU, sortedAllCV_meanDecreaseGini_byOTU
load(RFFileName)

# loads the DESeq2 results
DESeq2FileName <- "~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/autoGeneratedScriptOutput/GradientCoordinated_DESEQ2_results_remnantSites.r"

load(DESeq2FileName)

# This is still a structure in class DESeq2, but when I subset it into vectors using dataframe indexing, those
# vectors are class numeric as expected
DESeqResultsMatrix <- DESeqResults(DESeq_remnant_WE_Compare_pValueSorted)

# For the E/W differences, the log2 abundances are coded in the file with + for West, and - for East, but when plotting vertically with the dotplots,
# want - for West and + for East so the orientation is like a compass, so multiply the log2 values by -1 to flip the signs.
DESeqResultsMatrix$log2FoldChange <- DESeqResultsMatrix$log2FoldChange * -1

plotTitle <- "Remnant sites only\nLeft: West; Right: East"

# Call the producePlots function to make the 2 plots.
producePlots(coloring = "Precip")

#------------------------------------

#-------------------------------------
# This is for East/West differences FOR ONLY Disturbed sites

# This is importing the Random Forest results (taken from the DESeq2 script)
# Import the RF binary OTU results for Remnant Sites only to sort the DESeq2 results in the order of RF most importance

RFFileName <- "~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/autoGeneratedScriptOutput/GradientCoordinated_DisturbedSamples_5CV_RandomForest_2000Iter_RFData.r"

# loads sortedAllCV_meanDecreaseAcc_byOTU, sortedAllCV_meanDecreaseGini_byOTU
load(RFFileName)

# loads the DESeq2 results
DESeq2FileName <- "~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/autoGeneratedScriptOutput/GradientCoordinated_DESEQ2_results_disturbedSites.r"

load(DESeq2FileName)

# This is still a structure in class DESeq2, but when I subset it into vectors using dataframe indexing, those
# vectors are class numeric as expected
DESeqResultsMatrix <- DESeqResults(DESeq_disturbed_WE_Compare_pValueSorted)

# For the E/W differences, the log2 abundances are coded in the file with + for West, and - for East, but when plotting vertically with the dotplots,
# want - for West and + for East so the orientation is like a compass, so multiply the log2 values by -1 to flip the signs.
DESeqResultsMatrix$log2FoldChange <- DESeqResultsMatrix$log2FoldChange * -1

plotTitle <- "Disturbed sites only\nLeft: West; Right: East"

# Call the producePlots function to make the 2 plots.
producePlots(coloring = "Precip")

#------------------------------------
