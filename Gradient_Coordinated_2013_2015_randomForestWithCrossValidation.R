# Random forest classifers for all samples (predicting site history), and for disturbed samples only or 
# remnant samples only (predicting West/East). This runs all 3. The mean accuracy values output here
# are shown in Table 1. The results are saved for use with the DESEQ2 results in producing Figure 5.

# Working correctly 080117

library(randomForest)
library(caret)
# Need this for the confusion matrix function of caret
library(e1071)

# Import the OTU table
OTUTableGradient <- read.table("~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/GradientCoordinated_OTUTable_forR.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Import the AM fungal OTU numbers (identified from the phylogeny with reference sequences, Appendix 1, Figure S1).
AMFOTUNumbers <- scan("~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/GradientCoordinated_AMF_OTUNumbers_forR.txt")

# Remove the labels and other information columns to keep only the OTU table counts
OTUframeGradient <- data.frame(OTUTableGradient[,5:length(OTUTableGradient[1,])])

# This works correctly to subset only the AMF OTU columns by moving them by name (eg X1,X12) 
subsetOTUframeGradient_onlyAMFOTUs <- OTUframeGradient[,paste("X",AMFOTUNumbers,sep="")]

# Add the sample names back to the otu table data frame as rownames
sampleNames <- OTUTableGradient[,4]
rownames(subsetOTUframeGradient_onlyAMFOTUs) <- sampleNames

# Import environmental and soil chemical values for each seq barcode
gradientEnvAndSoilChemValues <- read.table("~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/GradientCoordinated_environ_soilNutrient_values_forR.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

subsetGradientEnvAndSoilChemValues_abbrv <- gradientEnvAndSoilChemValues[,-1:-3]
rownames(subsetGradientEnvAndSoilChemValues_abbrv) <- sampleNames

# Import the phylogeny-based genus attributions to each of the OTUs (from Figure S1). These are listed in 
# phylogenetic order e.g. phylogenetically closely related genera are listed immediately before/after one another)
genusAttributions <- read.table("~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/GradientCoordinated_OTU_GenusAttributions_forR.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

OTUNamesMatchOTUTable <- paste("X", genusAttributions$OTU_num, sep = "")

# Match is a great function to find the indices of the OTUs in the gradientFrame that match the 
# order of the phylogenetically-sorted OTU names in the OTUNamesMatch. This will be used for sorting the 
# DESeq2 output in phylogenetic order to see if there are consistent differences across the aridity gradient.
OTUTableIndicesInPhylogOrderForLabels <- match(OTUNamesMatchOTUTable, colnames(subsetOTUframeGradient_onlyAMFOTUs))

genusNamesForOTUMatch <- genusAttributions$Genus_Attribution

distGroup <- subsetGradientEnvAndSoilChemValues_abbrv$Disturbed_Remnant

# For the revised analysis 120416
# East/West geographic groups.
# West: Klemme, Stillwater, Hays, Ft. Riley (2013), Konza (2015)
# East: Lawrence (Rockefeller and Welda), SW MO, Morris, IL (2013), IL (2015)
geogrGroupWE <- c(rep("W", 7), rep("E", 2), rep("W", 4), rep("E", 8), rep("W", 17), rep("E", 13),
                  rep("W", 4), rep("E", 17), rep("W", 20), rep("E", 3), rep("W", 6), rep("E", 7))
geogrGroupWE_DF <- as.data.frame(geogrGroupWE)

subsetGradientEnvAndSoilChemValues_abbrv$distGroup <- as.factor(distGroup)

subsetGradientEnvAndSoilChemValues_abbrv$geogrGroupWE <- as.factor(geogrGroupWE)

# This is the main function for running the random forests. It is called 3 times, 1 for each
# data partition, below.
runRandomForest <- function(factorForTest, OTUTableToUse, cvFolds, plotFilePath, fileNameStem, randomForestFormula){

    runRandomForestFold <- function(foldNumber, factorToTest, OTUTableToUse, plotFilePath, fileNameStem, randomForestFormula){
        
        # Use the '-' indexing to select all rows except the row indices in testFolds[[x]] list.
        
        trainingSet <- OTUTableToUse[-cvFolds[[foldNumber]],]
        testSet <- OTUTableToUse[cvFolds[[foldNumber]],]
        
        print(paste("The fold number is:", foldNumber))
        
        set.seed(623)
        binaryCodedOTU_RF <- randomForest(randomForestFormula, data = trainingSet,
                                          importance = TRUE, ntree = 2000)
        
        pdf(paste(plotFilePath, fileNameStem, as.character(foldNumber), "_varImportPlot.pdf", sep = ""))
        # Plot the importance of the different variables (OTUs) in the accuracy of the random forest
        varImpPlot(binaryCodedOTU_RF)
        dev.off()
        
        # matrix of the importance values (decrease in acc., and Gini index) for each OTU - these are
        # what are plotted in the varImpPlot above.
        rfVarImportanceMatrix <- randomForest::importance(binaryCodedOTU_RF)
        
        # plot the error rates with the number of trees (don't fully understand this, or the different lines)
        pdf(paste(plotFilePath, fileNameStem, as.character(foldNumber), "_errorPlot.pdf", sep = ""))
        plot(binaryCodedOTU_RF)
        dev.off()
        
        # The prediction can be returned in terms of the binary 0/1 classification (response) generated 
        # from the majority vote (this is what I will end up using), the probability
        # of each sample being in either 0/1 (prob) from the votes, or the votes themselves from each tree in the forest
        # (vote). It's not possible to specify that it returns all 3, so need to run the prediction 
        # 3 separate times.
        predictedLabelsFromRF_response <- predict(binaryCodedOTU_RF, testSet, type = "response") 
        predictedLabelsFromRF_prob <- predict(binaryCodedOTU_RF, testSet, type = "prob")
        predictedLabelsFromRF_votes <- predict(binaryCodedOTU_RF, testSet, type = "vote", norm.votes = FALSE)
        
        testSetSampleNames <- sampleNames[as.numeric(names(predictedLabelsFromRF_response))]
        # The prediction label var name is passed to this function as a char string, so need to evaluate
        # the string as the var it represents here.
        testSetTrueLabels <- eval(parse(text = factorToTest))[as.numeric(names(predictedLabelsFromRF_response))]
        
        rFPredictActualCompare <- list("siteName" = testSetSampleNames, "actualLabels" = testSetTrueLabels,
                                       "rFPredictedLabels" = as.character(predictedLabelsFromRF_response),
                                       "rFPredictedLabelVoteFractions" = predictedLabelsFromRF_prob,
                                       "rFPredictedLabelVoteCounts" = predictedLabelsFromRF_votes,
                                       "rFOTUImportanceMatrix" = rfVarImportanceMatrix)
        
        # The 95%CI for the accuracy in the confusionMatrix is computed with a 2-sided binomial test with the null
        # being the proportion of the training data made up of the largest partition (largest class - in this case 'remnant').
        # The p-value given is 1-sided that the accuracy is greater than the proportion of the training data made up of the largest partition (largest class - in this case 'remnant').
        confusionMat <- caret::confusionMatrix(as.factor(rFPredictActualCompare$rFPredictedLabels), as.factor(rFPredictActualCompare$actualLabels))
        rFPredictActualCompare$confusionMatrix <- confusionMat
        return(rFPredictActualCompare)
    }
    
    allCrossValidationResults <- lapply(seq(1,length(cvFolds),1), function(x){runRandomForestFold(x, factorForTest, OTUTableToUse, plotFilePath, fileNameStem, randomForestFormula)})
    
    allCrossValidation_decreaseAccVals <- sapply(allCrossValidationResults, function(x){x$rFOTUImportanceMatrix[,3]})
    allCrossValidation_decreaseGiniVals <- sapply(allCrossValidationResults, function(x){x$rFOTUImportanceMatrix[,4]})
    allCrossValidation_confuseMatrixStats <- sapply(allCrossValidationResults, function(x){x$confusionMatrix$overall})
    
    # # -------------
    # # Toggle here for labelling by OTU (comment this out) or labelling by genus (run this block)
    # # Match is a great function to find the indices of the OTUs in the gradientFrame that match the 
    # # order of the phylogenetically-sorted OTU names in the OTUNamesMatch. This is to label the 
    # # OTUs on the y axis with their genus names.
    # OTUNameMatchGini <- match(rownames(allCrossValidation_decreaseGiniVals), OTUNamesMatchOTUTable)
    # # Replace the OTU numbers with the genus labels from the USEARCH genus match file loaded at the start of the script.
    # rownames(allCrossValidation_decreaseGiniVals) <- genusNamesForOTUMatch[OTUNameMatchGini]
    # 
    # OTUNameMatchDecAcc <- match(rownames(allCrossValidation_decreaseAccVals), OTUNamesMatchOTUTable)
    # # Replace the OTU numbers with the genus labels from the USEARCH genus match file loaded at the start of the script.
    # rownames(allCrossValidation_decreaseAccVals) <- genusNamesForOTUMatch[OTUNameMatchDecAcc]
    # # ---------------
    
    
    allCV_meanDecreaseAcc_byOTU <- rowMeans(allCrossValidation_decreaseAccVals)
    allCV_meanDecreaseGini_byOTU <- rowMeans(allCrossValidation_decreaseGiniVals)
    
    # Sort them in decreasing order for plotting
    sortedAllCV_meanDecreaseAcc_byOTU <- allCV_meanDecreaseAcc_byOTU[order(allCV_meanDecreaseAcc_byOTU, decreasing = TRUE)]
    sortedAllCV_meanDecreaseGini_byOTU <- allCV_meanDecreaseGini_byOTU[order(allCV_meanDecreaseGini_byOTU, decreasing = TRUE)]
    
    # Make the plots just like varImpPlot()
    pdf(paste(plotFilePath, fileNameStem, "_allCVItersMeanDecAccuracyPlot.pdf", sep = ""))
    dotchart(rev(sortedAllCV_meanDecreaseAcc_byOTU[1:50]), main = "Mean decrease in accuracy\n across all 10-fold cross validated\n random forests")
    dev.off()
     pdf(paste(plotFilePath, fileNameStem, "_allCVItersMeanDecGiniPlot.pdf", sep = ""))
    dotchart(rev(sortedAllCV_meanDecreaseGini_byOTU[1:50]), main = "Mean decrease in Gini score\n across all 10-fold cross validated\n random forests")
    dev.off()
    
    cvAccuracyStatsForPlot <- data.frame("accuracy" = allCrossValidation_confuseMatrixStats[1,],
                                         "lower_CI" = allCrossValidation_confuseMatrixStats[3,],
                                         "upper_CI" = allCrossValidation_confuseMatrixStats[4,],
                                         "nullAccuracy" = allCrossValidation_confuseMatrixStats[5,],
                                         "pVal" = allCrossValidation_confuseMatrixStats[6,],
                                         "CVFold" = seq(1,numCVFolds,1))
    
    
    pdf(paste(plotFilePath, fileNameStem, "_allCVItersAccuracyComparePlot.pdf", sep = ""), width = 11, height = 8.5)
    cvAccPlot <- ggplot(data = cvAccuracyStatsForPlot, aes(x = CVFold, y = accuracy)) +
        geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0, color = "gray") +
        geom_point(aes(x = CVFold, y = accuracy), size = 3) +
        geom_segment(aes(x = seq(0.5,numCVFolds - 0.5,1), y = nullAccuracy, xend = seq(1.5,numCVFolds + 0.5,1), yend = nullAccuracy), color = "red", linetype = "dashed") +
        ggtitle("Prediction accuracy on test set across all CV folds\n with binomial 95%CI and actual proportion of largest class (remnant sites) for the fold\n in red")
    print(cvAccPlot)
    dev.off()
    
    
    save(allCrossValidationResults, cvAccuracyStatsForPlot, sortedAllCV_meanDecreaseAcc_byOTU, sortedAllCV_meanDecreaseGini_byOTU,
         file = paste(plotFilePath, fileNameStem, "_2000Iter_RFData.r", sep = ""))

    # Aggregate the accuracies for the different folds
    allFoldAccuracyHolder <- sapply(allCrossValidationResults, function(x){return(x$confusionMatrix$overall[1])})
    
    # Calc the average accuracy over all the CV folds
    meanAccAllCVFolds <- mean(allFoldAccuracyHolder)
    
    print(paste("The mean accuracy over all the CV folds was:", meanAccAllCVFolds))
}

# These are the 3 different modules for running the random forests on each of the different data partitions.
# ==========================
# ------------------
# For classifying as disturbed/remnant using all sites (both west and east sites) with 10 fold CV:

binarySubsetOTUframeGradient_onlyAMFOTUs <- subsetOTUframeGradient_onlyAMFOTUs

# This re-codes the OTU table as presence/absence only (1/0)
binarySubsetOTUframeGradient_onlyAMFOTUs[binarySubsetOTUframeGradient_onlyAMFOTUs > 0] <- 1

# Add the distGroup
binarySubsetOTUframeGradient_onlyAMFOTUsWithLabels <- cbind(distGroupDF, binarySubsetOTUframeGradient_onlyAMFOTUs)

# This is important. The rownames of the data.frame used for the Random Forest are how the prediction
# output refers to the different samples. If I don't change the row names here in the subset matrix to be
# 1:end, then they are still the rownames from the original data matrix before I screened out any samples
rownames(binarySubsetOTUframeGradient_onlyAMFOTUsWithLabels) <- seq(1,nrow(binarySubsetOTUframeGradient_onlyAMFOTUsWithLabels),1)

OTUTableToUse <- binarySubsetOTUframeGradient_onlyAMFOTUsWithLabels

# Creates 10 test folds for decision tree accuracy testing that is aware of the disturbance labels while making the folds
# ie the 'disturbed' and 'remnant' samples are roughly evenly distributed among all the folds. This returns a list
# of length 10 with each list element being a vector of indices of the TEST partition (need to use '-' indexing to build
# decision trees on the remaining rows - ie the TRAINING partition). If set the seed immediately before running
# testFolds, it will always return identical folds each time it is run.

numCVFolds <- 10

set.seed(567)
cvFolds <- createFolds(binarySubsetOTUframeGradient_onlyAMFOTUsWithLabels$distGroup, k = numCVFolds, list = TRUE, returnTrain = FALSE)

# Look at the folds with their row numbers and corresponding disturbance factors to verify the distribution of
# sites among folds. Use 'print' to look at it, not 'View'

#cvFoldSampleNamesForVerify <- sapply(cvFolds, function(x){do.call(cbind, list(x,sampleNames[x], as.character(distGroupDF[x,1])))})

# This will be passed to the randomForest calc function to select the variable in the input data frame
# to predict on. This makes the randomForest function modular for testing both E/W and dist/remn.
factorForTest <- "distGroup"

randomForestCall <- paste("as.factor(", factorForTest, ") ~ .", sep = "")
randomForestFormula <- as.formula(randomForestCall)

plotFilePath <- "~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/autoGeneratedScriptOutput/"
fileNameStem <- "GradientCoordinated_allSamples_10CV_RandomForest"

print("*** Running random forest classifier for all samples (10-fold CV) ***")

runRandomForest(factorForTest, OTUTableToUse, cvFolds, plotFilePath, fileNameStem, randomForestFormula)

#-------------------

# ------------------
# For classifying as West/East using only Disturbed Sites with 5 fold CV:

binarySubsetOTUframeGradient_onlyAMFOTUs <- subsetOTUframeGradient_onlyAMFOTUs[distGroup == "Disturbed",]
geogrGroupWE_DF_dist <- geogrGroupWE_DF[distGroup == "Disturbed",]
sampleNames <- sampleNames[distGroup == "Disturbed"]

# This re-codes the OTU table as presence/absence only (1/0)
binarySubsetOTUframeGradient_onlyAMFOTUs[binarySubsetOTUframeGradient_onlyAMFOTUs > 0] <- 1

# Add the geographic group
binarySubsetOTUframeGradient_onlyAMFOTUsWithLabels <- cbind(geogrGroupWE_DF_dist, binarySubsetOTUframeGradient_onlyAMFOTUs)

# This is important. The rownames of the data.frame used for the Random Forest are how the prediction
# output refers to the different samples. If I don't change the row names here in the subset matrix to be
# 1:end, then they are still the rownames from the original data matrix before I screened out any samples
rownames(binarySubsetOTUframeGradient_onlyAMFOTUsWithLabels) <- seq(1,nrow(binarySubsetOTUframeGradient_onlyAMFOTUsWithLabels),1)

OTUTableToUse <- binarySubsetOTUframeGradient_onlyAMFOTUsWithLabels

# Creates 5 test folds for decision tree accuracy testing that is aware of the disturbance labels while making the folds
# ie the 'disturbed' and 'remnant' samples are roughly evenly distributed among all the folds. This returns a list
# of length 5 with each list element being a vector of indices of the TEST partition (need to use '-' indexing to build
# decision trees on the remaining rows - ie the TRAINING partition). If set the seed immediately before running
# testFolds, it will always return identical folds each time it is run.

numCVFolds <- 5

set.seed(567)
cvFolds <- createFolds(binarySubsetOTUframeGradient_onlyAMFOTUsWithLabels$geogrGroupWE, k = numCVFolds, list = TRUE, returnTrain = FALSE)

# This will be passed to the randomForest calc function to select the variable in the input data frame
# to predict on. This makes the randomForest function modular for testing both E/W and dist/remn.
factorForTest <- "geogrGroupWE_DF_dist"

randomForestCall <- paste("as.factor(", factorForTest, ") ~ .", sep = "")
randomForestFormula <- as.formula(randomForestCall)

plotFilePath <- "~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/autoGeneratedScriptOutput/"
fileNameStem <- "GradientCoordinated_DisturbedSamples_5CV_RandomForest"

print("*** Running random forest classifier for only disturbed sites (5-fold CV) ***")

runRandomForest(factorForTest, OTUTableToUse, cvFolds, plotFilePath, fileNameStem, randomForestFormula)
#-------------------


# ------------------
# For classifying as West/East using only Remnant Sites with 5 fold CV:

binarySubsetOTUframeGradient_onlyAMFOTUs <- subsetOTUframeGradient_onlyAMFOTUs[distGroup == "Remnant",]
geogrGroupWE_DF_remn <- geogrGroupWE_DF[distGroup == "Remnant",]
sampleNames <- sampleNames[distGroup == "Remnant"]

# This re-codes the OTU table as presence/absence only (1/0)
binarySubsetOTUframeGradient_onlyAMFOTUs[binarySubsetOTUframeGradient_onlyAMFOTUs > 0] <- 1

# Add the geographic group
binarySubsetOTUframeGradient_onlyAMFOTUsWithLabels <- cbind(geogrGroupWE_DF_remn, binarySubsetOTUframeGradient_onlyAMFOTUs)

# This is important. The rownames of the data.frame used for the Random Forest are how the prediction
# output refers to the different samples. If I don't change the row names here in the subset matrix to be
# 1:end, then they are still the rownames from the original data matrix before I screened out any samples
rownames(binarySubsetOTUframeGradient_onlyAMFOTUsWithLabels) <- seq(1,nrow(binarySubsetOTUframeGradient_onlyAMFOTUsWithLabels),1)

OTUTableToUse <- binarySubsetOTUframeGradient_onlyAMFOTUsWithLabels

# Creates 5 test folds for decision tree accuracy testing that is aware of the disturbance labels while making the folds
# ie the 'disturbed' and 'remnant' samples are roughly evenly distributed among all the folds. This returns a list
# of length 5 with each list element being a vector of indices of the TEST partition (need to use '-' indexing to build
# decision trees on the remaining rows - ie the TRAINING partition). If set the seed immediately before running
# testFolds, it will always return identical folds each time it is run.

numCVFolds <- 5

# This uses the same random number seed as the LDA (which didn't work well with set.seed(567))
set.seed(930)
cvFolds <- createFolds(binarySubsetOTUframeGradient_onlyAMFOTUsWithLabels$geogrGroupWE, k = numCVFolds, list = TRUE, returnTrain = FALSE)

# This will be passed to the randomForest calc function to select the variable in the input data frame
# to predict on. This makes the randomForest function modular for testing both E/W and dist/remn.
factorForTest <- "geogrGroupWE_DF_remn"

randomForestCall <- paste("as.factor(", factorForTest, ") ~ .", sep = "")
randomForestFormula <- as.formula(randomForestCall)

plotFilePath <- "~/Box Sync/R_code/Gradient_SequenceAndEnviron_Analysis/ScriptsForGitHubArchiving/autoGeneratedScriptOutput/"
fileNameStem <- "GradientCoordinated_RemnantSamples_5CV_RandomForest"

print("*** Running random forest classifier for only remnant sites (5-fold CV) ***")

runRandomForest(factorForTest, OTUTableToUse, cvFolds, plotFilePath, fileNameStem, randomForestFormula)
#-------------------