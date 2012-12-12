# R script to carry out  differential expression analysis using methods implemented in the edgeR Bioconductor Package
# Author       : Nyasha Chambwe
# Date Created : November 29, 2011 
# Date Edited  : July 20, 2012 --- modified script to integrate with GobyWeb -- developed based on edgeR  version 2.6.9
#--------------------------------------------------------------------------------------------------------------------------------------
#
# This version of the script is designed for running manually.
# Supply the command line options
#   input   - tsv input file for genes, optional
#   output      - tsv output file
#   graphOutput - png output file
#
#  Example command line:
#
#  R -f deAnalysisEdgeR.R --slave --quiet --no-restore --no-save --no-readline 
# --args input=input.tsv output=alloutput.tsv mdsPlotOutput=mds.png smearPlotOutput=smear.png sampleGroupMapping=sampleGroupMapping.tsv normalizationMethod="RLE" dispersionMethod="common"

# For discussion about the input and output files, see the documentation
# above the function "processFile".
#
#--------------------------------------------------------------------------------------------------------------------------------------
#
# Ascertain that all required R packages and their dependencies are loaded
require("edgeR", quietly=FALSE)
require("limma", quietly=FALSE)
require("Cairo", quietly=FALSE)
#
#--------------------------------------------------------------------------------------------------------------------------------------
# Functions for DIFF EXP processing and support
#-------------------------------------------------------------------------------------------------------------------------------------
# Given a sample to group mapping file and a vector of sampleNames return a vector of group membership for those samples in the order given by the sampleNames vector
#
groupNameFunction <- function(sampleGroupMapping, colNames){ 
  # Read the tab delimited sample to group mapping file
  sampleToGroupTable <- read.table(sampleGroupMapping, header=FALSE, stringsAsFactors=TRUE, sep="\t", check.names=FALSE)
  groupVector <- vector() 
  # Assign group labels to each sample based on the sample to group mapping file, in the order of the counts file seen in variable colNames
  for(col in colNames){
    sampleIndex <-which(sampleToGroupTable$V1==col)
    grp <- as.character(sampleToGroupTable[sampleIndex,]$V2)
    groupVector <- append(groupVector, grp)
  }
  return(groupVector)
}
#
#-------------------------------------------------------------------------------------------------------------------------------------
# * inputFile is created by Goby's "alignment-to-annotation-counts" mode
processInputFile <- function(inputFile, sampleGroupMapping){
  # Read the tab delimited counts input file
  countsTable <- read.table(inputFile, header=TRUE, stringsAsFactors=TRUE, sep="\t", check.names=FALSE) 
  # Replace row names with gene identifiers
  rownames(countsTable) <- countsTable$"element-id" 
  countsTable <- countsTable[,-1]  # slice out the element-id column
  countsTable <- countsTable[,-1]  # slice out the element-type column
  countsTable <- as.matrix(countsTable)
  sampleNames <- colnames(countsTable)
  
  # define a "group" vector that assigns each sample/library to a group/condition for comparison purposes
  targets <- groupNameFunction(sampleGroupMapping, sampleNames)
  
  # Setting up DGEList object
  d = DGEList(counts = countsTable, lib.size = colSums(countsTable), group = targets, remove.zeros = TRUE)
  return(d)
}
#
#-------------------------------------------------------------------------------------------------------------------------------------
# Function that performs normalization on a DGEList object according to a specific method 
performNormalization <- function(dgeObj, normalizationMethod, filterFlag){
  if(filterFlag=="true"){
    ########################################################################################################
    # Filtering low count genes (optional) some say necessary for speed, others say it biases results
    # Since it is not possible to achieve statistical significance with few total counts
    #in total for a tag, we filter out tags with n or fewer counts in total
    print("Filtering out low count tags.............................................................................")
    # Find the number of samples in the smallest group under consideration
    x <- min(table(dgeObj$samples$group))
    cpm.dgeObj <- cpm(dgeObj)
    filteredDgeObj <- dgeObj[rowSums(cpm.dgeObj>1) >= x, ]
    normDgeObj <- calcNormFactors(filteredDgeObj, method=normalizationMethod)
  } else{
    # Perform normalization without filtering
    print("Normlizing data WITHOUT filtering.............................................................................")
    normDgeObj <- calcNormFactors(dgeObj, method=normalizationMethod)
  }
  return(normDgeObj)
  
}
#
#-------------------------------------------------------------------------------------------------------------------------------------
# Function to generate a multi-dimensional scaling plot to show sample relationships
generateMDSPlot <-function(dgeObj, mdsPlotOutput){
  if (mdsPlotOutput != "") {
    # MDS-Plot for this analysis
    CairoPNG(mdsPlotOutput, width=700, height=700)
    mdsPlot <- plotMDS(dgeObj, top=500, labels=dgeObj$samples$group, cex=1, col="blue") # labels - character vector of group labels
    title(main = "Diagnostic Plot\n Sample Relationship based on Multi-Dimensional Scaling")
    dev.off()
  }
}
#
#-------------------------------------------------------------------------------------------------------------------------------------
#
generateSmearPlot <- function (dgeObj, deExact, smearPlotOutput) {
  # Determine the number of differentially expressed genes based on input parameters
  numDeTagsTable <- summary(decideTestsDGE(deExact, adjust.method="BH", p.value=0.05))
  numDeTags <- sum(numDeTagsTable[1,1], numDeTagsTable[3,1])
  deTags <- rownames(topTags(deExact, n = numDeTags)$table)
  
  if (smearPlotOutput != "") {
    # Generate a smear plot for the result with differentially expressed tags highlighed
    CairoPNG(smearPlotOutput, width=700, height=700)
    plotSmear(dgeObj, de.tags = deTags, main = "Smear Plot")
    abline(h = c(-2, 2), col = "dodgerblue", lwd = 2)
    dev.off()
  }
  
}  
#
#-------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------
# 
# Estimates the dispersion parameter for the negative binomial model
# Dispersion method: either "common" or "tagwise"
# If tagwise dispersion - common dispersion estimated first
estimateDispersion <- function(dgeObj, dispersionMethod){ 
  # Estimate Common Negative Binomial Dispersion by Conditional Maximum Likelihood
  d <- estimateCommonDisp(dgeObj)
  if(dispersionMethod=="tagwise"){  
    # Estimates tagwise dispersion values by an empirical Bayes method based on weighted conditional maximum likelihood.
    d <- estimateTagwiseDisp(d)
  }
  return(d)
}
#
#-------------------------------------------------------------------------------------------------------------------------------------
#
# Carry out statistical test for differential expression between two groups
estimateDifferentialExpression <- function(dgeObj, outputFile, smearPlotOutput, sampleGroupMapping){
  # Sanity check
  if(is.null(dgeObj$common.dispersion)){
    stop("The DGEList object provided must at least have a common dispersion estimate in order to carry out the exact test for differentialexpression. \n")}
  
  # Test for significant differential expression between two groups
  deExact <- exactTest(dgeObj)
  resultTable <- topTags(deExact, n=dim(dgeObj$counts)[1], adjust.method="BH", sort.by="p.value")
  
  # gets sample counts per million for each feature as well as the group means to append to result table
  sampleCountsTable <- cpm(dgeObj) 
  sampleNames<- colnames(sampleCountsTable)
  targets <- groupNameFunction(sampleGroupMapping, sampleNames)
  
  print("Getting group mean normalized counts per million ......................................................")
  transposedSampleCounts <- as.data.frame(t(sampleCountsTable) )
  transposedGroupTagMeans <- aggregate(transposedSampleCounts, list(GroupID=targets), mean)
  groupTagMeans <-  t(transposedGroupTagMeans)
  colnames(groupTagMeans) <- groupTagMeans[1,]
  groupTagMeans <-  groupTagMeans[-1,]
  
  countsMeanTable <-  merge(sampleCountsTable, groupTagMeans, by=0)
  row.names(countsMeanTable) <- countsMeanTable[,1]
  countsMeanTable <- countsMeanTable[,-1]
  
  finalResultTable <- merge(countsMeanTable, resultTable, by=0)
  names(finalResultTable)[1] <- "Feature ID"
  
  print("Writing stats results table ......................................................")
  write.table(finalResultTable, file=outputFile, sep="\t", quote=FALSE, row.names=FALSE)
  
  print("Generating Smear Plot ......................................................")
  generateSmearPlot(dgeObj, deExact, smearPlotOutput)
}
#-------------------------------------------------------------------------------------------------------------------------------------
#


runAnalysis <- function(inputFile, sampleGroupMapping, output, mdsPlotOutput,smearPlotOutput, annotationtype, normalizationMethod, dispersionMethod, filterFlag){
  print("processing counts files for edgeR analysis......................................................")
  result <- processInputFile(inputFile, sampleGroupMapping)
  print("carrying out count normalization......................................................")
  resultNormalized <- performNormalization(result, normalizationMethod, filterFlag)
  print("Generating diagnostic MDS plot ......................................................")
  generateMDSPlot(resultNormalized, mdsPlotOutput)
  print("Calculating estimates of dispersion ......................................................")
  resultNormalizedDispersed <- estimateDispersion(resultNormalized, dispersionMethod)
  print("Calculating differential expression statistics ......................................................")
  estimateDifferentialExpression(resultNormalizedDispersed, output, smearPlotOutput, sampleGroupMapping)
}

#--------------------------------------------------------------------------
# COMMAND LINE PARSING
#--------------------------------------------------------------------------

input <- ""
output <- ""
mdsPlotOutput <- ""
smearPlotOutput <- ""
sampleGroupMapping <- ""
elementType <- ""  
normalizationMethod <- ""
dispersionMethod <- ""
filterFlag <- ""

notused <- capture.output(commandArgs())
for (e in commandArgs()) {
  ta = strsplit(e,"=",fixed=TRUE)
  if(! is.na(ta[[1]][2])) {
    temp = ta[[1]][2]
    assign(ta[[1]][1],temp)
  } else {
    assign(ta[[1]][1],TRUE)
  }
}

if ((output == "") || (input == "")) {
  stop("This script requires input and output to be defined")
}
if(sampleGroupMapping==""){
  stop("This script requires a sample to group mapping file to be specified.")
}

#---------------------------------------------------------------------------------------------------
# Process given the specifed command line
#--------------------------------------------------------------------------

runAnalysis(input, sampleGroupMapping, output, mdsPlotOutput, smearPlotOutput, elementType, normalizationMethod, dispersionMethod, filterFlag)
