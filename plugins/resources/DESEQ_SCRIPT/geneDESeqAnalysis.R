# Copyright (c) 2011  by Cornell University and the Cornell Research
# Foundation, Inc.  All Rights Reserved.
#
# Permission to use, copy, modify and distribute any part of GobyWeb web
# application for next-generation sequencing data alignment and analysis,
# officially docketed at Cornell as D-5061 ("WORK") and its associated
# copyrights for educational, research and non-profit purposes, without
# fee, and without a written agreement is hereby granted, provided that
# the above copyright notice, this paragraph and the following three
# paragraphs appear in all copies.
#
# Those desiring to incorporate WORK into commercial products or use WORK
# and its associated copyrights for commercial purposes should contact the
# Cornell Center for Technology Enterprise and Commercialization at
# 395 Pine Tree Road, Suite 310, Ithaca, NY 14850;
# email:cctecconnect@cornell.edu; Tel: 607-254-4698;
# FAX: 607-254-5454 for a commercial license.
#
# IN NO EVENT SHALL THE CORNELL RESEARCH FOUNDATION, INC. AND CORNELL
# UNIVERSITY BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
# OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF
# WORK AND ITS ASSOCIATED COPYRIGHTS, EVEN IF THE CORNELL RESEARCH FOUNDATION,
# INC. AND CORNELL UNIVERSITY MAY HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
#
# THE WORK PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE CORNELL RESEARCH
# FOUNDATION, INC. AND CORNELL UNIVERSITY HAVE NO OBLIGATION TO PROVIDE
# MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  THE CORNELL
# RESEARCH FOUNDATION, INC. AND CORNELL UNIVERSITY MAKE NO REPRESENTATIONS AND
# EXTEND NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT
# NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
# PARTICULAR PURPOSE, OR THAT THE USE OF WORK AND ITS ASSOCIATED COPYRIGHTS
# WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.

# This script runs through a DESeq Analysis
# Nyasha Chambwe June 2, 2010 & Kevin C. Dorff July 19, 2010
# July 24, 2012: DESeq function removed: changed dispersion estimate function name from "estimateVarianceFunctions" to "estimateDispersions"

library("DESeq")
library("Cairo")

#--------------------------------------------------------------------------
# This version of the script is designed for running manually.
# Supply the command line options
#   input       - tsv input file
#   output      - tsv output file
#   graphOutput - png output file
#   elementType - the kind of element to keep for analysis, must match the second column of this file.
#  Example command line:
#
#  R -f geneDESeqAnalysis.R --slave --quiet --no-restore --no-save --no-readline --args exonInput=input-exon.tsv geneInput=input-gene.tsv output=alloutput.tsv graphOutput=output.png
#
# For discussion about the input and output files, see the documentation
# above the function "processFile".
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# FUNCTIONS FOR DIFF EXP PROCESSING AND SUPPORT
#--------------------------------------------------------------------------

#
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
# Given conditionsToCheck and groupNameToCheck, return groupSizes, the
# size of each group in the order of groupNameToCheck.
#
groupSizeFunction <- function(conditionsToCheck, groupNameToCheck) { 
  sum <- 0
  for (conditionToCheck in conditionsToCheck) {
    if (conditionToCheck == groupNameToCheck) {
      sum <- sum + 1
    }
  }
  sum
}

#
# This will place the suffix (file + "" + suffix) at the end
# of the graphOutputFile but before the file extension. If suffix is "GENE"
# and graphOutputFile is "this.png" this will output
# "thisGENE.png". If the file has no extension, such as "thispng"
# this will output "thispngGENE".
#
generateGraphFilename <- function(graphOutputFile, suffix) {
  parts <- array(unlist(strsplit(graphOutputFile, ".", fixed=TRUE)))
  extension <- tolower(parts[length(parts)])
  if (graphOutputFile == extension) {
    newFilename <- paste(graphOutputFile, "", suffix, sep="")
    newFilename
  } else {
    newFilename <- paste(substring(graphOutputFile, 1, nchar(graphOutputFile)-nchar(extension)-1),"", suffix , ".", extension, sep="")
    newFilename
  }
}

#
# Performs the DiffExp calculation using DESeq on inputFile and write
# the results to outputFile.
#
# * Process an input TSV file to an output file + optional graph file.
# * The input should contain data of EITHER "EXON" or "GENE" (specified in
#   elementType).
# * If appendOutputFile is FALSE, this will create a new file with column
#   headers.
# * If appendOutputFile is TRUE, this will append to an existing file and
#   not repeat the column headers.
# * The first two columns of the output file will be "element-id" and
#   "element-type".
# * The input file should be created by Goby's "alignment-to-annotation-counts"
#   mode with the option "--eval counts". That input file should have been
#   created for either "GENE" or "EXON" by using EITHER
#   "--include-annotation-types gene" or "--include-annotation-types exon"
#   (but don't use both). If you need both, create TWO input files and run
#   the second with appendOutputFile=TRUE.
# * One of the key features of running with "--eval counts" is that the
#   column names of the input file CONTAIN the group name for the column such
#   as a column named "count sample seqc-ilm-s1 [ILM]" we know are counts
#   for the sample "seqc-ilm-s1" and for the group "ILM".
# * If graphOutputFile, which should be a ".png" filename, a png graph will
#   be created at 700x700 pixels. The ACTUAL filename will include elementType,
#   so if elementType is "GENE" and graphOutputFile is "something.png" the
#   actual file that will be written will be "something-GENE.png"
#
processFile <- function(inputFile, sampleGroupMapping, outputFile, graphOutputFile, elementType, appendOutputFile) {
  # Read the tab delimited input file
  countsTable <- read.table(inputFile, header=TRUE, stringsAsFactors=TRUE, sep="\t", check.names=FALSE)
  head(countsTable)
  
  # replace row names with gene identifiers
  rownames(countsTable) <- countsTable$"element-id"
  
  countsTable <- countsTable[, -1]  # Slice out the element-id column
  countsTable <- countsTable[, -1]  # Slice out the element-type column
  colNames <- colnames(countsTable)
  
  # this vector describes how samples are assigned to a particular group
  conditions <- groupNameFunction(sampleGroupMapping, colNames)

  # The output marks columns "A" and "B". These are what "A" and "B" actually are
  groupNames <- unique(conditions)
  print(groupNames)
  groupSizes <- c()
  for (groupName in groupNames) {
    groupSizes <- append(groupSizes, groupSizeFunction(conditions, groupName))
  }
  
  # instantiate a new CountDataSet object
  countDataSet <- newCountDataSet(countsTable, conditions)
  
  # Alternatively could use the actual total numbers of reads or estimate them from data as below
  countDataSet <- estimateSizeFactors(countDataSet)
  sizeFactors(countDataSet)
  head(countDataSet)
  
  # Because this experiment has only one replicate per condition, we have to pool both samples in order to estimate variance
  
  if (sum(groupSizes) == length(groupNames)) {
    cat("Running with no replicates found\n")
    countDataSet <- estimateDispersions(countDataSet, pool=TRUE)
  } else {
    cat("Running with replicates found in at least one group\n")
    countDataSet <- estimateDispersions(countDataSet)
  }
  head(countDataSet)
  
  # This function tests for differences between the base means of two conditions (i.e., for differential expression in the case of RNA-Seq).
  result <- nbinomTest(countDataSet, groupNames[1], groupNames[2])
  head(result)
  
  #
  # Rename some of the columns
  #
  resultColNames <- colnames(result)
  newColNames <- c()
  for (resultColName in resultColNames) {
    if (resultColName == "id") {
      resultColName = "element-id"
    }
    baseGroupName <- substring(resultColName, 1, nchar(resultColName)-1)
    if (baseGroupName == "baseMean" || baseGroupName == "resVar") {
      groupCharacter <- substring(resultColName, nchar(resultColName), nchar(resultColName))
      index <- 1
      if (groupCharacter == "A") {
        index <- 1
      } else if (groupCharacter == "B") {
        index <- 2
      }
      newColNames <- append(newColNames, paste(baseGroupName, "-", groupNames[index], sep=""))
    } else {
      newColNames <- append(newColNames, resultColName)
    }
  }
  colnames(result) <- newColNames
  head(result)
  
  # Insert the element-type column
  numRows <- length(rownames(result))
  numCols <- length(colnames(result))
  result <- data.frame("element-id"=result[1:numRows,1],"element-type"=c(elementType),result[1:numRows,2:numCols],check.names=FALSE)
  head(result)

  write.table(result, file=outputFile, append=appendOutputFile, col.names=!appendOutputFile, sep="\t", quote=FALSE, row.names=FALSE)
  
  if (graphOutputFile != "") {
    graphOutputFile <- generateGraphFilename(graphOutputFile, elementType)
    # MA-Plot for this analysis
    CairoPNG(graphOutputFile, width=700, height=700)
    pValueThreshold <- 0.1
    plot(result$baseMean, result$log2FoldChange, xlab="Mean of Counts", ylab="Log2 Fold Change", log="x", pch=19, cex= .25, col=ifelse(result$padj < pValueThreshold, "red", "black"))
    title(main = list(paste("MA-Plot ", elementType, groupNames[1], "vs", groupNames[2]), cex=1.5, col="black", font=3))
    dev.off()
  }
}

#--------------------------------------------------------------------------
# COMMAND LINE PARSING
#--------------------------------------------------------------------------

input <- ""
output <- ""
graphOutput <- ""
sampleGroupMapping <- ""
elementType <- ""

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

#--------------------------------------------------------------------------
# Process given the specifed command line
#--------------------------------------------------------------------------

appendOutputFile <- FALSE
processFile(input, sampleGroupMapping, output, graphOutput, elementType, appendOutputFile)
appendOutputFile <- TRUE

