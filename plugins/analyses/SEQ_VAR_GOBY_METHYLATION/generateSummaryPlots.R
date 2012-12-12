# This script takes methyl stats output files and renders summary QC plots for display on GobyWeb Results page
# Author: Nyasha Chambwe 
# Date: June 27, 2012
#--------------------------------------------------------------------------
require("ggplot2")
require("Cairo")
#--------------------------------------------------------------------------
# This version of the script is designed for running manually.
# Supply the command line options
#   depthsFile   - tsv input file for genes, optional
#   conversionFile   - tsv input file for genes, optional
#   coverageGraphOutput  - png output file
#   conversionGraphOutput - png output file
#
#  Example command line:
#  R -f generateSummaryPlots.R --slave --quiet --no-restore --no-save --args depthsFile=depths.tsv conversionFile=conversion-rates.tsv coverageGraphOutput=coverage.png conversionGraphOutput=conversion.png 

# 
#--------------------------------------------------------------------------
# Function to calculate the sample average density per depth bin over the positive and negative strand 
average.densities.per.sample <- function(depthsFile){
  # copy half the depths File
  sampleDepthsFile <- depthsFile[ 1:(dim(depthsFile)[1]/2), ]
  # remove strand  and num-observed-in-bin columns
  sampleDepthsFile <- sampleDepthsFile[,-c(1,7)]
  index.positive.strand <- 1
  index.negative.strand <- (dim(depthsFile)[1]/2) + 1
  # average density over positive and negative strand
  for(i in index.negative.strand:(dim(depthsFile)[1])){
    pos.density <- depthsFile[index.positive.strand,6]
    neg.density <- depthsFile[index.negative.strand,6]
    avg.density<- (pos.density+ neg.density)/2
    sampleDepthsFile[index.positive.strand,5] <- avg.density
    index.positive.strand= index.positive.strand + 1
  }
  return(sampleDepthsFile)
}
#--------------------------------------------------------------------------
# Function to standardize the appearance of a ggplot by formatting axis appearance, adding plot titles, removing grid lines and default grey panels
formatPlotAppearance<- function(plotObject, title){
  # Format x and y axis appearance
  plotObject <- plotObject + opts(axis.title.x= theme_text(face= "bold", hjust= 0.5,  lineheight = 20, colour="black")) + opts(axis.title.y= theme_text(face= "bold", lineheight = 20, colour="black", angle=90))
  # Add plot title, format font and size
  plotObject <- plotObject + opts(title=title)
  plotObject <- plotObject + opts(plot.title = theme_text(size=2, face="bold",  hjust= 0.5,  lineheight =5,  colour="black"))
  # Remove grid lines
  plotObject <- plotObject + opts(panel.grid.major=theme_blank(), panel.grid.minor=theme_blank())
  plotObject <- plotObject + theme_bw()
  return(plotObject)
}
#--------------------------------------------------------------------------
# Function to draw sample coverage plots
# input: dataTable - *depths.tsv file 
# input: title - character string for plot title
drawSamplesDepthPlots <- function(dataTable, title){
  depthsPlot <- ggplot(data=dataTable, aes(x = depthMidPoint, y=density))
  p <- depthsPlot + geom_point() + geom_smooth() 
  p <- p + scale_x_log10()
  # Change x and y axis labels
  p <- p + xlab("Coverage") + ylab("Density")
  p <- formatPlotAppearance(p, title) 
  p <- p + facet_wrap(~ sample, scales="free_x", ncol=2)
  p <- p + opts(strip.text.x= theme_text(size=9))
               # strip.background = theme_rect(colour="white", fill="lightgrey"))
  return(p)
}
#--------------------------------------------------------------------------
# Function to calculate the percent non-CpG conversion for each sample
# input: *conversion-rates.tsv
# headers: strand sample num-converted-not-in-CpG-context num-not-in-CpG-context percent-converted-in-non-CpG-context
average.non.cpg.conversion <- function(conversionRatesFile){
  # copy half the conversion rates file
  sampleConversionRatesFile <- conversionRatesFile[ 1:(dim(conversionRatesFile)[1]/2), ]
  index.positive.strand <- 1
  index.negative.strand <- (dim(conversionRatesFile)[1]/2) + 1
  # average conversion counts over positive and negative strand
  for(i in index.negative.strand:(dim(conversionRatesFile)[1])){
    pos.num.converted <- conversionRatesFile[index.positive.strand,3]
    neg.num.converted <- conversionRatesFile[index.negative.strand,3]
    pos.num.noncpg <- conversionRatesFile[index.positive.strand,4]
    neg.num.noncpg <- conversionRatesFile[index.negative.strand,4]
    sum.num.converted <- pos.num.converted + neg.num.converted
    sum.num.noncpg <- pos.num.noncpg + neg.num.noncpg
    perc.conversion<- (sum.num.converted/sum.num.noncpg)*100
    sampleConversionRatesFile[index.positive.strand,5] <- perc.conversion
    index.positive.strand= index.positive.strand + 1
  }
  # remove strand column
  sampleConversionRatesFile <- sampleConversionRatesFile[,-1]
  return(sampleConversionRatesFile)
}
#--------------------------------------------------------------------------
# Function to draw histogram for the Non-CpG conversion rate across all samples
# input: *conversion-rates.tsv
drawConversionPlot <- function(dataTable, title){
  conversionPlot <- ggplot(data=dataTable, aes(x=percentConvertedInNonCpGcontext))
  q <- conversionPlot + geom_histogram(binwidth=5) + 
  scale_x_continuous(limits=c(0,100), breaks=seq(0,100,10)) + coord_flip()
  # Change x and y axis labels
  q <- q + ylab("Count") + xlab("Conversion in Non-CpG Context (%)")
  q <- formatPlotAppearance(q, title)
  return(q)
}

#--------------------------------------------------------------------------
# COMMAND LINE PARSING
#--------------------------------------------------------------------------

depthsFile   <- ""
conversionFile  <- ""
coverageGraphOutput <- ""
conversionGraphOutput <-""

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

if ((depthsFile == "") || (conversionFile == "") || coverageGraphOutput == "" || conversionGraphOutput=="") {
  stop("This script requires two graph output files to be defined and tab delimited files of sample depths and conversion rates to be specified")
}

#--------------------------------------------------------------------------

depths <- read.table(depthsFile, header=TRUE, stringsAsFactors=TRUE, sep="\t", check.names=FALSE)
num.depth.bins <- dim(table(depths[,3]))
num.samples <- (dim(depths)[1]/num.depth.bins)/2
sampleDepths <- average.densities.per.sample(depths)
colnames(sampleDepths) <- c("sample","depth-bin", "depthMidPoint", "log2OfDepthMidpoint", "density")
plot.width=700
num.rows=ceiling(num.samples/2)
plot.height=num.rows*233
  
CairoPNG("coverage.png", width=plot.width, height=plot.height)
drawSamplesDepthPlots(sampleDepths, "Read Coverage Across Samples")  
dev.off()

#--------------------------------------------------------------------------
conversionRates <- read.table(conversionFile, header=TRUE, stringsAsFactors=TRUE, sep="\t", check.names=FALSE)
sampleConversionRates <- average.non.cpg.conversion(conversionRatesFile=conversionRates)
colnames(sampleConversionRates) <- c("sample","numConvertedNotInCpGcontext", "numNotInCpGcontext", "percentConvertedInNonCpGcontext")
CairoPNG("conversion.png", width=350, height=350)
drawConversionPlot(sampleConversionRates, "Bisulfite Conversion Rates")
dev.off()
