#=====================================================================================
#  Code chunk 1
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "E:/R_WGCNA_test/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"; #need to be changed
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
GTEx_Data <- read.delim("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", row.names=1);
dim(GTEx_Data); #show the dimension of this data
names(GTEx_Data); # show the columns names 

#=====================================================================================
#  Code chunk 2
#=====================================================================================

GTEx_Data1 <- GTEx_Data[!duplicated(GTEx_Data$Description),];# remove rows with the same symbol names
rownames(GTEx_Data1) <- GTEx_Data1$Description;
row.names(GTEx_Data1)
datExpr0 = as.data.frame(t(GTEx_Data1[-c(1)]));

#=====================================================================================
#  Code chunk 3
#=====================================================================================

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK


#=====================================================================================
#  Code chunk 4
#=====================================================================================


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
datExpr0 = log2(datExpr0+0.00001)

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================

# Plot a line to show the cut
abline(h = 305000, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 350000, minSize = 1)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================

# read traits data
datTraits = read.csv("traits.csv",row.names = 1);



#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.

jpeg(filename = "Sample dendrogram and trait heatmap.jpg", pointsize=14,
     width=1200,height = 1000, quality = 200);
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap");
dev.off();


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


save(datExpr, datTraits, file = "tissueSpecific-01-dataInput.RData")

