
# load library and set environment
library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
library(magrittr)      # provides the %>% operator
library(WGCNA)
library(pheatmap)
setwd("F:/data science/R/WGCNA-PhosphProteomeData")
getwd()

# data preprocessing and explore.
df_all <- read.csv("data/phosphoproteomicdata.csv", header = TRUE) #read initial data
df_trim <- df_all[,c(0:16,ncol(df_all))] #data trimming to make data frame with samples are rows and phosite are columns

#duplicated phosite from different protein isoforms need to be distinguished by adding a prefix.
df_trim$ProteinAndSite <- make.unique(df_trim$ProteinAndSite) # make all the values unique adding a suffix
rownames(df_trim) <- df_trim$ProteinAndSite
df_trim$ProteinAndSite <- NULL

#check and confirm data type
sapply(df_trim, class) #check the types of values
class(df_trim) #check the types of dataframe.

# Matrix and numeric types are needed for the following steps
df_trim <- as.matrix(df_trim) # change to matrix
storage.mode(df_trim) <- "numeric" #change values to numeric

#this df is the final data before starting the WGCNA analysis
df <- t(df_trim) #transpose data frame to make rows are samples and columns are phosites
rownames(df) <- gsub("\\.", " ", rownames(df)) #replace . with space to keep alligned
dim(df)

#read meta data containing group information and other measurement
df_meta <- read.csv("data/isoprostanandmetadata.csv", row.names=1, check.names = FALSE)

#QC steps
gsg <- goodSamplesGenes(df) #check if there are missing values or zero variant values across genes and samples.
sampleTree <- hclust(dist(df), method = "average") #Clustering samples based on distance 

#Visualizing sampletree
#Setting the graphical parameters
par(cex = 1);
par(mar = c(1,5,3,1))

#Plotting the cluster dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", cex.lab = 1.2,
     cex.axis = 1.2, cex.main = 1.5)

#An alternative way to check for outliers are pca
pca_df <- prcomp(df, scale. = TRUE)
pca_table <- data.frame(pca_df$x)

par(cex = 1);
par(mar = c(5,5,3,3))
plot(pca_table$PC1, pca_table$PC2,
     xlab = "PC1",
     ylab = "PC2",
     main = "PCA of samples",
     pch = 19)
text(pca_table$PC1, pca_table$PC2, labels = rownames(pca_table), pos =3, cex =0.5)

# you can also color samples by groups
pca_table$group <- df_meta$groups
ggplot(pca_table, aes(x=PC1, y=PC2, color=group))+
  geom_point(size=3)+
  geom_text(aes(label = pca_table$group), vjust = -1, size=3)+
  theme_minimal()+
  labs(title = "PCA of samples colored by group")

#since there are too many genes and some of them are not variable across samples, you can pick up the most variable one for further analysis
phosite_var <- apply(df, 2, var)
top_var_phosite <- names(sort(phosite_var, decreasing = TRUE)[1:2000])
df_topvar <- df[,top_var_phosite]

#histogram plot of variance
hist(phosite_var,
     breaks = 1000,
     main = "Distribution of gene variance",
     xlab = "Variance",
     col = "skyblue")

#density plot of variance
plot(density(phosite_var),
     main = "Density of gene variance",
     xlab = "Variance",
     col = "darkblue",
     lwd = 2)

#variance vs mean plot
phosite_mean <- colMeans(df)
plot(phosite_mean, phosite_var,
     xlab = "Mean",
     ylab = "Variance",
     main = "Variance vs. Mean Expression",
     pch = 16, col = rgb(0, 0, 1, 0.3)) # you can see lower level of express has high variance generally

#variance vs mean plot top2000
phosite_topvar_var <- apply(df_topvar, 2, var)
phosite_topvar_mean <- colMeans(df_topvar)
plot(phosite_topvar_mean, phosite_topvar_var,
     xlab = "Mean",
     ylab = "Variance",
     main = "Variance vs. Mean Expression",
     pch = 16, col = rgb(0, 0, 1, 0.3))


# as soon as the data passed QC, the following steps are formal WGCNA analysis
#pick up threshold to fit a scale free topology network.
powers = c(c(1:10), seq(from = 12, to = 30, by = 2))
sft = pickSoftThreshold(
  df_topvar,             # <= Input data
  powerVector = powers,
  networkType = "unsigned",
  verbose = 5
)

sft.data <- sft$fitIndices

par(mfrow = c(1,2));
plot(sft.data$Power,
     sft.data$SFT.R.sq,
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit R^2",
     main = paste("Scale independence"),
     ylim = c(0,1)
)
text(sft.data$Power,
     sft.data$SFT.R.sq,
     labels = powers, cex = 0.9, col = "red"
)
abline(h = 0.80, col = "red")

plot(sft.data$Power,
     sft.data$mean.k.,
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     main = "Mean connectivity"
)
text(sft.data$Power,
     sft.data$mean.k.,
     labels = powers,
     cex = cex1, col = "red")

# pick 18 as power and run correlation and adjacency
picked_power = 18
#temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(df_topvar,
                          TOMType = "unsigned",
                          power = picked_power,                # <= power here
                          MergeCutHeight = 0.25,
                          minModuleSize = 20,
                          numericLabels = F,
                          randomSeed = 1234,
                          verbose = 3)
#cor <- temp_cor
# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

# this is another way to plot the dendrogrma and modules
diss=1-TOMsimilarityFromExpr(df_topvar, power = 18)
colnames(diss) =rownames(diss) =colnames(df_topvar)

hier=hclust(as.dist(diss), method="average" )
plotDendroAndColors(hier, mergedColors[netwk$blockGenes[[1]]], "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

#plot similarity matrix heatmap
diag(diss) = NA
sizeGrWindow(7,7)
TOMplot(diss, hier, as.character(mergedColors[netwk$blockGenes[[1]]]))

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(df_topvar, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Module to trait association
trait_data <- df_meta[,-ncol(df_meta)]
module.trait.correlation = cor(MEs0, trait_data, use = "p")
module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nrow(df_topvar)) 

# Will display correlations and their p-values
textMatrix = paste(signif(module.trait.correlation, 2), "\n(",
                   signif(module.trait.Pvalue, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.correlation)

# Create a logical mask for significant + strong correlations
highlight_mask = (module.trait.Pvalue < 0.01) & (abs(module.trait.correlation) > 0.5)

# Blank out all entries that do not meet the criteria
textMatrix[!highlight_mask] = ""

# Optional: Add asterisk to highlight significant entries
textMatrix[highlight_mask] = paste0("*", textMatrix[highlight_mask])

# Plot the heatmap with only significant strong correlations labeled
pdf("module_trait_heatmap.pdf", width = 10, height = 8)
par(mar = c(6, 8.5, 3, 1))
labeledHeatmap(
  Matrix = module.trait.correlation,
  xLabels = names(trait_data),
  yLabels = names(MEs0),
  ySymbols = names(MEs0),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 1,
  zlim = c(-1,1),
  main = paste("Module-trait relationships")
)
dev.off()

# pick out a few modules of interest here
modules_of_interest = c("brown", "turquoise", "yellow")

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)
# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)
subexpr <- df_topvar[, colnames(df_topvar) %in% submod$gene_id]

rownames(submod) <- submod$gene_id 

#make a long dataframe for plot
subdf_long = as.data.frame(t(subexpr)) %>%
  mutate(gene_id = row.names(.)) %>%
  pivot_longer(-gene_id) %>%
  mutate(module = submod[gene_id,]$colors)

subdf_long %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")

#Turquoise module is associated with most isoprostan traits, what are turquoise proteins
turquoise_gene_id <- subset(module_df,colors %in% "turquoise")$gene_id
turquoise_expr <- df_topvar[, colnames(df_topvar) %in% turquoise_gene_id]

#calculate correlation between phosite in turquoise module and isoprostan trait
gene.trait.correlation = cor(turquoise_expr, trait_data, use = "p")
gene.trait.Pvalue = corPvalueStudent(gene.trait.correlation, nrow(df)) 

#the module trait correlation show that some genes in turquoise module are strongly related to isoprostan levels
jpeg("turquoisemodule gene and triat cor.jpeg", width = 1000, height = 1500, res=200)
#par(mar = c(6, 8.5, 3, 1))
pheatmap(gene.trait.correlation,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         fontsize_row = 6,
         cellwidth = 15,
         cellheight = 5,
         clustering_method = "complete",  # or "ward.D", "average", etc.
         display_numbers = FALSE)          # Optional: show correlation values
dev.off()



