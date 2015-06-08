# ----------------------------------------------------------------------
# FLS R group 2015-06-09: Basic heatmaps of micriobiome data
# ----------------------------------------------------------------------



# ----------------------------------------------------------------------
# INTRO
# ----------------------------------------------------------------------
# Such choice, many defaults (wow)
#   - heatmap()
#   - heatmap.2()
#   - aheatmap()
#   - heatplot()
#   - ggplot2's geom_tile()
# These are mostly pretty similar, but defaults vary a lot

# OTUs for fecal 16S sequences from mothers and their monozygotic twins. Reyes et al. (2013)
# Were the abundance values not already proportional, we could do:
# reyes_prop <- reyes/rowSums(reyes)
# Required packages: vegan, RColorBrewer, gplots, ggplot2, reshape, stringr



reyes <- read.csv('16s_taxa_all.rdp_9.scaled.tsv', sep='\t')
# reyes <- read.csv('virome_taxa_all.lambda.scaled.tsv', sep='\t')
rownames(reyes) <- reyes$TAXON
reyes <- reyes[, -1]
reyes <- as.matrix(reyes)



# ----------------------------------------------------------------------
# HEATMAP()
# ----------------------------------------------------------------------



# Minimum viable heatmap with woeful defaults
heatmap(reyes)
heatmap(reyes, Colv=NA, Rowv=NA)

# Some colour
library(RColorBrewer)
jBuPuFun <- colorRampPalette(brewer.pal(n = 9, "BuPu"))
paletteSize <- 256
jBuPuPalette <- jBuPuFun(paletteSize)

# Discard taxa to improve readability?
# ...High pass filter? (discard rows with low row-wise maxima)
taxon_max_abundances <- apply(reyes, 1, max)
reyes_hpf <- reyes[which(taxon_max_abundances > 0.001),]
heatmap(reyes_hpf, Colv=NA, Rowv=NA, col=jBuPuPalette)
# low_abundance_taxa = names(which(max_taxon_abundances < 0.01))
# reyes_cull <- reyes[, -which(names(reyes) %in% low_abundance_taxa)]
# ...Or alternatively discard rows with low variance (probably uninteresting)
taxon_variances <- apply(reyes, 1, var)
reyes_most_variable <- reyes[taxon_variances > quantile(taxon_variances, 0.8),]
heatmap(reyes_most_variable, Colv=NA, Rowv=NA, col=jBuPuPalette)

# Now, but using an appropriate distance measure for clustering taxa
# Vegan library Provides Bray-Curtis dissimilarity measure
library(vegan)
bray_dist_mat <- vegdist(reyes_most_variable, method='bray')
bray_clusters <- hclust(bray_dist_mat, 'aver')
heatmap(as.matrix(reyes_most_variable), Rowv=as.dendrogram(bray_clusters), Colv=NA, col=jBuPuPalette)
# Clustering rows AND columns too
bray_dist_mat_t <- vegdist(t(reyes_most_variable), method='bray')
bray_clusters_t <- hclust(bray_dist_mat_t, 'aver')
heatmap(as.matrix(reyes_most_variable),
        Rowv=as.dendrogram(bray_clusters),
        Colv=as.dendrogram(bray_clusters_t),
        col=jBuPuPalette)



# ----------------------------------------------------------------------
# HEATMAP.2()
# ----------------------------------------------------------------------



# gplots' heatmap.2() is more commonly used (AFAIAA)
#   Can indicate experimental groups / covariates
#   Histograms
# Assign each row to one of two made-up categories
example_categories <- round(runif(n = 55, min = 1, max = 2)) # 1 or 2
example_categories <- replace(example_categories, which(example_categories == 1), "lightgreen")
example_categories <- replace(example_categories, which(example_categories == 2), "darkgreen")
cbind(row.names(reyes_most_variable), example_categories)
# Plot
library(gplots)
heatmap.2(as.matrix(reyes_most_variable),
          Rowv=as.dendrogram(bray_clusters),
          Colv=as.dendrogram(bray_clusters_t),
          RowSideColors=example_categories,
          trace='none',
          col=jBuPuPalette)




# ----------------------------------------------------------------------
# GGPLOT2()
# ----------------------------------------------------------------------
# Basic heatmap using geom_tile()


library(ggplot2)
library(reshape)
# library(stringr)

# Wrangling... ggplot likes melted dataframes 
reyes_most_variable_m <- melt(as.matrix(reyes_most_variable))
# reyes_most_variable_m$timepoint <- str_sub(reyes_most_variable_m$X2,-1,-1)
# reyes_most_variable_m$individual <- str_sub(reyes_most_variable_m$X2,1,-3)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

gg <- (ggplot(reyes_most_variable_m,
              aes(x = X2, y = X1, fill = value))
              + geom_tile()
              + scale_fill_gradientn(colours = myPalette(100))
              + scale_x_discrete(expand = c(0, 0))
              + scale_y_discrete(expand = c(0, 0))
              + coord_equal()
              + theme_bw())
print(gg)

# Mention phyloseq, heatplus
# Credit to http://www.molecularecologist.com/2013/08/making-heatmaps-with-r-for-microbiome-analysis/


