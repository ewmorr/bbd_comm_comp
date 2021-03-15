require(tidyverse)
source("~/ggplot_theme.txt")

source("~/repo/bbd_comm_comp/sum_trees/read_ASV_dat.LULU_tab.r")


#remove samples with less than 1000 seqs
full_metadata.ge1K = full_metadata %>%
filter(total_seqs >= 1000)

climate_vars = full_metadata.ge1K %>%
    dplyr::select(
    "HDD4.mean_nongrowing", "freezeThaw.mean_nongrowing", "ppt.mean_nongrowing", "HDD4.mean_growing", "freezeThaw.mean_growing","ppt.mean_growing"
)
tree_vars = full_metadata.ge1K %>%
dplyr::select(
"RaisedCanker", "dbh", "Wax", "TreeCond", "NeoFruiting"
)

#Scale the data (this is also perormed for the HMSC modeling to standardize variance)
climate_vars = apply(climate_vars,2,scale)
tree_vars = apply(tree_vars,2,scale)

#########
#run PCA#

#run the PCA on normalized data
climate.pca = princomp(climate_vars)
tree.pca = princomp(tree_vars)

#look at the result. The "Proportion variance explained" row indicates variance explained by each axis (i.e. percent variance explained)
summary(climate.pca)
summary(tree.pca)

#The first three and four axes explain >92% variance
#For climate cumu variance Comp.1 == 0.4884098, 2 == 0.7937840, 3 == 0.9209703. 4 == 0.99053291, 5 == 0.997836893, 6 == 1.000000000
#For tree vars cumu variance Comp.1 == 0.3563361, 2 == 0.6084786, 3 == 0.7979404, 4 == 0.9267141, 5 == 1.00000000

#can also be veiwed using the plot function
plot(climate.pca)
plot(tree.pca)

#and a simple biplot can be produced
biplot(climate.pca)
biplot(tree.pca)

########
#ggplot#

#### I like to plot with ggplot, and to do so requires extracting the axis scores. You would use the same approach to pull out axis scores for regression
#for example the scores are here

#to make a plottable data.frame including your metadata
climate.pca.df = data.frame(sample = full_metadata.ge1K$sample, PC1 = climate.pca$scores[,1], PC2 = climate.pca$scores[,2], PC3 = climate.pca$scores[,3], PC4 = climate.pca$scores[,4])
climate.pca_loadings = data.frame(sample = full_metadata.ge1K$sample, label = names(climate.pca$loadings[,1]), PC1 = climate.pca$loadings[,1], PC2 = climate.pca$loadings[,2], PC3 = climate.pca$loadings[,3], PC4 = climate.pca$loadings[,4] )

tree.pca.df = data.frame(sample = full_metadata.ge1K$sample, PC1 = tree.pca$scores[,1], PC2 = tree.pca$scores[,2], PC3 = tree.pca$scores[,3], PC4 = tree.pca$scores[,4])
tree.pca_loadings = data.frame(sample = full_metadata.ge1K$sample, label = names(tree.pca$loadings[,1]), PC1 = tree.pca$loadings[,1], PC2 = tree.pca$loadings[,2], PC3 = tree.pca$loadings[,3], PC4 = tree.pca$loadings[,4] )

saveRDS(climate.pca.df, "intermediate_RDS/climate.pca.df.rds")
saveRDS(climate.pca_loadings, "intermediate_RDS/climate.pca_loadings.rds")
saveRDS(tree.pca.df, "intermediate_RDS/tree.pca.df.rds")
saveRDS(tree.pca_loadings, "intermediate_RDS/tree.pca_loadings.rds")



pdf("covariate_PCAs/climate.pdf", width = 8, height = 6)

ggplot() +
geom_point(data = climate.pca.df, aes(x = PC1, y = PC2)) +
geom_text(data = climate.pca_loadings, aes(x = PC1*2.5, y = PC2*2.5, label = label)) +
geom_segment(data = climate.pca_loadings, aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2)) +
my_gg_theme

ggplot() +
geom_point(data = climate.pca.df, aes(x = PC1, y = PC3)) +
geom_text(data = climate.pca_loadings, aes(x = PC1*2.5, y = PC3*2.5, label = label)) +
geom_segment(data = climate.pca_loadings, aes(x = 0, y = 0, xend = PC1*2, yend = PC3*2)) +
my_gg_theme

ggplot() +
geom_point(data = climate.pca.df, aes(x = PC2, y = PC3)) +
geom_text(data = climate.pca_loadings, aes(x = PC2*2.5, y = PC3*2.5, label = label)) +
geom_segment(data = climate.pca_loadings, aes(x = 0, y = 0, xend = PC2*2, yend = PC3*2)) +
my_gg_theme

ggplot() +
geom_point(data = climate.pca.df, aes(x = PC1, y = PC4)) +
geom_text(data = climate.pca_loadings, aes(x = PC1*2.5, y = PC4*2.5, label = label)) +
geom_segment(data = climate.pca_loadings, aes(x = 0, y = 0, xend = PC1*2, yend = PC4*2)) +
my_gg_theme

ggplot() +
geom_point(data = climate.pca.df, aes(x = PC2, y = PC4)) +
geom_text(data = climate.pca_loadings, aes(x = PC2*2.5, y = PC4*2.5, label = label)) +
geom_segment(data = climate.pca_loadings, aes(x = 0, y = 0, xend = PC2*2, yend = PC4*2)) +
my_gg_theme

ggplot() +
geom_point(data = climate.pca.df, aes(x = PC3, y = PC4)) +
geom_text(data = climate.pca_loadings, aes(x = PC3*2.5, y = PC4*2.5, label = label)) +
geom_segment(data = climate.pca_loadings, aes(x = 0, y = 0, xend = PC3*2, yend = PC4*2)) +
my_gg_theme

dev.off()


pdf("covariate_PCAs/tree.pdf", width = 8, height = 6)

ggplot() +
geom_point(data = tree.pca.df, aes(x = PC1, y = PC2)) +
geom_text(data = tree.pca_loadings, aes(x = PC1*2.5, y = PC2*2.5, label = label)) +
geom_segment(data = tree.pca_loadings, aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2)) +
my_gg_theme

ggplot() +
geom_point(data = tree.pca.df, aes(x = PC1, y = PC3)) +
geom_text(data = tree.pca_loadings, aes(x = PC1*2.5, y = PC3*2.5, label = label)) +
geom_segment(data = tree.pca_loadings, aes(x = 0, y = 0, xend = PC1*2, yend = PC3*2)) +
my_gg_theme

ggplot() +
geom_point(data = tree.pca.df, aes(x = PC2, y = PC3)) +
geom_text(data = tree.pca_loadings, aes(x = PC2*2.5, y = PC3*2.5, label = label)) +
geom_segment(data = tree.pca_loadings, aes(x = 0, y = 0, xend = PC2*2, yend = PC3*2)) +
my_gg_theme

ggplot() +
geom_point(data = tree.pca.df, aes(x = PC1, y = PC4)) +
geom_text(data = tree.pca_loadings, aes(x = PC1*2.5, y = PC4*2.5, label = label)) +
geom_segment(data = tree.pca_loadings, aes(x = 0, y = 0, xend = PC1*2, yend = PC4*2)) +
my_gg_theme

ggplot() +
geom_point(data = tree.pca.df, aes(x = PC2, y = PC4)) +
geom_text(data = tree.pca_loadings, aes(x = PC2*2.5, y = PC4*2.5, label = label)) +
geom_segment(data = tree.pca_loadings, aes(x = 0, y = 0, xend = PC2*2, yend = PC4*2)) +
my_gg_theme

ggplot() +
geom_point(data = tree.pca.df, aes(x = PC3, y = PC4)) +
geom_text(data = tree.pca_loadings, aes(x = PC3*2.5, y = PC4*2.5, label = label)) +
geom_segment(data = tree.pca_loadings, aes(x = 0, y = 0, xend = PC3*2, yend = PC4*2)) +
my_gg_theme

dev.off()

