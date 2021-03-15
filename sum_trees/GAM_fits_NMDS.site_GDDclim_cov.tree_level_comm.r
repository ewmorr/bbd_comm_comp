require(tidyverse)
require(vegan)
require(gridExtra)
require(RColorBrewer)

source("~/ggplot_theme.txt")

#read data
source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")
#this pulls in objects:
#asv_tab
#asv_tax
#id_bench_map

#joins metadata files to get metadata object:
#full_metadata

#creates negatives and controls only asv_tab (long format):
#asv_tab.negatives.long

#creates new object with lowest informative taxonomic level and a character asv_tax with unknown instead of NA
#asv_informative_taxa
#asv_tax.char

#and calc duration_infectiones neonectria occurence in objects:
#Nf_v_Nd.long.metadata
#Nf_v_Nd.bin.metadata

#######
#PLOTS#
######
#NMDS#

##########################
#1000 seqs per sample min#
##########################

#rarefied table
asv_tab.gt1K.rare = readRDS(file = "intermediate_RDS/asv_tab.gt1K.rare.tree_sum.rds")
#NMDS
asv_tab.gt1K.rare.mds = readRDS(file = "intermediate_RDS/asv_tab.gt1K.rare.mds.tree_sum.rds")


#Add metadata for ordisruf
asv_tab.gt1K.rare.mds.metadata = left_join(
data.frame(sample = rownames(asv_tab.gt1K.rare.mds$points), asv_tab.gt1K.rare.mds$points),
full_metadata, by = "sample")

asv_tab.gt1K.rare.mds.metadata.neoOcurence = left_join(asv_tab.gt1K.rare.mds.metadata, Nf_v_Nd.bin, by = "sample") %>% left_join(., data.frame(sample = names(colSums(asv_tab)), total_seqs = colSums(asv_tab)))

asv_tab.gt1K.rare.mds.metadata.neoOcurence$Site = factor(asv_tab.gt1K.rare.mds.metadata.neoOcurence$Site, levels = c("MEN1", "MES1", "ADN1", "ADS1", "CW1", "TSP1", "GK1j", "ASH2", "MI1", "WF1"))

surf_vars = asv_tab.gt1K.rare.mds.metadata.neoOcurence %>% select(ends_with("growing"))
surf_vars[, 25:28] = c(surf_vars[,1]/surf_vars[7], surf_vars[,2]/surf_vars[8], surf_vars[,5]/surf_vars[7], surf_vars[,6]/surf_vars[8])
surf_vars[, 29:32] = c(surf_vars[,9]/surf_vars[21], surf_vars[,10]/surf_vars[22], surf_vars[,17]/surf_vars[21], surf_vars[,18]/surf_vars[22])

var_names = c(
    "Growing season GDD4",
    "Non-growing season GDD4",
    "Growing season freeze-thaw",
    "Non-growing season freeze-thaw",
    "Growing season precipitation",
    "Non-growing season precipitation",
    "Growing season length",
    "Non-growing season length",
    "Growing season GDD4 10-yr mean",
    "Non-growing season GDD4 10-yr mean",
    "Growing season freeze-thaw 10-yr mean",
    "Non-growing season freeze-thaw 10-yr mean",
    "Growing season GDD4 10-yr st. dev.",
    "Non-growing season GDD4 10-yr st. dev.",
    "Growing season freeze-thaw 10-yr st. dev.",
    "Non-growing season freeze-thaw 10-yr st. dev.",
    "Growing season precipitation 10-yr mean",
    "Non-growing season precipitation 10-yr mean",
    "Growing season precipitation 10-yr st. dev.",
    "Non-growing season precipitation 10-yr st. dev.",
    "Growing season length 10-yr mean",
    "Non-growing season length 10-yr mean",
    "Growing season length 10-yr st. dev.",
    "Non-growing season length 10-yr st. dev.",
    "Growing season GDD4 per day",
    "Non-growing season GDD4 per day",
    "Growing precipitation per day",
    "Non-growing precipitation per day",
    "Growing season GDD4 per day 10-yr mean",
    "Non-growing season GDD4 per day 10-yr mean",
    "Growing precipitation per day 10-yr mean",
    "Non-growing precipitation per day 10-yr mean"
)

#loop all of this
surf_list = vector(mode = "list", length = length(surf_vars))

for(i in 1:length(surf_vars)){
    surf_list[[i]]$ordisurf = ordisurf(asv_tab.gt1K.rare.mds ~ surf_vars[,i])
    surf_list[[i]]$anova = anova(surf_list[[i]]$ordisurf)
    surf_list[[i]]$coefs$r2 = signif(surf_list[[i]]$anova$r.sq, 2)
     surf_list[[i]]$coefs$p = signif(surf_list[[i]]$anova$s.pv,2)
    surf_list[[i]]$grid = expand.grid(
        x = surf_list[[i]]$ordisurf$grid$x,
        y = surf_list[[i]]$ordisurf$grid$y
    )
    surf_list[[i]]$grid$z = as.vector(surf_list[[i]]$ordisurf$grid$z)
    surf_list[[i]]$grid = data.frame(na.omit(surf_list[[i]]$grid))
}

coefs_tab = data.frame(
    variable = vector(mode = "character", length = length(surf_list)),
    r2 = vector(mode = "numeric", length = length(surf_list)),
    p.val = vector(mode = "numeric", length = length(surf_list)),
    stringsAsFactors = F
)

for(i in 1:length(surf_list)){
    coefs_tab$variable[i] = var_names[i]
    coefs_tab$r2[i] = surf_list[[i]]$coefs$r2
    coefs_tab$p.val[i] = surf_list[[i]]$coefs$p
}

write.table(coefs_tab, file = "GAM_fits/GDD_coefs.txt", row.names = F, quote = F, sep = "\t")


plot_list = vector(mode = "list", length = length(surf_vars))

for(i in 1:length(surf_list)){
    plot_list[[i]] = ggplot() +
    geom_point(data = asv_tab.gt1K.rare.mds.metadata.neoOcurence, aes(MDS1, MDS2, fill = Site), size = 1.5, shape = 21) +
    stat_contour(data = surf_list[[i]]$grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
    my_gg_theme +
    theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
    scale_color_gradient(name = "GAM val", low = "light grey", high = "black") +
    scale_fill_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI"), guide = F) +
    labs(
        title = paste(var_names[i],"\nR2 =", surf_list[[i]]$coefs$r2,", P =", surf_list[[i]]$coefs$p ),
        x = "NMDS1",
        y = "NMDS2"
    )

}

pdf("GAM_fits/NMDS_1Kmin_GAM_GDD.sum_trees.pdf", width = 50, height = 16)
grid.arrange(
    plot_list[[1]],
    plot_list[[2]],
    plot_list[[3]],
    plot_list[[4]],
    plot_list[[5]],
    plot_list[[6]],
    plot_list[[7]],

    plot_list[[9]],
    plot_list[[10]],
    plot_list[[11]],
    plot_list[[12]],
    plot_list[[17]],
    plot_list[[18]],
    plot_list[[21]],

    plot_list[[13]],
    plot_list[[14]],
    plot_list[[15]],
    plot_list[[16]],
    plot_list[[19]],
    plot_list[[20]],
    plot_list[[23]],
    nrow = 3
)
dev.off()


pdf("GAM_fits/NMDS_1Kmin_GAM_GDD.per_day.sum_trees.pdf", width = 30, height = 12)
grid.arrange(
plot_list[[25]],
plot_list[[26]],
plot_list[[27]],
plot_list[[28]],
plot_list[[29]],
plot_list[[30]],
plot_list[[31]],
plot_list[[32]],
nrow = 2
)
dev.off()



#################################
#1000 seqs per sample min binary#
#################################

#rarefied table
asv_tab.gt1K.rare = readRDS(file = "intermediate_RDS/asv_tab.gt1K.rare.tree_sum.rds")
#NMDS
asv_tab.gt1K.rare.mds.bin = readRDS(file = "intermediate_RDS/asv_tab.gt1K.rare.mds.bin.tree_sum.rds")


#Add metadata for ordisruf
asv_tab.gt1K.rare.mds.metadata = left_join(
data.frame(sample = rownames(asv_tab.gt1K.rare.mds.bin$points), asv_tab.gt1K.rare.mds.bin$points),
full_metadata, by = "sample")

asv_tab.gt1K.rare.mds.metadata.neoOcurence = left_join(asv_tab.gt1K.rare.mds.metadata, Nf_v_Nd.bin, by = "sample") %>% left_join(., data.frame(sample = names(colSums(asv_tab)), total_seqs = colSums(asv_tab)))

asv_tab.gt1K.rare.mds.metadata.neoOcurence$Site = factor(asv_tab.gt1K.rare.mds.metadata.neoOcurence$Site, levels = c("MEN1", "MES1", "ADN1", "ADS1", "CW1", "TSP1", "GK1j", "ASH2", "MI1", "WF1"))

#surf_vars and names defined above

#loop all of this
surf_list = vector(mode = "list", length = length(surf_vars))

for(i in 1:length(surf_vars)){
    surf_list[[i]]$ordisurf = ordisurf(asv_tab.gt1K.rare.mds.bin ~ surf_vars[,i])
    surf_list[[i]]$anova = anova(surf_list[[i]]$ordisurf)
    surf_list[[i]]$coefs$r2 = signif(surf_list[[i]]$anova$r.sq, 2)
    surf_list[[i]]$coefs$p = signif(surf_list[[i]]$anova$s.pv,2)
    surf_list[[i]]$grid = expand.grid(
    x = surf_list[[i]]$ordisurf$grid$x,
    y = surf_list[[i]]$ordisurf$grid$y
    )
    surf_list[[i]]$grid$z = as.vector(surf_list[[i]]$ordisurf$grid$z)
    surf_list[[i]]$grid = data.frame(na.omit(surf_list[[i]]$grid))
}

coefs_tab = data.frame(
variable = vector(mode = "character", length = length(surf_list)),
r2 = vector(mode = "numeric", length = length(surf_list)),
p.val = vector(mode = "numeric", length = length(surf_list)),
stringsAsFactors = F
)

for(i in 1:length(surf_vars)){
    coefs_tab$variable[i] = var_names[i]
    coefs_tab$r2[i] = surf_list[[i]]$coefs$r2
    coefs_tab$p.val[i] = surf_list[[i]]$coefs$p
}

write.table(coefs_tab, file = "GAM_fits/GDD_coefs.bin.txt", row.names = F, quote = F, sep = "\t")


plot_list = vector(mode = "list", length = length(surf_vars))

for(i in 1:length(surf_list)){
    plot_list[[i]] = ggplot() +
    geom_point(data = asv_tab.gt1K.rare.mds.metadata.neoOcurence, aes(MDS1, MDS2, fill = Site), size = 1.5, shape = 21) +
    stat_contour(data = surf_list[[i]]$grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
    my_gg_theme +
    theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
    scale_color_gradient(name = "GAM val", low = "light grey", high = "black") +
    scale_fill_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI"), guide = F) +
    labs(
    title = paste(var_names[i],"\nR2 =", surf_list[[i]]$coefs$r2,", P =", surf_list[[i]]$coefs$p ),
    x = "NMDS1",
    y = "NMDS2"
    )
    
}

#Site and neo occurence for reference

pdf("GAM_fits/NMDS_1Kmin-bin_GAM_GDD.sum_trees.pdf", width = 50, height = 16)
grid.arrange(
plot_list[[1]],
plot_list[[2]],
plot_list[[3]],
plot_list[[4]],
plot_list[[5]],
plot_list[[6]],
plot_list[[7]],

plot_list[[9]],
plot_list[[10]],
plot_list[[11]],
plot_list[[12]],
plot_list[[17]],
plot_list[[18]],
plot_list[[21]],

plot_list[[13]],
plot_list[[14]],
plot_list[[15]],
plot_list[[16]],
plot_list[[19]],
plot_list[[20]],
plot_list[[23]],
nrow = 3
)
dev.off()


pdf("GAM_fits/NMDS_1Kmin-bin_GAM_GDD.per_day.sum_trees.pdf", width = 30, height = 12)
grid.arrange(
plot_list[[25]],
plot_list[[26]],
plot_list[[27]],
plot_list[[28]],
plot_list[[29]],
plot_list[[30]],
plot_list[[31]],
plot_list[[32]],
nrow = 2
)
dev.off()
