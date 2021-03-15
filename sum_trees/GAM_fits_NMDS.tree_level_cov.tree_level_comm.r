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

surf_vars = asv_tab.gt1K.rare.mds.metadata.neoOcurence %>%
    select("dbh", "NeoFruiting", "RaisedCanker", "TarrySpots", "Wax", "Xylococcus", "TreeCond")

surf.dbh = ordisurf(asv_tab.gt1K.rare.mds ~ surf_vars[,1])
surf.neo_fruiting = ordisurf(asv_tab.gt1K.rare.mds ~ surf_vars[,2])
surf.cankers = ordisurf(asv_tab.gt1K.rare.mds ~ surf_vars[,3])
surf.tarry_spots = ordisurf(asv_tab.gt1K.rare.mds ~ surf_vars[,4])
surf.wax = ordisurf(asv_tab.gt1K.rare.mds ~ surf_vars[,5])
surf.xylococcus = ordisurf(asv_tab.gt1K.rare.mds ~ surf_vars[,6])
surf.tree_cond = ordisurf(asv_tab.gt1K.rare.mds ~ surf_vars[,7])

surf.dbh.anova = anova(surf.dbh)
surf.neo_fruiting.anova = anova(surf.neo_fruiting)
surf.cankers.anova = anova(surf.cankers)
surf.tarry_spots.anova = anova(surf.tarry_spots)
surf.wax.anova = anova(surf.wax)
surf.xylococcus.anova = anova(surf.xylococcus)
surf.tree_cond.anova = anova(surf.tree_cond)

#Extract plottable

#dbh
surf.dbh.grid = expand.grid(x = surf.dbh$grid$x, y = surf.dbh$grid$y)
surf.dbh.grid$z = as.vector(surf.dbh$grid$z)
surf.dbh.grid = data.frame(na.omit(surf.dbh.grid))

p.dbh = ggplot() +
geom_point(data = asv_tab.gt1K.rare.mds.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = dbh), size = 2.5) +
stat_contour(data = surf.dbh.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "DBH GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "DBH obs", low = "light grey", high = "black") +
labs(title = paste("DBH (cm)\nR2 =",signif(surf.dbh.anova$r.sq, 2),", P =", signif(surf.dbh.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#neo_fruiting
surf.neo_fruiting.grid = expand.grid(x = surf.neo_fruiting$grid$x, y = surf.neo_fruiting$grid$y)
surf.neo_fruiting.grid$z = as.vector(surf.neo_fruiting$grid$z)
surf.neo_fruiting.grid = data.frame(na.omit(surf.neo_fruiting.grid))

p.neo_fruiting = ggplot() +
geom_point(data = asv_tab.gt1K.rare.mds.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = NeoFruiting), size = 2.5) +
stat_contour(data = surf.neo_fruiting.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "neo_fruiting GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "neo_fruiting obs", low = "light grey", high = "black") +
labs(title = paste("Neonectria fruiting (0-5)\nR2 =",signif(surf.neo_fruiting.anova$r.sq, 2),", P =", signif(surf.neo_fruiting.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#cankers
surf.cankers.grid = expand.grid(x = surf.cankers$grid$x, y = surf.cankers$grid$y)
surf.cankers.grid$z = as.vector(surf.cankers$grid$z)
surf.cankers.grid = data.frame(na.omit(surf.cankers.grid))

p.cankers = ggplot() +
geom_point(data = asv_tab.gt1K.rare.mds.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = RaisedCanker), size = 2.5) +
stat_contour(data = surf.cankers.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "cankers GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "cankers obs", low = "light grey", high = "black") +
labs(title = paste("Cankers\nR2 =",signif(surf.cankers.anova$r.sq, 2),", P =", signif(surf.cankers.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#tarry_spots
surf.tarry_spots.grid = expand.grid(x = surf.tarry_spots$grid$x, y = surf.tarry_spots$grid$y)
surf.tarry_spots.grid$z = as.vector(surf.tarry_spots$grid$z)
surf.tarry_spots.grid = data.frame(na.omit(surf.tarry_spots.grid))

p.tarry_spots = ggplot() +
geom_point(data = asv_tab.gt1K.rare.mds.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = TarrySpots), size = 2.5) +
stat_contour(data = surf.tarry_spots.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "tarry_spots GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "tarry_spots obs", low = "light grey", high = "black") +
labs(title = paste("Tarry spots\nR2 =",signif(surf.tarry_spots.anova$r.sq, 2),", P =", signif(surf.tarry_spots.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#wax
surf.wax.grid = expand.grid(x = surf.wax$grid$x, y = surf.wax$grid$y)
surf.wax.grid$z = as.vector(surf.wax$grid$z)
surf.wax.grid = data.frame(na.omit(surf.wax.grid))

p.wax = ggplot() +
geom_point(data = asv_tab.gt1K.rare.mds.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = Wax), size = 2.5) +
stat_contour(data = surf.wax.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "wax GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "wax obs", low = "light grey", high = "black") +
labs(title = paste("Wax (0-5)\nR2 =",signif(surf.wax.anova$r.sq, 2),", P =", signif(surf.wax.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#xylococcus
surf.xylococcus.grid = expand.grid(x = surf.xylococcus$grid$x, y = surf.xylococcus$grid$y)
surf.xylococcus.grid$z = as.vector(surf.xylococcus$grid$z)
surf.xylococcus.grid = data.frame(na.omit(surf.xylococcus.grid))

p.xylococcus = ggplot() +
geom_point(data = asv_tab.gt1K.rare.mds.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = Xylococcus), size = 2.5) +
stat_contour(data = surf.xylococcus.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "xylococcus GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "xylococcus obs", low = "light grey", high = "black") +
labs(title = paste("Xylococcus\nR2 =",signif(surf.xylococcus.anova$r.sq, 2),", P =", signif(surf.xylococcus.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#tree_cond
surf.tree_cond.grid = expand.grid(x = surf.tree_cond$grid$x, y = surf.tree_cond$grid$y)
surf.tree_cond.grid$z = as.vector(surf.tree_cond$grid$z)
surf.tree_cond.grid = data.frame(na.omit(surf.tree_cond.grid))

p.tree_cond = ggplot() +
geom_point(data = asv_tab.gt1K.rare.mds.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = TreeCond), size = 2.5) +
stat_contour(data = surf.tree_cond.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "tree_cond GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "tree_cond obs", low = "light grey", high = "black") +
labs(title = paste("Tree condition (0-5)\nR2 =",signif(surf.tree_cond.anova$r.sq, 2),", P =", signif(surf.tree_cond.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#Site and neo occurence for reference

p.site = ggplot() +
geom_point(data = asv_tab.gt1K.rare.mds.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = Site), size = 3) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both")) +
scale_color_gradient(name = "ppt", low = "light grey", high = "black") +
scale_fill_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")) +
labs(title = "Sites and\nneonectria occurence", x = "NMDS1", y = "NMDS2") +
guides(fill = guide_legend(override.aes=list(shape=21,size=3)), shape = guide_legend(override.aes=list(size=3)))


pdf("GAM_fits/NMDS_neo_occurence_by_disease_severity_1Kmin_GAM.sum_trees.tree_level_cov.pdf", width = 22, height = 18)
grid.arrange(p.site,p.dbh,p.neo_fruiting,p.wax,p.cankers,p.tarry_spots,p.xylococcus,p.tree_cond,nrow = 3)
dev.off()




#################################
#1000 seqs per sample min binary#
#################################

#rarefied table
asv_tab.gt1K.rare = readRDS(file = "intermediate_RDS/asv_tab.gt1K.rare.tree_sum.rds")
#NMDS
asv_tab.gt1K.rare.mds.bin = readRDS(file = "intermediate_RDS/asv_tab.gt1K.rare.mds.bin.tree_sum.rds")


#Add metadata for ordisruf
asv_tab.gt1K.rare.mds.bin.metadata = left_join(
data.frame(sample = rownames(asv_tab.gt1K.rare.mds.bin$points), asv_tab.gt1K.rare.mds.bin$points),
full_metadata, by = "sample")

asv_tab.gt1K.rare.mds.bin.metadata.neoOcurence = left_join(asv_tab.gt1K.rare.mds.bin.metadata, Nf_v_Nd.bin, by = "sample") %>% left_join(., data.frame(sample = names(colSums(asv_tab)), total_seqs = colSums(asv_tab)))

asv_tab.gt1K.rare.mds.bin.metadata.neoOcurence$Site = factor(asv_tab.gt1K.rare.mds.bin.metadata.neoOcurence$Site, levels = c("MEN1", "MES1", "ADN1", "ADS1", "CW1", "TSP1", "GK1j", "ASH2", "MI1", "WF1"))

surf_vars = asv_tab.gt1K.rare.mds.metadata.neoOcurence %>%
select("dbh", "NeoFruiting", "RaisedCanker", "TarrySpots", "Wax", "Xylococcus", "TreeCond")

surf.dbh = ordisurf(asv_tab.gt1K.rare.mds.bin ~ surf_vars[,1])
surf.neo_fruiting = ordisurf(asv_tab.gt1K.rare.mds.bin ~ surf_vars[,2])
surf.cankers = ordisurf(asv_tab.gt1K.rare.mds.bin ~ surf_vars[,3])
surf.tarry_spots = ordisurf(asv_tab.gt1K.rare.mds.bin ~ surf_vars[,4])
surf.wax = ordisurf(asv_tab.gt1K.rare.mds.bin ~ surf_vars[,5])
surf.xylococcus = ordisurf(asv_tab.gt1K.rare.mds.bin ~ surf_vars[,6])
surf.tree_cond = ordisurf(asv_tab.gt1K.rare.mds.bin ~ surf_vars[,7])

surf.dbh.anova = anova(surf.dbh)
surf.neo_fruiting.anova = anova(surf.neo_fruiting)
surf.cankers.anova = anova(surf.cankers)
surf.tarry_spots.anova = anova(surf.tarry_spots)
surf.wax.anova = anova(surf.wax)
surf.xylococcus.anova = anova(surf.xylococcus)
surf.tree_cond.anova = anova(surf.tree_cond)

#Extract plottable

#dbh
surf.dbh.grid = expand.grid(x = surf.dbh$grid$x, y = surf.dbh$grid$y)
surf.dbh.grid$z = as.vector(surf.dbh$grid$z)
surf.dbh.grid = data.frame(na.omit(surf.dbh.grid))

p.dbh = ggplot() +
geom_point(data = asv_tab.gt1K.rare.mds.bin.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = dbh), size = 2.5) +
stat_contour(data = surf.dbh.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "DBH GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "DBH obs", low = "light grey", high = "black") +
labs(title = paste("DBH (cm)\nR2 =",signif(surf.dbh.anova$r.sq, 2),", P =", signif(surf.dbh.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#neo_fruiting
surf.neo_fruiting.grid = expand.grid(x = surf.neo_fruiting$grid$x, y = surf.neo_fruiting$grid$y)
surf.neo_fruiting.grid$z = as.vector(surf.neo_fruiting$grid$z)
surf.neo_fruiting.grid = data.frame(na.omit(surf.neo_fruiting.grid))

p.neo_fruiting = ggplot() +
geom_point(data = asv_tab.gt1K.rare.mds.bin.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = NeoFruiting), size = 2.5) +
stat_contour(data = surf.neo_fruiting.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "neo_fruiting GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "neo_fruiting obs", low = "light grey", high = "black") +
labs(title = paste("Neonectria fruiting (0-5)\nR2 =",signif(surf.neo_fruiting.anova$r.sq, 2),", P =", signif(surf.neo_fruiting.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#cankers
surf.cankers.grid = expand.grid(x = surf.cankers$grid$x, y = surf.cankers$grid$y)
surf.cankers.grid$z = as.vector(surf.cankers$grid$z)
surf.cankers.grid = data.frame(na.omit(surf.cankers.grid))

p.cankers = ggplot() +
geom_point(data = asv_tab.gt1K.rare.mds.bin.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = RaisedCanker), size = 2.5) +
stat_contour(data = surf.cankers.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "cankers GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "cankers obs", low = "light grey", high = "black") +
labs(title = paste("Cankers\nR2 =",signif(surf.cankers.anova$r.sq, 2),", P =", signif(surf.cankers.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#tarry_spots
surf.tarry_spots.grid = expand.grid(x = surf.tarry_spots$grid$x, y = surf.tarry_spots$grid$y)
surf.tarry_spots.grid$z = as.vector(surf.tarry_spots$grid$z)
surf.tarry_spots.grid = data.frame(na.omit(surf.tarry_spots.grid))

p.tarry_spots = ggplot() +
geom_point(data = asv_tab.gt1K.rare.mds.bin.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = TarrySpots), size = 2.5) +
stat_contour(data = surf.tarry_spots.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "tarry_spots GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "tarry_spots obs", low = "light grey", high = "black") +
labs(title = paste("Tarry spots\nR2 =",signif(surf.tarry_spots.anova$r.sq, 2),", P =", signif(surf.tarry_spots.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#wax
surf.wax.grid = expand.grid(x = surf.wax$grid$x, y = surf.wax$grid$y)
surf.wax.grid$z = as.vector(surf.wax$grid$z)
surf.wax.grid = data.frame(na.omit(surf.wax.grid))

p.wax = ggplot() +
geom_point(data = asv_tab.gt1K.rare.mds.bin.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = Wax), size = 2.5) +
stat_contour(data = surf.wax.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "wax GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "wax obs", low = "light grey", high = "black") +
labs(title = paste("Wax (0-5)\nR2 =",signif(surf.wax.anova$r.sq, 2),", P =", signif(surf.wax.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#xylococcus
surf.xylococcus.grid = expand.grid(x = surf.xylococcus$grid$x, y = surf.xylococcus$grid$y)
surf.xylococcus.grid$z = as.vector(surf.xylococcus$grid$z)
surf.xylococcus.grid = data.frame(na.omit(surf.xylococcus.grid))

p.xylococcus = ggplot() +
geom_point(data = asv_tab.gt1K.rare.mds.bin.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = Xylococcus), size = 2.5) +
stat_contour(data = surf.xylococcus.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "xylococcus GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "xylococcus obs", low = "light grey", high = "black") +
labs(title = paste("Xylococcus\nR2 =",signif(surf.xylococcus.anova$r.sq, 2),", P =", signif(surf.xylococcus.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#tree_cond
surf.tree_cond.grid = expand.grid(x = surf.tree_cond$grid$x, y = surf.tree_cond$grid$y)
surf.tree_cond.grid$z = as.vector(surf.tree_cond$grid$z)
surf.tree_cond.grid = data.frame(na.omit(surf.tree_cond.grid))

p.tree_cond = ggplot() +
geom_point(data = asv_tab.gt1K.rare.mds.bin.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = TreeCond), size = 2.5) +
stat_contour(data = surf.tree_cond.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "tree_cond GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "tree_cond obs", low = "light grey", high = "black") +
labs(title = paste("Tree condition (0-5)\nR2 =",signif(surf.tree_cond.anova$r.sq, 2),", P =", signif(surf.tree_cond.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#Site and neo occurence for reference

p.site = ggplot() +
geom_point(data = asv_tab.gt1K.rare.mds.bin.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = Site), size = 3) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both")) +
scale_color_gradient(name = "ppt", low = "light grey", high = "black") +
scale_fill_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")) +
labs(title = "Sites and\nneonectria occurence", x = "NMDS1", y = "NMDS2") +
guides(fill = guide_legend(override.aes=list(shape=21,size=3)), shape = guide_legend(override.aes=list(size=3)))


pdf("GAM_fits/NMDS.bin_neo_occurence_by_disease_severity_1Kmin_BINARY_GAM.sum_trees.tree_level_cov.pdf", width = 22, height = 18)
grid.arrange(p.site,p.dbh,p.neo_fruiting,p.wax,p.cankers,p.tarry_spots,p.xylococcus,p.tree_cond,nrow = 3)
dev.off()


##########################
#5000 seqs per sample min#
##########################

#rarefied table
asv_tab.gt5K.rare = readRDS(file = "intermediate_RDS/asv_tab.gt5K.rare.tree_sum.rds")
#NMDS
asv_tab.gt5K.rare.mds = readRDS(file = "intermediate_RDS/asv_tab.gt5K.rare.mds.tree_sum.rds")


#Add metadata for ordisruf
asv_tab.gt5K.rare.mds.metadata = left_join(
data.frame(sample = rownames(asv_tab.gt5K.rare.mds$points), asv_tab.gt5K.rare.mds$points),
full_metadata, by = "sample")

asv_tab.gt5K.rare.mds.metadata.neoOcurence = left_join(asv_tab.gt5K.rare.mds.metadata, Nf_v_Nd.bin, by = "sample") %>% left_join(., data.frame(sample = names(colSums(asv_tab)), total_seqs = colSums(asv_tab)))

asv_tab.gt5K.rare.mds.metadata.neoOcurence$Site = factor(asv_tab.gt5K.rare.mds.metadata.neoOcurence$Site, levels = c("MEN1", "MES1", "ADN1", "ADS1", "CW1", "TSP1", "GK1j", "ASH2", "MI1", "WF1"))

surf_vars = asv_tab.gt5K.rare.mds.metadata.neoOcurence %>%
select("dbh", "NeoFruiting", "RaisedCanker", "TarrySpots", "Wax", "Xylococcus", "TreeCond")

surf.dbh = ordisurf(asv_tab.gt5K.rare.mds ~ surf_vars[,1])
surf.neo_fruiting = ordisurf(asv_tab.gt5K.rare.mds ~ surf_vars[,2])
surf.cankers = ordisurf(asv_tab.gt5K.rare.mds ~ surf_vars[,3])
surf.tarry_spots = ordisurf(asv_tab.gt5K.rare.mds ~ surf_vars[,4])
surf.wax = ordisurf(asv_tab.gt5K.rare.mds ~ surf_vars[,5])
surf.xylococcus = ordisurf(asv_tab.gt5K.rare.mds ~ surf_vars[,6])
surf.tree_cond = ordisurf(asv_tab.gt5K.rare.mds ~ surf_vars[,7])

surf.dbh.anova = anova(surf.dbh)
surf.neo_fruiting.anova = anova(surf.neo_fruiting)
surf.cankers.anova = anova(surf.cankers)
surf.tarry_spots.anova = anova(surf.tarry_spots)
surf.wax.anova = anova(surf.wax)
surf.xylococcus.anova = anova(surf.xylococcus)
surf.tree_cond.anova = anova(surf.tree_cond)

#Extract plottable

#dbh
surf.dbh.grid = expand.grid(x = surf.dbh$grid$x, y = surf.dbh$grid$y)
surf.dbh.grid$z = as.vector(surf.dbh$grid$z)
surf.dbh.grid = data.frame(na.omit(surf.dbh.grid))

p.dbh = ggplot() +
geom_point(data = asv_tab.gt5K.rare.mds.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = dbh), size = 2.5) +
stat_contour(data = surf.dbh.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "DBH GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "DBH obs", low = "light grey", high = "black") +
labs(title = paste("DBH (cm)\nR2 =",signif(surf.dbh.anova$r.sq, 2),", P =", signif(surf.dbh.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#neo_fruiting
surf.neo_fruiting.grid = expand.grid(x = surf.neo_fruiting$grid$x, y = surf.neo_fruiting$grid$y)
surf.neo_fruiting.grid$z = as.vector(surf.neo_fruiting$grid$z)
surf.neo_fruiting.grid = data.frame(na.omit(surf.neo_fruiting.grid))

p.neo_fruiting = ggplot() +
geom_point(data = asv_tab.gt5K.rare.mds.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = NeoFruiting), size = 2.5) +
stat_contour(data = surf.neo_fruiting.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "neo_fruiting GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "neo_fruiting obs", low = "light grey", high = "black") +
labs(title = paste("Neonectria fruiting (0-5)\nR2 =",signif(surf.neo_fruiting.anova$r.sq, 2),", P =", signif(surf.neo_fruiting.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#cankers
surf.cankers.grid = expand.grid(x = surf.cankers$grid$x, y = surf.cankers$grid$y)
surf.cankers.grid$z = as.vector(surf.cankers$grid$z)
surf.cankers.grid = data.frame(na.omit(surf.cankers.grid))

p.cankers = ggplot() +
geom_point(data = asv_tab.gt5K.rare.mds.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = RaisedCanker), size = 2.5) +
stat_contour(data = surf.cankers.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "cankers GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "cankers obs", low = "light grey", high = "black") +
labs(title = paste("Cankers\nR2 =",signif(surf.cankers.anova$r.sq, 2),", P =", signif(surf.cankers.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#tarry_spots
surf.tarry_spots.grid = expand.grid(x = surf.tarry_spots$grid$x, y = surf.tarry_spots$grid$y)
surf.tarry_spots.grid$z = as.vector(surf.tarry_spots$grid$z)
surf.tarry_spots.grid = data.frame(na.omit(surf.tarry_spots.grid))

p.tarry_spots = ggplot() +
geom_point(data = asv_tab.gt5K.rare.mds.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = TarrySpots), size = 2.5) +
stat_contour(data = surf.tarry_spots.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "tarry_spots GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "tarry_spots obs", low = "light grey", high = "black") +
labs(title = paste("Tarry spots\nR2 =",signif(surf.tarry_spots.anova$r.sq, 2),", P =", signif(surf.tarry_spots.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#wax
surf.wax.grid = expand.grid(x = surf.wax$grid$x, y = surf.wax$grid$y)
surf.wax.grid$z = as.vector(surf.wax$grid$z)
surf.wax.grid = data.frame(na.omit(surf.wax.grid))

p.wax = ggplot() +
geom_point(data = asv_tab.gt5K.rare.mds.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = Wax), size = 2.5) +
stat_contour(data = surf.wax.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "wax GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "wax obs", low = "light grey", high = "black") +
labs(title = paste("Wax (0-5)\nR2 =",signif(surf.wax.anova$r.sq, 2),", P =", signif(surf.wax.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#xylococcus
surf.xylococcus.grid = expand.grid(x = surf.xylococcus$grid$x, y = surf.xylococcus$grid$y)
surf.xylococcus.grid$z = as.vector(surf.xylococcus$grid$z)
surf.xylococcus.grid = data.frame(na.omit(surf.xylococcus.grid))

p.xylococcus = ggplot() +
geom_point(data = asv_tab.gt5K.rare.mds.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = Xylococcus), size = 2.5) +
stat_contour(data = surf.xylococcus.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "xylococcus GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "xylococcus obs", low = "light grey", high = "black") +
labs(title = paste("Xylococcus\nR2 =",signif(surf.xylococcus.anova$r.sq, 2),", P =", signif(surf.xylococcus.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#tree_cond
surf.tree_cond.grid = expand.grid(x = surf.tree_cond$grid$x, y = surf.tree_cond$grid$y)
surf.tree_cond.grid$z = as.vector(surf.tree_cond$grid$z)
surf.tree_cond.grid = data.frame(na.omit(surf.tree_cond.grid))

p.tree_cond = ggplot() +
geom_point(data = asv_tab.gt5K.rare.mds.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = TreeCond), size = 2.5) +
stat_contour(data = surf.tree_cond.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "tree_cond GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "tree_cond obs", low = "light grey", high = "black") +
labs(title = paste("Tree condition (0-5)\nR2 =",signif(surf.tree_cond.anova$r.sq, 2),", P =", signif(surf.tree_cond.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#Site and neo occurence for reference

p.site = ggplot() +
geom_point(data = asv_tab.gt5K.rare.mds.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = Site), size = 3) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both")) +
scale_color_gradient(name = "ppt", low = "light grey", high = "black") +
scale_fill_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")) +
labs(title = "Sites and\nneonectria occurence", x = "NMDS1", y = "NMDS2") +
guides(fill = guide_legend(override.aes=list(shape=21,size=3)), shape = guide_legend(override.aes=list(size=3)))


pdf("GAM_fits/NMDS_neo_occurence_by_disease_severity_5Kmin_GAM.sum_trees.tree_level_cov.pdf", width = 22, height = 18)
grid.arrange(p.site,p.dbh,p.neo_fruiting,p.wax,p.cankers,p.tarry_spots,p.xylococcus,p.tree_cond,nrow = 3)
dev.off()




#################################
#5000 seqs per sample min binary#
#################################

#rarefied table
asv_tab.gt5K.rare = readRDS(file = "intermediate_RDS/asv_tab.gt5K.rare.tree_sum.rds")
#NMDS
asv_tab.gt5K.rare.mds.bin = readRDS(file = "intermediate_RDS/asv_tab.gt5K.rare.mds.bin.tree_sum.rds")


#Add metadata for ordisruf
asv_tab.gt5K.rare.mds.bin.metadata = left_join(
data.frame(sample = rownames(asv_tab.gt5K.rare.mds.bin$points), asv_tab.gt5K.rare.mds.bin$points),
full_metadata, by = "sample")

asv_tab.gt5K.rare.mds.bin.metadata.neoOcurence = left_join(asv_tab.gt5K.rare.mds.bin.metadata, Nf_v_Nd.bin, by = "sample") %>% left_join(., data.frame(sample = names(colSums(asv_tab)), total_seqs = colSums(asv_tab)))

asv_tab.gt5K.rare.mds.bin.metadata.neoOcurence$Site = factor(asv_tab.gt5K.rare.mds.bin.metadata.neoOcurence$Site, levels = c("MEN1", "MES1", "ADN1", "ADS1", "CW1", "TSP1", "GK1j", "ASH2", "MI1", "WF1"))

surf_vars = asv_tab.gt5K.rare.mds.metadata.neoOcurence %>%
select("dbh", "NeoFruiting", "RaisedCanker", "TarrySpots", "Wax", "Xylococcus", "TreeCond")

surf.dbh = ordisurf(asv_tab.gt5K.rare.mds.bin ~ surf_vars[,1])
surf.neo_fruiting = ordisurf(asv_tab.gt5K.rare.mds.bin ~ surf_vars[,2])
surf.cankers = ordisurf(asv_tab.gt5K.rare.mds.bin ~ surf_vars[,3])
surf.tarry_spots = ordisurf(asv_tab.gt5K.rare.mds.bin ~ surf_vars[,4])
surf.wax = ordisurf(asv_tab.gt5K.rare.mds.bin ~ surf_vars[,5])
surf.xylococcus = ordisurf(asv_tab.gt5K.rare.mds.bin ~ surf_vars[,6])
surf.tree_cond = ordisurf(asv_tab.gt5K.rare.mds.bin ~ surf_vars[,7])

surf.dbh.anova = anova(surf.dbh)
surf.neo_fruiting.anova = anova(surf.neo_fruiting)
surf.cankers.anova = anova(surf.cankers)
surf.tarry_spots.anova = anova(surf.tarry_spots)
surf.wax.anova = anova(surf.wax)
surf.xylococcus.anova = anova(surf.xylococcus)
surf.tree_cond.anova = anova(surf.tree_cond)

#Extract plottable

#dbh
surf.dbh.grid = expand.grid(x = surf.dbh$grid$x, y = surf.dbh$grid$y)
surf.dbh.grid$z = as.vector(surf.dbh$grid$z)
surf.dbh.grid = data.frame(na.omit(surf.dbh.grid))

p.dbh = ggplot() +
geom_point(data = asv_tab.gt5K.rare.mds.bin.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = dbh), size = 2.5) +
stat_contour(data = surf.dbh.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "DBH GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "DBH obs", low = "light grey", high = "black") +
labs(title = paste("DBH (cm)\nR2 =",signif(surf.dbh.anova$r.sq, 2),", P =", signif(surf.dbh.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#neo_fruiting
surf.neo_fruiting.grid = expand.grid(x = surf.neo_fruiting$grid$x, y = surf.neo_fruiting$grid$y)
surf.neo_fruiting.grid$z = as.vector(surf.neo_fruiting$grid$z)
surf.neo_fruiting.grid = data.frame(na.omit(surf.neo_fruiting.grid))

p.neo_fruiting = ggplot() +
geom_point(data = asv_tab.gt5K.rare.mds.bin.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = NeoFruiting), size = 2.5) +
stat_contour(data = surf.neo_fruiting.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "neo_fruiting GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "neo_fruiting obs", low = "light grey", high = "black") +
labs(title = paste("Neonectria fruiting (0-5)\nR2 =",signif(surf.neo_fruiting.anova$r.sq, 2),", P =", signif(surf.neo_fruiting.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#cankers
surf.cankers.grid = expand.grid(x = surf.cankers$grid$x, y = surf.cankers$grid$y)
surf.cankers.grid$z = as.vector(surf.cankers$grid$z)
surf.cankers.grid = data.frame(na.omit(surf.cankers.grid))

p.cankers = ggplot() +
geom_point(data = asv_tab.gt5K.rare.mds.bin.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = RaisedCanker), size = 2.5) +
stat_contour(data = surf.cankers.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "cankers GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "cankers obs", low = "light grey", high = "black") +
labs(title = paste("Cankers\nR2 =",signif(surf.cankers.anova$r.sq, 2),", P =", signif(surf.cankers.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#tarry_spots
surf.tarry_spots.grid = expand.grid(x = surf.tarry_spots$grid$x, y = surf.tarry_spots$grid$y)
surf.tarry_spots.grid$z = as.vector(surf.tarry_spots$grid$z)
surf.tarry_spots.grid = data.frame(na.omit(surf.tarry_spots.grid))

p.tarry_spots = ggplot() +
geom_point(data = asv_tab.gt5K.rare.mds.bin.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = TarrySpots), size = 2.5) +
stat_contour(data = surf.tarry_spots.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "tarry_spots GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "tarry_spots obs", low = "light grey", high = "black") +
labs(title = paste("Tarry spots\nR2 =",signif(surf.tarry_spots.anova$r.sq, 2),", P =", signif(surf.tarry_spots.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#wax
surf.wax.grid = expand.grid(x = surf.wax$grid$x, y = surf.wax$grid$y)
surf.wax.grid$z = as.vector(surf.wax$grid$z)
surf.wax.grid = data.frame(na.omit(surf.wax.grid))

p.wax = ggplot() +
geom_point(data = asv_tab.gt5K.rare.mds.bin.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = Wax), size = 2.5) +
stat_contour(data = surf.wax.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "wax GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "wax obs", low = "light grey", high = "black") +
labs(title = paste("Wax (0-5)\nR2 =",signif(surf.wax.anova$r.sq, 2),", P =", signif(surf.wax.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#xylococcus
surf.xylococcus.grid = expand.grid(x = surf.xylococcus$grid$x, y = surf.xylococcus$grid$y)
surf.xylococcus.grid$z = as.vector(surf.xylococcus$grid$z)
surf.xylococcus.grid = data.frame(na.omit(surf.xylococcus.grid))

p.xylococcus = ggplot() +
geom_point(data = asv_tab.gt5K.rare.mds.bin.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = Xylococcus), size = 2.5) +
stat_contour(data = surf.xylococcus.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "xylococcus GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "xylococcus obs", low = "light grey", high = "black") +
labs(title = paste("Xylococcus\nR2 =",signif(surf.xylococcus.anova$r.sq, 2),", P =", signif(surf.xylococcus.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#tree_cond
surf.tree_cond.grid = expand.grid(x = surf.tree_cond$grid$x, y = surf.tree_cond$grid$y)
surf.tree_cond.grid$z = as.vector(surf.tree_cond$grid$z)
surf.tree_cond.grid = data.frame(na.omit(surf.tree_cond.grid))

p.tree_cond = ggplot() +
geom_point(data = asv_tab.gt5K.rare.mds.bin.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = TreeCond), size = 2.5) +
stat_contour(data = surf.tree_cond.grid, aes(x = x, y = y , z = z, colour = (..level..)), size = 1) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both"), guide = F) +
scale_color_gradient(name = "tree_cond GAM", low = "light grey", high = "black") +
scale_fill_gradient(name = "tree_cond obs", low = "light grey", high = "black") +
labs(title = paste("Tree condition (0-5)\nR2 =",signif(surf.tree_cond.anova$r.sq, 2),", P =", signif(surf.tree_cond.anova$s.pv,2) ), x = "NMDS1", y = "NMDS2")

#Site and neo occurence for reference

p.site = ggplot() +
geom_point(data = asv_tab.gt5K.rare.mds.bin.metadata.neoOcurence, aes(MDS1, MDS2, shape = occurence, fill = Site), size = 3) +
my_gg_theme +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
scale_shape_manual(name = "Neonectria\noccurence", values = c(21,22,23,24), labels = c("Nd" = "N. ditissima", "Nf" = "N. faginata", "none" = "none", "both" = "both")) +
scale_color_gradient(name = "ppt", low = "light grey", high = "black") +
scale_fill_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")) +
labs(title = "Sites and\nneonectria occurence", x = "NMDS1", y = "NMDS2") +
guides(fill = guide_legend(override.aes=list(shape=21,size=3)), shape = guide_legend(override.aes=list(size=3)))


pdf("GAM_fits/NMDS.bin_neo_occurence_by_disease_severity_5Kmin_BINARY_GAM.sum_trees.tree_level_cov.pdf", width = 22, height = 18)
grid.arrange(p.site,p.dbh,p.neo_fruiting,p.wax,p.cankers,p.tarry_spots,p.xylococcus,p.tree_cond,nrow = 3)
dev.off()

