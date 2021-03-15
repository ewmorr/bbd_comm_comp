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

#and calcuduration_infectiones neonectria occurence in objects:
#Nf_v_Nd.long.metadata
#Nf_v_Nd.bin.metadata

#######
#PLOTS#

######
#NMDS#

##########################
#1000 seqs per sample min#

#rarefied table
asv_tab.gt1K.rare = readRDS(file = "intermediate_RDS/asv_tab.gt1K.rare.tree_sum.rds")
asv_tab.gt1K.rare = asv_tab.gt1K.rare[,colSums(asv_tab.gt1K.rare) > 0]
#NMDS
asv_tab.gt1K.rare.mds = readRDS(file = "intermediate_RDS/asv_tab.gt1K.rare.mds.tree_sum.rds")

#filter metadata
full_metadata.1K = full_metadata %>% filter(sample %in% rownames(asv_tab.gt1K.rare))
#Add metadata

asv_tab.gt1K.rare.mds.metadata = left_join(
data.frame(sample = rownames(asv_tab.gt1K.rare.mds$points), asv_tab.gt1K.rare.mds$points),
full_metadata.1K, by = "sample")

#PLOT
#climate/continental scale vars

p1 = ggplot(asv_tab.gt1K.rare.mds.metadata, aes(MDS1, MDS2, fill = as.factor(n_plugs))) +
geom_point(size = 2.5, shape = 22, color = "black") +
my_gg_theme +
scale_fill_brewer(palette = "Dark2") +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
labs(fill = "no. plugs\nper tree")

pdf("n_plug_effects/full_NMDS_n_plugs.1K.pdf", width = 8, height = 6)
print(p1)
dev.off()

n_plugs.adonis.1K = adonis(log10(asv_tab.gt1K.rare+1) ~ full_metadata.1K$Site * full_metadata.1K$n_plugs)
n_plugs.adonis.bin.1K = adonis(asv_tab.gt1K.rare ~ full_metadata.1K$n_plugs, binary = T)

#####################
#Richness by duration_infection etc#

sample_rarefied_richness = rarefy(t(asv_tab[,colSums(asv_tab) > 1000]), sample = 1000)
sample_rarefied_richness.metadata = left_join(
    data.frame(sample = names(sample_rarefied_richness), richness = sample_rarefied_richness),
    full_metadata,
    by = "sample"
)

summary(aov(richness ~ n_plugs*Site, data = sample_rarefied_richness.metadata))

anova(lm(richness ~ n_plugs*Site, data = sample_rarefied_richness.metadata))

p1 = ggplot(sample_rarefied_richness.metadata, aes(as.factor(n_plugs), richness)) +
geom_boxplot(outlier.size = 3) +
geom_point(position = position_jitter(), alpha = 0.5) +
my_gg_theme +
labs(x = "No. plugs per tree", y = "rarefied richness (1K seqs)")

sample_rarefied_richness = rarefy(t(asv_tab[,colSums(asv_tab) > 5000]), sample = seq(100, 5000, 100))
sample_rarefied_richness = as.data.frame(sample_rarefied_richness)
sample_rarefied_richness$sample = asv_tab[,colSums(asv_tab) > 5000] %>% colnames

sample_rarefied_richness.long = sample_rarefied_richness %>% pivot_longer(cols = -sample, names_to = "subsample", values_to = "richness")
sample_rarefied_richness.long$subsample = sub("N", "", sample_rarefied_richness.long$subsample) %>% as.numeric

sample_rarefied_richness.long = left_join(sample_rarefied_richness.long, full_metadata %>% select(sample, n_plugs), by = "sample")

p2= ggplot(sample_rarefied_richness.long, aes(subsample, richness, color = as.factor(n_plugs))) +
geom_line(aes(group = sample), alpha = 0.75) +
scale_color_brewer(palette = "Dark2") +
labs(color = "No. plugs\nper tree", y = "rarefied richness") +
my_gg_theme +
theme(legend.title = element_text(size = 18))

pdf("n_plug_effects/rarefied_richness_n_plugs.pdf", width = 10, height = 10)
grid.arrange(p1,p2, ncol = 1)
dev.off()

##########################
#within site comm. cmomp.#

site_levels = full_metadata$Site %>% unique

nmds_list = vector("list", 10)

for(i in 1:10){
    row_names = full_metadata %>% filter(Site == site_levels[i]) %>% select(sample)
    temp_tab = asv_tab.gt1K.rare[rownames(asv_tab.gt1K.rare) %in% row_names$sample, ]
    temp_tab = temp_tab[,colSums(temp_tab) > 0]
    temp_nmds = metaMDS(log10(temp_tab +1), autotransform = F, k = 2)
    nmds_list[[i]] = temp_nmds
}

nmds_list.points = vector("list", 10)

for(i in 1:10){
    nmds_list.points[[i]] = nmds_list[[i]]$points
}

nmds_list.points.df = do.call("rbind", nmds_list.points)

nmds_list.points.df.metadata = left_join(
    data.frame(sample = rownames(nmds_list.points.df), nmds_list.points.df),
    full_metadata %>% select(sample, n_plugs, Site),
    by = "sample"
)

p1 = ggplot(nmds_list.points.df.metadata, aes(MDS1, MDS2,  fill = as.factor(n_plugs))) +
facet_wrap(~Site, scale = "free", ncol = 4) +
geom_point(size = 2.5, shape = 22, color = "black") +
my_gg_theme +
scale_fill_brewer(palette = "Dark2") +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
labs(fill = "no. plugs\nper tree")

pdf("n_plug_effects/NMDS_within_site_plugs.pdf", width = 10, height = 6)
print(p1)
dev.off()


p1 = ggplot(asv_tab.gt1K.rare.mds.metadata, aes(MDS1, MDS2, fill = as.factor(n_plugs))) +
geom_point(size = 2.5, shape = 22, color = "black") +
facet_wrap(~Site, scale = "free", ncol = 4) +
my_gg_theme +
scale_fill_brewer(palette = "Dark2") +
theme(legend.title = element_text(size = 17, hjust = 0), legend.text = element_text(size = 15)) +
labs(fill = "no. plugs\nper tree")
#The full NMDS provides similar qualitative groupings and is esier to see
pdf("n_plug_effects/NMDS_within_site_plugs.pdf", width = 10, height = 6)
print(p1)
dev.off()


#PERMANOVA

adonis_list = vector("list")

for(i in 1:10){
    row_names = full_metadata %>% filter(Site == site_levels[i]) %>% select(sample)
    temp_tab = asv_tab.gt1K.rare[rownames(asv_tab.gt1K.rare) %in% row_names$sample, ]
    temp_tab = temp_tab[,colSums(temp_tab) > 0]
    temp_meta = full_metadata %>% filter(sample %in% rownames(temp_tab)) %>% select(sample, n_plugs)
    print(max(temp_meta$n_plugs))
    if(max(temp_meta$n_plugs)>1){
         rownames(temp_meta) = temp_meta$sample
        temp_meta = temp_meta[match(rownames(temp_tab), rownames(temp_meta)),]
        temp_adonis = adonis(log10(temp_tab +1) ~ temp_meta$n_plugs)
        adonis_list[[site_levels[i]]] = temp_adonis
    }
}

adonis_tab = data.frame(
    Site = names(adonis_list),
    R2 = vector(mode = "numeric", length = length(adonis_list)),
    p.val = vector(mode = "numeric", length = length(adonis_list))
)

for(i in 1:length(adonis_list)){
    site_nam = names(adonis_list)[i]
    adonis_tab$Site[i] = names(adonis_list)[i]
    adonis_tab$R2[i] = adonis_list[[site_nam]]$aov.tab$R2[1]
    adonis_tab$p.val[i] = adonis_list[[site_nam]]$aov.tab$`Pr(>F)`[1]
}

write.table(adonis_tab, file = "n_plug_effects/within_site_adonis.txt", row.names = F, sep = "\t", quote = F)
