require(tidyverse)
require(data.table)
get.sample.name <- function(fname) strsplit(basename(fname), "_[ATCG]{8}(-|\\.)[ATCG]{8}", perl = T)[[1]][1]

#read data
asv_tax = read.table("~/GARNAS_neonectria_barcoding_files_cat_03242020/dada2_out/ASVs_taxonomy.tsv", header = T)
asv_tab = read.table("~/GARNAS_neonectria_barcoding_files_cat_03242020/dada2_out/ASVs_counts.tsv", sep = "\t", header = T, row.names = 1)
track.long = read.csv("~/GARNAS_neonectria_barcoding_files_cat_03242020/dada2_processing_tables_figs/read_processing_tracking.csv", row.names = 1)
id_bench_map = read.table("~/GARNAS_neonectria_barcoding_files_cat_03242020/sample_data/sample_mapping.txt", header = T)
metadata_map = read.table("~/GARNAS_neonectria_barcoding_files_cat_03242020/sample_data/metadata.txt", header = T)
survey_dat = read.table("~/GARNAS_neonectria_barcoding_files_cat_03242020/sample_data/trees_site_survey_data.txt", header = T, sep = "\t")
neo_cov = read.table("~/GARNAS_neonectria_barcoding_files_cat_03242020/sample_data/plug_neonectria_coverage.txt", header = T)
site_info = read.csv("~/GARNAS_neonectria_barcoding_files_cat_03242020/sample_data/site_info.csv")
site_means = read.table("~/GARNAS_neonectria_barcoding_files_cat_03242020/sample_data/BBD_survey_transect_data.site_mean.txt", header = T)
site_climate = read.table("~/GARNAS_neonectria_barcoding_files_cat_03242020/sample_data/sites_climate.txt", header = T)
site_climate.GDD = read.table(file = "~/GARNAS_neonectria_barcoding_files_cat_03242020/sample_data/site_info.GDD.txt", header = T)
#Some gymnastics to get sample labels in the correct format. Could fix this upstream...
colnames(asv_tab) = unname(sapply(colnames(asv_tab), get.sample.name))
track.long$sample <- unname(sapply(as.character(track.long$sample), get.sample.name))
#track.long$sample = paste0("X", track.long$sample)

#table joins
metadata_ordered = full_join(metadata_map, id_bench_map)
survey_dat.neo_cov = full_join(survey_dat, neo_cov, by = c("Site", "Tree", "Plug")) %>%
    left_join(., site_info, by = "Site") %>%
    left_join(., site_means, by = "Site") %>%
    left_join(., site_climate, by = c("Site","lat","lon")) %>%
    left_join(., site_climate.GDD, by = "Site")

full_metadata = full_join(metadata_ordered, survey_dat.neo_cov, by = c("Site", "Tree", "Plug"))
if("seq.rep" %in% colnames(full_metadata)){
    full_metadata = full_metadata %>% filter(seq.rep != "y")
}

#########################################
#Filter plants and animals out of tables#

plant_asvs = filter(data.frame(asv_tax, asv_name = rownames(asv_tax)), Kingdom == "k__Viridiplantae")$asv_name
animal_asvs = filter(data.frame(asv_tax, asv_name = rownames(asv_tax)), Kingdom == "k__Metazoa")$asv_name

asv_tab = subset(asv_tab, !rownames(asv_tab) %in% plant_asvs)
asv_tab = subset(asv_tab, !rownames(asv_tab) %in% animal_asvs)

asv_tax = subset(asv_tax, !rownames(asv_tax) %in% plant_asvs)
asv_tax = subset(asv_tax, !rownames(asv_tax) %in% animal_asvs)

################################################
#Negative & control samples table with taxonomy#

#make negatives only asv_tab (long format)
asv_tab.negatives = semi_join(
    data.frame(sample = rownames(t(asv_tab)), t(asv_tab)),
    full_metadata %>% filter(bench.control != "n")
)
rownames(asv_tab.negatives) = asv_tab.negatives$sample
asv_tab.negatives$sample = NULL
asv_tab.negatives = t(asv_tab.negatives)
#join with taxonomy
asv_tab.negatives.asvnames = data.frame(ASV = rownames(asv_tab.negatives), asv_tab.negatives)

asv_tab.negatives.long = melt(asv_tab.negatives.asvnames[rowSums(asv_tab.negatives) > 0,] %>% data.table,
id = "ASV", variable.name = "sample", value.name = "count") %>%
data.frame


#################################################
#Process asv_tax for lowest informative taxonomy#

asv_tax.char = apply(asv_tax, 2, as.character)
asv_tax.char[is.na(asv_tax.char)] = "unknown"
rownames(asv_tax.char) = rownames(asv_tax)

get_informative_taxa = function(x){
    found_info = 0
    for(i in length(x):2){
        if(x[i] != "unknown"){
            if(i == length(x)){
                return(paste(as.character(x[i-1]), sub("s__", "", as.character(x[i]) ), sep = " " ))
            }
            else{
                return(as.character(x[i]))
            }
            found_info = 1
            break
        }
    }
    if(found_info == 0){return(x[1])}
}

asv_informative_taxa = vector(mode = "character", length = length(rownames(asv_tax.char)))
asv_informative_taxa = apply(asv_tax.char, 1, function(x) get_informative_taxa(x))

#############################
#Nf and Nd counts (sum ASVs)#

#Also read mapping data in order to avg mapping counts and ASV counts for detection
#read the mapping data
asv_tab.mapping = read.table("~/GARNAS_neonectria_barcoding_files_cat_03242020/neo_map/all_seqs.derep.otu_tab.txt", header = T)
rownames(asv_tab.mapping) = asv_tab.mapping$OTUID
asv_tab.mapping$OTUID = NULL
#remove negative controls from both tables
asv_tab = asv_tab[,!colnames(asv_tab) %in% colnames(asv_tab.negatives)]
asv_tab.mapping = asv_tab.mapping[,!colnames(asv_tab.mapping) %in% colnames(asv_tab.negatives)]
#BP225 did not make it trhough DADA2 so has to be removed manually
asv_tab.mapping$BP225 = NULL
#############################
#Nf and Nd counts (sum ASVs)#

Nf_asvs = filter(data.frame(asv_tax, asv_name = rownames(asv_tax)), Species == "s__faginata")$asv_name
Nd_asvs = filter(data.frame(asv_tax, asv_name = rownames(asv_tax)), Species == "s__ditissima")$asv_name

Nf_counts = subset(asv_tab, rownames(asv_tab) %in% Nf_asvs) %>% colSums
Nd_counts = subset(asv_tab, rownames(asv_tab) %in% Nd_asvs) %>% colSums

Nf_v_Nd = full_join(
data.frame(Nf = Nf_counts, sample = names(Nf_counts)),
data.frame(Nd = Nd_counts, sample = names(Nd_counts)),
by = "sample")

Nf_v_Nd.long = gather(Nf_v_Nd, "spp", "Neonectria_count.asv", -sample)

#also do mapping file

Nf_counts.mapping = subset(asv_tab.mapping, rownames(asv_tab.mapping) %in% Nf_asvs) %>% colSums
Nd_counts.mapping = subset(asv_tab.mapping, rownames(asv_tab.mapping) %in% Nd_asvs) %>% colSums

Nf_v_Nd.mapping = full_join(
data.frame(Nf = Nf_counts.mapping, sample = names(Nf_counts.mapping)),
data.frame(Nd = Nd_counts.mapping, sample = names(Nd_counts.mapping)),
by = "sample")

Nf_v_Nd.mapping.long = gather(Nf_v_Nd.mapping, "spp", "Neonectria_count.mapping", -sample)

#join mapping and ASV counts
Nf_Nd.asv_mapping_comp = full_join(Nf_v_Nd.long, Nf_v_Nd.mapping.long, by = c("sample", "spp"))
#vsearch does not write files with zero counts, so add zero for NAs
Nf_Nd.asv_mapping_comp[is.na(Nf_Nd.asv_mapping_comp)] = 0
#add sample sequence totals from asv_tab
Nf_Nd.asv_mapping_comp.w_counts = left_join(Nf_Nd.asv_mapping_comp,
data.frame(sample = colnames(asv_tab), total_seqs = colSums(asv_tab))
)

#CAL AVG
Nf_Nd.asv_mapping_comp$Neo_avg = (Nf_Nd.asv_mapping_comp$Neonectria_count.asv + Nf_Nd.asv_mapping_comp$Neonectria_count.mapping) / 2

#remake count vectors
Nf_counts = as.vector((filter(Nf_Nd.asv_mapping_comp, spp == "Nf") %>% dplyr::select(Neo_avg))[[1]] )
names(Nf_counts) = (filter(Nf_Nd.asv_mapping_comp, spp == "Nf") %>% dplyr::select(sample))[[1]]

Nd_counts = as.vector((filter(Nf_Nd.asv_mapping_comp, spp == "Nd") %>% dplyr::select(Neo_avg))[[1]] )
names(Nd_counts) = (filter(Nf_Nd.asv_mapping_comp, spp == "Nd") %>% dplyr::select(sample))[[1]]

#############
#SUM by tree#
#############

site_tree_label = full_metadata %>% filter(bench.control == "n" & locus == "ITS2") %>% dplyr::select(c("Site", "Tree", "metadata.label"))
site_tree_label$Site.tree = interaction(site_tree_label$Site, site_tree_label$Tree)

#########################
#Make sure colnames match

b = names(Nf_counts)
a = site_tree_label$metadata.label

not_in_tab = unique(a[! a %in% b])
site_tree_label = site_tree_label %>% filter(!metadata.label %in% not_in_tab)

#New vectors indexed by tree instead of metadata.label
Nf_counts.tree_sum = vector(mode = "numeric", length = site_tree_label$Site.tree %>% unique %>% length)
names(Nf_counts.tree_sum) = site_tree_label$Site.tree %>% unique

Nd_counts.tree_sum = vector(mode = "numeric", length = site_tree_label$Site.tree %>% unique %>% length)
names(Nd_counts.tree_sum) = site_tree_label$Site.tree %>% unique

#add vals
for(i in 1:length(site_tree_label$Site.tree)){
    #print(paste(site_tree_label$Site.tree[i], "plus", site_tree_label$metadata.label[i]))
    
    Nf_counts.tree_sum[names(Nf_counts.tree_sum) == site_tree_label$Site.tree[i]] = Nf_counts.tree_sum[names(Nf_counts.tree_sum) == site_tree_label$Site.tree[i]] + Nf_counts[names(Nf_counts) == site_tree_label$metadata.label[i]]
    
    Nd_counts.tree_sum[names(Nd_counts.tree_sum) == site_tree_label$Site.tree[i]] = Nd_counts.tree_sum[names(Nd_counts.tree_sum) == site_tree_label$Site.tree[i]] + Nd_counts[names(Nd_counts) == site_tree_label$metadata.label[i]]
}

#check that sums are equal
#sum(Nf_counts.tree_sum) - sum(Nf_counts)
#sum(Nd_counts.tree_sum) - sum(Nd_counts)

Nf_counts = Nf_counts.tree_sum
Nd_counts = Nd_counts.tree_sum

Nf_v_Nd = full_join(
data.frame(Nf = Nf_counts, sample = names(Nf_counts)),
data.frame(Nd = Nd_counts, sample = names(Nd_counts)),
by = "sample")

Nf_v_Nd.long = gather(Nf_v_Nd, "spp", "Neonectria_count.asv", -sample)

###########################
#Nf and Nd occurence (0/1)#

Nf_v_Nd.bin = full_join(
data.frame(Nf = as.numeric(as.matrix(Nf_counts) > 0), sample = names(Nf_counts)),
data.frame(Nd = as.numeric(as.matrix(Nd_counts) > 0), sample = names(Nd_counts)),
by = "sample")

for(i in 1:length(Nf_v_Nd.bin$sample)){
    if(Nf_v_Nd.bin$Nf[i] == 1 & Nf_v_Nd.bin$Nd[i] == 1){
        Nf_v_Nd.bin$occurence[i] = "both"
    }
    if(Nf_v_Nd.bin$Nf[i] == 1 & Nf_v_Nd.bin$Nd[i] == 0){
        Nf_v_Nd.bin$occurence[i] = "Nf"
    }
    if(Nf_v_Nd.bin$Nf[i] == 0 & Nf_v_Nd.bin$Nd[i] == 1){
        Nf_v_Nd.bin$occurence[i] = "Nd"
    }
    if(Nf_v_Nd.bin$Nf[i] == 0 & Nf_v_Nd.bin$Nd[i] == 0){
        Nf_v_Nd.bin$occurence[i] = "none"
    }
}

##############
#Add metadata#

####################################################################################
#First transform full_metadata to contain tree level values and sample as Site.tree#

full_metadata.site_tree = full_metadata %>% filter(bench.control == "n" & locus == "ITS2" & !metadata.label %in% not_in_tab)
full_metadata.site_tree$sample = paste(as.character(full_metadata.site_tree$Site), as.character(full_metadata.site_tree$Tree), sep = ".")

#add number of plugs sampled per tree
plugs_per_sample = full_metadata.site_tree %>% group_by(sample) %>% dplyr::select(sample, Plug) %>% summarize(n_plugs = n())
full_metadata.site_tree = full_metadata.site_tree %>% dplyr::select(-Plug, -metadata.label, -bench.control, -locus, -Pool, -Site_tree_plug, -NeoCoverage, -cambium, -phloem.depth)
full_metadata.site_tree = left_join(full_metadata.site_tree %>% distinct, plugs_per_sample, by = "sample")

#then add metadata to dfs
Nf_v_Nd.long.metadata = left_join(Nf_v_Nd.long, full_metadata.site_tree, by = "sample")

Nf_v_Nd.bin.metadata = left_join(Nf_v_Nd.bin, full_metadata.site_tree, by = "sample")

##########################################################
#Read in LULU filtered (93% minimum_simlarity) ASV table#

asv_tab = read.table("~/GARNAS_neonectria_barcoding_files_cat_03242020/LULU/asv_tab.LULU_93.txt", header = T)
#remove negatives
asv_tab = asv_tab[,!colnames(asv_tab) %in% colnames(asv_tab.negatives)]
#remove zero sum asvs (though there actually are none)
asv_tab = asv_tab[rowSums(asv_tab) > 1,]

########################
#SUM ASV counts by tree#
########################

site_tree_label = full_metadata %>% filter(bench.control == "n" & locus == "ITS2") %>% dplyr::select(c("Site", "Tree", "metadata.label"))

site_tree_label$Site.tree = interaction(site_tree_label$Site, site_tree_label$Tree)

#########################
#Make sure colnames match

b = colnames(asv_tab)
a = site_tree_label$metadata.label

not_in_tab = unique(a[! a %in% b])
site_tree_label = site_tree_label %>% filter(!metadata.label %in% not_in_tab)

####################################################
#new table indexed by tree instead of metadata.label
asv_tab.tree_sum = data.frame(matrix(
data = 0,
nrow = length(rownames(asv_tab)),
ncol = site_tree_label$Site.tree %>% unique %>% length
))

rownames(asv_tab.tree_sum) = rownames(asv_tab)
colnames(asv_tab.tree_sum) = site_tree_label$Site.tree %>% unique

#add cols
for(i in 1:length(site_tree_label$Site.tree)){
    #tree.col = asv_tab.tree_sum[,colnames(asv_tab.tree_sum) == site_tree_label$Site.tree[i]]
    print(paste(site_tree_label$Site.tree[i], "plus", site_tree_label$metadata.label[i]))
    asv_tab.tree_sum[,colnames(asv_tab.tree_sum) == site_tree_label$Site.tree[i]] = asv_tab.tree_sum[,colnames(asv_tab.tree_sum) == site_tree_label$Site.tree[i]] + asv_tab[,colnames(asv_tab) == site_tree_label$metadata.label[i]]
}

asv_tab = asv_tab.tree_sum

#filter out taxa that only occur in one sample
asv_tab = asv_tab[rowSums(asv_tab > 0) > 1,]
#Also filter asv_tax based on retained ASVs
asv_tax = asv_tax[rownames(asv_tax) %in% rownames(asv_tab),]

############################################
#replace original matdata w/ transformed df#

full_metadata = full_metadata.site_tree
total_seqs = data.frame(sample = colnames(asv_tab), total_seqs = colSums(asv_tab))
full_metadata = full_join(full_metadata, total_seqs, by = "sample")
Nf_v_Nd.long.metadata = full_join(Nf_v_Nd.long.metadata, total_seqs, by = "sample")
Nf_v_Nd.bin.metadata = full_join(Nf_v_Nd.bin.metadata, total_seqs, by = "sample")
