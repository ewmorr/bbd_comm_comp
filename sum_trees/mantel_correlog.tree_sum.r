require(tidyverse)
require(vegan)
#require(gridExtra)
#require(RColorBrewer)

source("~/ggplot_theme.txt")

#read data
source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r")


#rarefied table
asv_tab.gt1K.rare = readRDS(file = "intermediate_RDS/asv_tab.gt1K.rare.tree_sum.rds")

full_metadata = full_metadata %>% filter(sample %in% rownames(asv_tab.gt1K.rare))

site_coords = data.frame(full_metadata$lat, full_metadata$lon)

convert_lat_lon_to_cartesian = function(lat, lon){
    #https://stackoverflow.com/questions/1185408/converting-from-longitude-latitude-to-cartesian-coordinates
    lat_rad = lat*pi/180
    lon_rad = lon*pi/180
    
    R = 6371
    x = R * cos(lat_rad) * cos(lon_rad)
    y = R * cos(lat_rad) * sin(lon_rad)
    z = R *sin(lat_rad)
    return(c(x,y,z))
}

site_cart_coords = data.frame(x = vector("numeric", 10), y = vector("numeric", 10), z = vector("numeric", 10))

for(i in 1:length(site_coords$full_metadata.lat)){
    site_cart_coords[i,1:3] = convert_lat_lon_to_cartesian(site_coords[i,1], site_coords[i,2])
}

rownames(site_cart_coords) = rownames(asv_tab.gt1K.rare)

mantel.correlog(
    D.eco = vegdist(log10(asv_tab.gt1K.rare+1), method = "bray"),
    D.geo = vegdist(site_cart_coords, method = "euclidean"),
    n.class = 0#,
    #break.pts = c(200,400,600,800,1000,1200,1400,1600,1800)
) %>% plot

site.mantel.correlog = mantel.correlog(
D.eco = vegdist(log10(asv_tab.gt1K.rare+1), method = "bray"),
D.geo = vegdist(site_cart_coords, method = "euclidean"),
#mult = "hochberg",
n.class = 0#,
#break.pts = c(200,400,600,800,1000,1200,1400,1600,1800)
)

mantel(
xdis = vegdist(log10(asv_tab.gt1K.rare+1), method = "bray"),
ydis = vegdist(site_cart_coords, method = "euclidean")
)

require(ggplot2)

mantel.plot.df = data.frame(distance = site.mantel.correlog[[1]][,1], Mantel.cor = site.mantel.correlog[[1]][,3], P = site.mantel.correlog[[1]][,5]
)

#mantel.plot.df[is.na(mantel.plot.df)] = 0

for(i in 1:length(mantel.plot.df$P)){
    if(is.na(mantel.plot.df$P[i])){next}
    if(mantel.plot.df$P[i] < 0.05 & mantel.plot.df$Mantel.cor[i] > 0){
        mantel.plot.df$P.fac[i] = "P < 0.05"
    }else{
        mantel.plot.df$P.fac[i] = "P >= 0.05"
    }
}

p1 = ggplot(mantel.plot.df, aes(distance, Mantel.cor)) +
geom_line() +
geom_point(aes(fill = P.fac), size = 4, shape = 21) +
scale_fill_manual(values = c("P < 0.05" = "black", "P >= 0.05" = "grey")) +
scale_x_continuous(breaks = c(0,200,400,600,800,1000,1200,1400,1600), limits = c(0,1700)) +
my_gg_theme +
labs(x = "Distance (km)", y = "Mantel r")
#scale_color_manual

pdf("prelim_figs/mantel_correlogram.trees.pdf", width = 8, height = 4)
print(p1)
dev.off()

