
source("~/repo/neonectria_barcoding_012220/sum_trees/read_ASV_dat.LULU_tab.r") #when we use the version in the bbd_comm_comp repo the script breaks for some reason. Or it might be loqding tidyverse in this script instead of through the sourced script... Will have to test

library(Hmsc)
library(parallel)
library(corrplot)
library(MASS)
set.seed(6)

###############################
#Y matrix

#remove samples with less than 1000 seqs
lt_1K_samps = full_metadata %>%
    filter(total_seqs < 1000) %>%
    dplyr::select("sample")

#Species matrix filtered to only those samples with >= 1000 sequences and THEN species that occur in 10% of samples
Ytable = t(asv_tab[,!colnames(asv_tab) %in% lt_1K_samps$sample])
yprob = as.matrix(Ytable)
yprob[yprob > 0] = 1
yprob = yprob[,colSums(yprob) >= nrow(yprob)*0.1]

######################
#SHOULD also run this with rarefied data and no total_seqs correction to look for effects
######################

#filter and order metadata
full_metadata.sorted = left_join(
data.frame(sample = rownames(Ytable)),
full_metadata,
by = "sample"
) %>% left_join(
., Nf_v_Nd.bin
)


###Remove NF and ND asvs and replace with mapping based occurence
###Going to skip this step, because we are interested in comm comp at large and other species may have multiple ASVs
#yprob = yprob[,!colnames(yprob) %in% Nf_asvs]
#yprob = yprob[,!colnames(yprob) %in% Nd_asvs]

#yprob = cbind(yprob, Nf = full_metadata.sorted$Nf)
#yprob = cbind(yprob, Nd = full_metadata.sorted$Nd)

#####################
#Define random levels
xycoords = matrix(c(site_info$lat, site_info$lon),ncol=2)

colnames(xycoords) = c("x-coordinate","y-coordinate")
rownames(xycoords) = site_info$Site

studyDesign = data.frame(sample = as.factor(full_metadata.sorted$sample), site_spatial = as.factor(full_metadata.sorted$Site))

rL.sample = HmscRandomLevel(units = as.factor(full_metadata.sorted$sample))
rL.site_spatial = HmscRandomLevel(sData = xycoords)

#################
#covariate matrix

#For now including ALL potential climate variables. This is fine because we deal with covariance in the grouping arguments to variance partitioning. NEED TO THINK CAREFULLY about which disease metrics to include. For example, it might not make sense to include infection duration here because this is not a BBD specific analysis where DI has some hypothesized value. Instead include actual measurements of disease state and tree health. What about conflicting measurements? For example, DBH could potentially have and interaction with TreeCond and Neonectria fruiting, where large DBH with high fruting have high tree cond, but large DBH with low fruting have low TreeCond. Probably would mean that DBH should be excluded from the disease state group (but still in the model?) OR an interaction term should be included, but which is the right choice??? Similar potential problems with wax density
#On rethinking the tree level vars prblem, it shouldn't matter whether there is a postive or negative effect and/or interaction, we are looking for total variance explained
#leaving out NeoFruting for now because it is a symptom of one of the particular ASVs (and we know it explains a lot of variance in Nf)
XData = full_metadata.sorted %>% dplyr::select("HDD4.mean_nongrowing", "freezeThaw.mean_nongrowing", "ppt.mean_nongrowing", "HDD4.mean_growing", "freezeThaw.mean_growing","ppt.mean_growing", "RaisedCanker", "dbh", "Wax", "TreeCond", "total_seqs")
XData[,colnames(XData) == "total_seqs"] = log(XData[,colnames(XData) == "total_seqs"])

#Scaling covariates
XData = apply(XData,2,scale)
XData = data.frame(XData)

m.spatial.full_covars = Hmsc(Y=yprob, XData=XData, XFormula=~.,
studyDesign=studyDesign, ranLevels=list("sample"=rL.sample, "site_spatial" = rL.site_spatial),distr="probit")

#set priors for total_seqs so that it is used as an offset for sampling effort
#This is failing with error saying V0 is not positive definite...
#K = ncol(XData)
#m.spatial <- setPriors(m.spatial,
#V0 = diag(c(rep(1, K-1), 0.00001)),
#f0 = K+1,
#mGamma = c(rep(0, K-1), 1),
#UGamma = diag(c(rep(1, K-1), 1e-8))
#)

#Set MCMC sampling paramteres and run
nChains = 2
test.run = F
if (test.run){
    # with this option, the vignette runs fast but results are not reliable
    thin = 1
    samples = 10
    transient = 5
    verbose = 0
} else {
    # with this option, the vignette evaluates slow but it reproduces the results of
    # the .pdf version
    thin = 10
    samples = 1000
    transient = 1000
    verbose = 0
}


m.spatial.full_covars = sampleMcmc(m.spatial.full_covars, thin = thin, samples = samples, transient = transient,
nChains = nChains, verbose = verbose,updater=list(GammaEta=FALSE))
saveRDS(m.spatial.full_covars, "intermediate_RDS/full_community_ge_1K_seqs_per_sample_gt_10_freq.BIN.HSMCS_MCMC.rds")
######################
#MCMC convergence

#ess == effective sample size
#psrf == potential scale reduction factor (closeto 1 indicates good MCMC convergence)
#beta == species niches (as realted to covarites)
#gamma == influence of traits on species niches
#omega == residual species association
#rho == phylogenetic signal

mpost = convertToCodaObject(m.spatial.full_covars)
par(mfrow=c(3,2))
ess.beta = effectiveSize(mpost$Beta)
psrf.beta = gelman.diag(mpost$Beta, multivariate=FALSE)$psrf
hist(ess.beta)
hist(psrf.beta)
ess.gamma = effectiveSize(mpost$Gamma)
psrf.gamma = gelman.diag(mpost$Gamma, multivariate=FALSE)$psrf
hist(ess.gamma)
hist(psrf.gamma)
ns = 50
sppairs = matrix(sample(x = 1:ns^2, size = 100))
tmp = mpost$Omega[[1]]
for (chain in 1:length(tmp)){
    tmp[[chain]] = tmp[[chain]][,sppairs]
}
ess.omega = effectiveSize(tmp)
psrf.omega = gelman.diag(tmp, multivariate=FALSE)$psrf
hist(ess.omega)
hist(psrf.omega)

##########################
#Model fit and partioning

preds = computePredictedValues(m.spatial.full_covars)
MF = evaluateModelFit(hM=m.spatial.full_covars, predY=preds)
hist(MF$TjurR2, xlim = c(0,1), main=paste0("Mean = ", round(mean(MF$TjurR2),2)))
#Above is R2 for all species and the mean of species. Important

#SHOULD MAKE A PLOT OF THIS

#######################
#Variance partioning
#the group variable assigns X covariates to different covariate groups
#so first look at the design matrix

head(m.spatial.full_covars$X)
#For our real data we fit an intercept and three continuous variables, so they can each be assigned separate groups
#If instead we had a categorical variable the levels could be assigned to a single group along with the intercept

VP = computeVariancePartitioning(m.spatial.full_covars, group = c(1, rep(2,6), rep(3,4), 4), groupnames = c("intercept", "site climate", "tree state", "total sequences"))
plotVariancePartitioning(m.spatial.full_covars, VP = VP)

saveRDS(VP, "intermediate_RDS/full_community_ge_1K_seqs_per_sample_gt_10_freq.BIN.VP.rds")

####################
#Plot variance partioning
#If included a square of continuous variable, negative response would indicate intermediate niche optimum (i.e., abundance goes up initially but then goes down)
postBeta = getPostEstimate(m.spatial.full_covars, parName = "Beta")
plotBeta(m.spatial.full_covars, post = postBeta, param = "Support",
plotTree = F, supportLevel = 0.95, split=.4, spNamesNumbers = c(T,F))
#This can also be mapped on a tree with plotTree = T, but then tree must be included in model

############################################
#transform VP and postBeta object for ggplot
#VP cpntains R2 vals and postBeta contains support (i.e. alpha)

VP.vals = data.frame(VP$vals)

VP.vals$variable = c("Intercept","Site climate", "Tree state", "Total sequences", "Tree-level\nrandom effect", "Site-level spatial\nrandom effect")
VP.vals.long = VP.vals %>%
    pivot_longer(-variable, names_to = "ASV", values_to = "R2")


#Transform R2 based on positive or negative response
postBeta.mean = data.frame(postBeta$mean)
#colnames(postBeta.mean) = c("Nf", "Nd")
postBeta.mean$variable = c("intercept","GDD nongrowing", "freeze-thaw nongrowing", "precip. nongrowing", "GDD growing", "freeze-thaw growing", "precip. growing", "cankers", "DBH", "beech scale", "crown dieback", "total sequences")
postBeta.mean.long = postBeta.mean %>%
pivot_longer(-variable, names_to = "ASV", values_to = "mean")



postBeta.support = data.frame(postBeta$support)
#colnames(postBeta.support) = c("Nf", "Nd")
postBeta.support$variable = c("intercept","GDD nongrowing", "freeze-thaw nongrowing", "precip. nongrowing", "GDD growing", "freeze-thaw growing", "precip. growing", "cankers", "DBH", "beech scale", "crown dieback", "total sequences")
postBeta.support.long = postBeta.support %>%
pivot_longer(-variable, names_to = "ASV", values_to = "support")

#It's not really necessary to extract both postive and negative support values because they are just inverse of each other. If support < 0.05 the negative slope is significant and vice versa
postBeta.supportNeg = data.frame(postBeta$supportNeg)
#colnames(postBeta.supportNeg) = c("Nf", "Nd")
postBeta.supportNeg$variable = c("intercept","GDD nongrowing", "freeze-thaw nongrowing", "precip. nongrowing", "GDD growing", "freeze-thaw growing", "precip. growing", "cankers", "DBH", "beech scale", "crown dieback", "total sequences")
postBeta.supportNeg.long = postBeta.supportNeg %>%
pivot_longer(-variable, names_to = "ASV", values_to = "supportNeg")


postBeta.mean.support = full_join(postBeta.support.long, postBeta.mean.long)

postBeta.mean.support$P.val = vector(mode = "character", length = length(postBeta.mean.support$support))

for(i in 1:length(postBeta.mean.support$P.val)){
    if(postBeta.mean.support$support[i] > 0.95 || postBeta.mean.support$support[i] < 0.05){
        postBeta.mean.support$P.val[i] = "P<0.05"
    }else{
        postBeta.mean.support$P.val[i] = "n.s."
    }
}

require(RColorBrewer)
source("~/ggplot_theme.txt")

postBeta.mean.support$variable = factor(postBeta.mean.support$variable, levels = c("intercept","GDD nongrowing", "freeze-thaw nongrowing", "precip. nongrowing", "GDD growing", "freeze-thaw growing", "precip. growing", "cankers", "DBH", "beech scale", "crown dieback", "total sequences"))

###############################
#Plots
#These are slightly different than the plots produced for the Neonectria paper becasue we have grouped the VP into similar types of variables (i.e., site level and tree level)

p1 = ggplot(postBeta.mean.support %>% filter(variable != "intercept" & P.val != "n.s."), aes(ASV, variable, fill = mean)) +
#geom_tile(size = 1, height = 0.975, width = 0.975) +
geom_tile(color = "black") +
scale_fill_gradient2(low = "#2c7bb6", high = "#d7191c", mid = "white", midpoint = 0, breaks = c(-1.5,0,1.5)) +
my_gg_theme +
labs(
    x = "",
    y = "",
    fill = "Slope"
) +
theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
   axis.text = element_text(size = 12),
    axis.text.x =element_blank()
)

pdf("HMSC_figs/full_comm_ge1K_seqs_per_sample_gt_33_freq.bin.spatial.sig_params.pdf", width = 10, height = 4)
p1
dev.off()


display.brewer.all(colorblindFriendly = T)

five_color_YlRed = c("#ffffb2", "#fecc5c", "#fd8d3c", "#f03b20", "#bd0026")
five_color_OrPR = c("#542788", "#fee0b6", "#998ec3", "#f1a340", "#d8daeb")

#Can sort by R2 values of specific group by first extracting that variable and then sorting the df
VP.vals.long$variable = factor(VP.vals.long$variable, levels = c("Intercept","Site climate", "Tree state", "Total sequences", "Tree-level\nrandom effect", "Site-level spatial\nrandom effect"))

VP.vals.long.site_sort = VP.vals.long %>% filter(variable == "Site climate")
VP.vals.long.site_sort = VP.vals.long.site_sort[order(VP.vals.long.site_sort$R2),]
VP.vals.long$ASV = factor(VP.vals.long$ASV, levels = VP.vals.long.site_sort$ASV)

VP.vals.vars_mean_R2 = VP.vals.long %>% filter(variable != "Intercept") %>% group_by(variable) %>% dplyr::summarize(mean = mean(R2))

p2 = ggplot(VP.vals.long %>% filter(variable != "Intercept"), aes(ASV, R2, fill = variable)) +
geom_bar(position="stack", stat="identity") +
#scale_fill_brewer(palette = "Paired")
#scale_fill_manual(values = five_color_YlRed)
scale_fill_manual(
values = rev(five_color_OrPR),
labels = c(paste(VP.vals.vars_mean_R2$variable, round(VP.vals.vars_mean_R2$mean, 2), sep = ": "))
) +
scale_y_continuous(expand = c(0.005,0.005)) +
scale_x_discrete(expand = c(0.005,0.005)) +
my_gg_theme +
#labs(y = expression(paste("Proportion of explained variance (R"^2, ")"))) +
labs(y = "Proportion of explained variance") +
theme(
axis.text.x = element_blank(),
legend.key.height=unit(2,"cm")
)

p2

pdf("HMSC_figs/full_comm_ge1K_seqs_per_sample_gt_33_freq.bin.spatial.VP.pdf", width = 12, height = 6)
p2
dev.off()

MF$TjurR2 %>% range
MF$TjurR2 %>% mean
MF$TjurR2 %>% sd

#scale_fill_continuous(limits = c(-1,1))

####################
####################
#Estimated residual var between spp

OmegaCor = computeAssociations(m.spatial.full_covars)
supportLevel = 0.95 #this is alpha = 0.05

#This is tree level (i.e., index [[1]]
toPlot = ((OmegaCor[[1]]$support>supportLevel)
+ (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
corrplot(toPlot, method = "color", order = "hclust", type = "upper", diag = T,
col=colorRampPalette(c("blue","white","red"))(500),
#tl.cex=.6, tl.col="black",
#tl.cex=0.1, tl.col="white",
tl.pos = "n",
bg = "black",
#title=paste("random effect level:", m.spatial.full_covars$rLNames[1]),
mar=c(2,0,1,0))

#This is site level (i.e., index [[2]]
toPlot = ((OmegaCor[[2]]$support>supportLevel)
+ (OmegaCor[[2]]$support<(1-supportLevel))>0)*OmegaCor[[2]]$mean
corrplot(toPlot, method = "color", order = "hclust", type = "lower", diag = F, add = T,
col=colorRampPalette(c("blue","white","red"))(500),
#tl.cex=0.01, tl.col="white",
tl.pos = "n",
bg = "black",
#title=paste("random effect level:", m.spatial.full_covars$rLNames[2]),
mar=c(2,0,1,0))

tree_level_support = OmegaCor[[1]]$support
tree_level_support[tree_level_support > 0.5] = 1 - tree_level_support[tree_level_support > 0.5]
site_level_support = OmegaCor[[2]]$support
site_level_support[site_level_support > 0.5] = 1 - site_level_support[site_level_support > 0.5]

corrplot(OmegaCor[[1]]$mean, method = "color", order = "hclust", type = "upper", diag = F,
col=colorRampPalette(c("blue","white","red"))(500),
#p.mat = tree_level_support, sig.level = 0.25,
#insig = "blank", #OR
#insig = "pch" , pch = ".", pch.cex = 0.05, pch.col = "white",
tl.pos = "n",
bg = "grey",
#title=paste("random effect level:", m.spatial.full_covars$rLNames[1]),
mar=c(2,0,1,0))

corrplot(OmegaCor[[2]]$mean, method = "color", order = "hclust", type = "lower", diag = F, add = T,
col=colorRampPalette(c("blue","white","red"))(500),
#p.mat = site_level_support, sig.level = 0.25,
#insig = "blank", #OR
#insig = "pch" , pch = ".", pch.cex = 0.05, pch.col = "white",
tl.pos = "n",
bg = "grey",
#title=paste("random effect level:", m.spatial.full_covars$rLNames[2]),
mar=c(2,0,1,0))


pdf("HMSC_figs/full_comm_ge1K_seqs_per_sample_gt_33_freq.bin.spatial.residual_corr.site_lower_tree_upper.pdf", width = 6, height = 6)

corrplot(OmegaCor[[1]]$mean, method = "color", order = "hclust", type = "upper", diag = T,
col=colorRampPalette(c("blue","white","red"))(500),
p.mat = tree_level_support, sig.level = 0.25,
#insig = "blank", #OR
insig = "pch" , pch = "+", pch.cex = 0.5, pch.col = "white",
tl.pos = "n",
addrect = 2, rect.col = "black", rect.lwd = 5,
#bg = "grey",
#title=paste("random effect level:", m.spatial.full_covars$rLNames[1]),
mar=c(2,0,1,0))

corrplot(OmegaCor[[2]]$mean, method = "color", order = "hclust", type = "lower", diag = F, add = T,
col=colorRampPalette(c("blue","white","red"))(500),
p.mat = site_level_support, sig.level = 0.25,
#insig = "blank", #OR
insig = "pch" , pch = "+", pch.cex = 0.0005, pch.col = "white",
tl.pos = "n",
#bg = "grey",
#title=paste("random effect level:", m.spatial.full_covars$rLNames[2]),
mar=c(2.5,0,1,0))

dev.off()


pdf("HMSC_figs/full_comm_ge1K_seqs_per_sample_gt_33_freq.bin.spatial.residual_corr.site_lower_tree_upper.no_sig.pdf", width = 6, height = 6)

corrplot(OmegaCor[[1]]$mean, method = "color", order = "hclust", type = "upper", diag = T,
col=colorRampPalette(c("blue","white","red"))(500),
#p.mat = tree_level_support, sig.level = 0.25,
#insig = "blank", #OR
#insig = "pch" , pch = ".", pch.cex = 0.0005, pch.col = "white",
tl.pos = "n",
addrect = 2, rect.col = "black", rect.lwd = 5,
#bg = "grey",
#title=paste("random effect level:", m.spatial.full_covars$rLNames[1]),
mar=c(2,0,1,0))

corrplot(OmegaCor[[2]]$mean, method = "color", order = "hclust", type = "lower", diag = F, add = T,
col=colorRampPalette(c("blue","white","red"))(500),
#p.mat = site_level_support, sig.level = 0.25,
#insig = "blank", #OR
#insig = "pch" , pch = ".", pch.cex = 0.0005, pch.col = "white",
tl.pos = "n",
#bg = "grey",
#title=paste("random effect level:", m.spatial.full_covars$rLNames[2]),
mar=c(2.5,0,1,0))

dev.off()



require(scales)

pdf("HMSC_figs/full_comm_ge1K_seqs_per_sample_gt_33_freq.bin.spatial.residual_corr.site_lower_tree_upper.no_sig.pdf", width = 6, height = 6)



corrplot(rescale(OmegaCor[[1]]$support, from = c(-1,1)), method = "color", order = "hclust", type = "upper", diag = T,
col=colorRampPalette(c("blue","white","red"))(500),
#p.mat = tree_level_support, sig.level = 0.25,
#insig = "blank", #OR
#insig = "pch" , pch = ".", pch.cex = 0.0005, pch.col = "white",
tl.pos = "n",
#addgrid.col = "black",
addrect = 2, rect.col = "black", rect.lwd = 5,
#bg = "grey",
#title=paste("random effect level:", m.spatial.full_covars$rLNames[1]),
mar=c(2,0,1,0))

tree_corrs_sig = OmegaCor[[1]]$mean
tree_corrs_sig[tree_level_support > 0.05] = NA

corrplot(OmegaCor[[1]]$mean, method = "color", order = "hclust", type = "upper", diag = T, add = T,
col=colorRampPalette(c("blue","white","red"))(500),
p.mat = tree_level_support, sig.level = 0.25,
insig = "blank", #OR
#insig = "pch" , pch = ".", pch.cex = 0.0005, pch.col = "white",
tl.pos = "n",
#addrect = 2, rect.col = "black", rect.lwd = 5,
#bg = "grey",
#title=paste("random effect level:", m.spatial.full_covars$rLNames[1]),
mar=c(2,0,1,0))


corrplot(OmegaCor[[2]]$mean, method = "color", order = "hclust", type = "lower", diag = F, add = T,
col=colorRampPalette(c("blue","white","red"))(500),
#p.mat = site_level_support, sig.level = 0.25,
#insig = "blank", #OR
#insig = "pch" , pch = ".", pch.cex = 0.0005, pch.col = "white",
tl.pos = "n",
#bg = "grey",
#title=paste("random effect level:", m.spatial.full_covars$rLNames[2]),
mar=c(2.5,0,1,0))

dev.off()


######################################################
#GGPLOT and TIDYVERSE methods for corr plotting below
######################################################

supportLevel = 0.75

require(purrr))
#Tree level corr table manipulations

#Functions for lower triangle and hclust
get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}
reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
}

reorder_cormat(OmegaCor[[1]]$mean) %>% corrplot

################
#This is nice because it uses map to apply across the list, but might be easier to order the matrix if performing functions one at a time
tree_ASV_pairwise_resid = OmegaCor[[1]] %>%
    map(~get_lower_tri(.x)) %>%
    map(~rownames_to_column(data.frame(.x), var="ASVs1")) %>%
    map(~pivot_longer(.x, -ASVs1, "ASVs2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    mutate(
        sig_support =
            ifelse(support > .95 | support < 0.05, "0.95<P<0.05",
                ifelse(support >= 0.05 & support < 0.25, "0.75<P<0.25", "n.s.")
            )
    )
tree_ASV_pairwise_resid = tree_ASV_pairwise_resid %>% mutate(mean_if_sig = ifelse(sig_support == "n.s.", 0, mean))
###################

tree_ASV_pairwise_resid.ordered.lower_tri = OmegaCor[[1]]$mean %>%
    reorder_cormat %>%
    get_lower_tri %>%
    data.frame %>%
    rownames_to_column(var = "ASVs1") %>%
    pivot_longer(., -ASVs1, "ASVs2") %>%
    filter(!is.na(value))

tree_ASV_pairwise_support.lower_tri = OmegaCor[[1]]$support %>%
get_lower_tri %>%
data.frame %>%
rownames_to_column(var = "ASVs1") %>%
pivot_longer(., -ASVs1, "ASVs2") %>%
filter(!is.na(value))

tree_ASV_pairwise_resid.ordered.lower_tri.support =
    left_join(
        tree_ASV_pairwise_resid.ordered.lower_tri,
        tree_ASV_pairwise_support.lower_tri,
        by = c("ASVs1", "ASVs2")
    ) %>%
    mutate(
       sig_support =
            ifelse(value.y > .95 | value.y < 0.05, "0.95<P<0.05",
                ifelse(value.y >= 0.05 & value.y < 0.25, "0.75<P<0.25", "n.s.")
        )
    ) %>%
    filter(!is.na(sig_support)) %>%
    mutate(mean_if_sig = ifelse(sig_support == "n.s.", 0, value.x))


tree_ASV_pairwise_resid.ordered.lower_tri.support$ASVs1 = factor(tree_ASV_pairwise_resid.ordered.lower_tri.support$ASVs1, levels = unique(tree_ASV_pairwise_resid.ordered.lower_tri.support$ASVs1))

tree_ASV_pairwise_resid.ordered.lower_tri.support$ASVs2 = factor(tree_ASV_pairwise_resid.ordered.lower_tri.support$ASVs2, levels = unique(tree_ASV_pairwise_resid.ordered.lower_tri.support$ASVs2))

ggplot(tree_ASV_pairwise_resid.ordered.lower_tri.support, aes(ASVs1, ASVs2, fill = mean_if_sig, color = sig_support)) +
geom_tile() +
scale_fill_gradient2(low = "#2c7bb6", high = "#d7191c", mid = "white", midpoint = 0) +
scale_color_manual(values = c("0.95<P<0.05" = "black", "0.75<P<0.25" = "white", "n.s." = "white"), guide = F) +
my_gg_theme +
theme(
axis.text = element_blank()
)



