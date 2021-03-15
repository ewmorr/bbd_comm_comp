## Repo for analysis of fungal community composition and diversity in bark endophyte communities of American beech (12/18/2020)

- Bioinformatic analyses and ASV calling were performed with code from the repo neonectria_barcoding_012220
- This repo picks up from after ASV calling, LULU post-processing
- The majority of analyses will use data after combining samples at the tree level (a random subset of trees had multiple plugs collected; the effects of combining multiple plugs ), but some analyses (e.g., analysis of community composition changes at different scales) will use data at the plug level
- data tables are read from the local directory `GARNAS_neonectria_barcoding_files_cat_03242020` and results written to the folder `GARNAS_bbd_comm_comp_12182020`
- Many initial analyses on community composition have already been performed and reside in `GARNAS_neonectria_barcoding_files_cat_03242020`
- The scripts for those results have been copied to this folder in most cases, but retain the original file path for writing results

### Q1: What are the primary drivers of community composition––are tree level or site level covariates stronger predictors.
- Q1a: How much variance is explained by meaured covariates versus random effects (site level and tree level)
- Q1b: are traits (such as from Polme et al. 2021) associated with variance
- See Ovaskainen et al. 2017 example 1 in the paper for HMSC approach to addressing these
-
    - Tried HMSC with all climated vars (growing and nongrowing season etc) and tree level vars but excluding Neonectria fruiting. Site level vars appear to explain more variation on average, and there is much more residual variation at the tree level. Try adding Neonectria fruting. This var explains a lot of variation in Nf, but may be important for tree level "biotic state"

### Q2 What are the primary dirvers of diversity––are tree evel or site level covariates stronger predictos.
- Can use a similar approach as Q1 but using univariate linear mixed effects models



