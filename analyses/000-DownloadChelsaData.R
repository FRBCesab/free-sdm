################################################################################
###                                                                          ###
###                      DOWNLOAD CHELSA CLIMATE LAYERS                      ###
###                                                                          ###
###--------------------------------------------------------------------------###
###                                                                          ###
### AUTHORS : Nicolas Casajus                                                ###
### DATE    : 2019/06/21                                                     ###
###                                                                          ###
###--------------------------------------------------------------------------###
###                                                                          ###
### - This R script downloads climate normals from the Chelsa database for   ###
###   19 bioclimatic variables, several horizons, GCMs and RCPs.             ###
###                                                                          ###
###--------------------------------------------------------------------------###
###                                                                          ###
### > sessionInfo()                                                          ###
###                                                                          ###
### R version 3.5.3 (2019-03-11)                                             ###
### Platform: x86_64-apple-darwin15.6.0 (64-bit)                             ###
### Running under: macOS Mojave 10.14.5                                      ###
###                                                                          ###
### locale:                                                                  ###
### [1] fr_FR.UTF-8/fr_FR.UTF-8/fr_FR.UTF-8/C/fr_FR.UTF-8/fr_FR.UTF-8        ###
###                                                                          ###
### attached base packages:                                                  ###
### [1] stats     graphics  grDevices utils     datasets  methods   base     ###
###                                                                          ###
################################################################################


### Load Addings                  -------------------

source(paste0(root_path, "R/list_chelsa_cmip5.R"))
source(paste0(root_path, "R/wget_chelsa.R"))


### Variables Definition          -------------------

out_path <- paste0(root_path, "data/raw/climate")
biovars  <- c(1, 7, 12, 15)
horizons <- c("1979-2013", "2041-2060", "2061-2080")
gcm_ids  <- c("CESM1-BGC", "MIROC5", "CMCC-CMS", "IPSL-CM5A-LR", "MPI-ESM-MR")
rcp_ids  <- "rcp85"


### Download Climate Layers       -------------------

wget_chelsa(out_path, biovars, horizons, gcm_ids, rcp_ids)
