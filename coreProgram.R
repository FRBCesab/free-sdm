
### Path to Git Project           -------------------
root_path <- "~/free-sdm"


### Username (Initials)           -------------------
### (Possible values: NC MG NL)
user      <- "NC"


### Number of Cores               -------------------
n_cores   <- 16

### Do Not Run                    -------------------
# source(paste0(root_path, "/analyses/000-DownloadChelsaData.R"))
# source(paste0(root_path, "/analyses/001-UpscaleChelsaData.R"))
# source(paste0(root_path, "/analyses/002-SplitSpecies.R"))


### Install Required Packages     -------------------
source(paste0(root_path, "/parameters/000-ListOfPackages.R"))


### Load Biomod Parameters        -------------------
source(paste0(root_path, "/parameters/001-BiomodParameters.R"))


### Run Core Program              -------------------
### (Model + Projections + SRC)
source(paste0(root_path, "/analyses/003-RunNicheModels.R"))


### That's it!                    -------------------
