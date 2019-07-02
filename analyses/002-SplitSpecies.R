################################################################################
###                                                                          ###
###                             SPLIT SPECIES LIST                           ###
###                                                                          ###
###--------------------------------------------------------------------------###
###                                                                          ###
### AUTHORS : Nicolas Casajus                                                ###
### DATE    : 2019/07/02                                                     ###
###                                                                          ###
###--------------------------------------------------------------------------###
###                                                                          ###
### - This R script splits species list in 3 users (old-school               ###
###   parallelization).                                                      ###
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


### Import Species Synthesis      -------------------

splist <- read.delim(
  paste(
    root_path,
    "data",
    "raw",
    "species",
    "species_list_all.txt",
    sep = .Platform$file.sep
  )
)


### Split Species List            -------------------

n_of_clusters <- 3
cluster_names <- c("NC", "MG", "NL")

species_list <- list()
for (i in 1:n_of_clusters) {
  species_list[[i]] <- vector()
  names(species_list)[i] <- cluster_names[i]
}

occs <- list()
occs[[1]] <- c(10, 100)
occs[[2]] <- c(101, 1000)
occs[[3]] <- c(1001, 10000)
occs[[4]] <- c(10001, 50000)

for (i in 1:length(occs)) {

  species <- as.character(
    splist[which(
      splist[ , "occurrence"] >= occs[[i]][1] &
      splist[ , "occurrence"] <= occs[[i]][2]
    ),
    "species"]
  )

  n_sp <- floor(length(species) / n_of_clusters)

  species <- sample(x = species, size = length(species), replace = FALSE)

  for (j in 1:length(cluster_names)) {

    if (j < n_of_clusters) {

      species_list[[j]] <- c(species_list[[j]], species[1:n_sp])
      species <- species[-c(1:n_sp)]

    } else {

      species_list[[j]] <- c(species_list[[j]], species[1:length(species)])
    }
  }
}

for (i in 1:n_of_clusters) {


### Export Species List           -------------------

  species <- splist[splist[ , "species"] %in% species_list[[i]], ]

  save(
    species,
    file = paste(
      root_path,
      "data",
      "ready",
      "species",
       paste0("species_list_", names(species_list)[i]),
       sep = .Platform$file.sep
    )
  )
}
