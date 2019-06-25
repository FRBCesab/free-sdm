################################################################################
###                                                                          ###
###                       UPSCALE CHELSA CLIMATE LAYERS                      ###
###                                                                          ###
###--------------------------------------------------------------------------###
###                                                                          ###
### AUTHORS : Nicolas Casajus                                                ###
### DATE    : 2019/06/25                                                     ###
###                                                                          ###
###--------------------------------------------------------------------------###
###                                                                          ###
### - This R script projects and upscales climate normals from the Chelsa    ###
###   database.                                                              ###
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
### other attached packages:                                                 ###
### [1] raster_2.8-19    sp_1.3-1                                            ###
###                                                                          ###
################################################################################


### Project Path                  -------------------

root_path <- "~/Desktop/free-sdm"


### Load Addings                  -------------------

library(sp)
library(raster)


### Import Grid Area              -------------------

grille <- raster(
  paste(
    root_path,
    "data",
    "grid",
    "grid_area_with_cell_ids.tif",
    sep = .Platform$file.sep
  )
)


### Import Clipper (NA for Bio15) -------------------

clipper <- readRDS(
  paste(
    root_path,
    "data",
    "grid",
    "clipper.rds",
    sep = .Platform$file.sep
  )
)


### List Climate Informations     -------------------

fls <- list.files(
  path        = paste0(root_path, "/data/raw/climate"),
  recursive   = TRUE,
  pattern     = "\\.tif$",
  full.names  = FALSE
)


### Get Variables Names           -------------------

layers <- gsub("\\/", "_", fls)
layers <- gsub("\\.tif$", "", layers)
layers <- gsub("CHELSA_bio10_", "Bio", layers)

layers <- unlist(
  lapply(
    strsplit(layers, "_"),
    function(x) paste0(c(x[length(x)], x[-length(x)]), collapse = "_")
  )
)


### Get Projections Names         -------------------

proj <- unique(substr(layers, 7, nchar(layers)))


### Core Program                  -------------------

for (i in 1:length(proj)) {

  pos <- grep(proj[i], layers)

  cat(paste0("### Upscaling layers for ", proj[i], " ...\n"))

  for (j in 1:length(pos)) {


### Import Biovar j for proj i    -------------------

    xxx <- raster(paste0(root_path, "/data/raw/climate/", fls[pos[j]]))
    names(xxx) <- substr(layers[pos[j]], 1, 5)


### Set Projection System         -------------------

    projection(xxx) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

    cat(paste0("    > ", names(xxx), "\n"))


### Correct NA Encoding           -------------------

    ocean <- which(xxx[] < -30000)
    xxx[][ocean] <- NA


### Project and Upscale Layer     -------------------

    xxx <- projectRaster(from = xxx, to = grille)


### Convert Units for Celsius     -------------------

    if (names(xxx) %in% c("Bio01", "Bio07")) {

      xxx[] <- xxx[] / 10
    }


### Convert NA to 0 (Bio15)       -------------------

    if (names(xxx) == "Bio15") {

      checks <- unlist(
        cellFromPolygon(
          object  = xxx,
          p       = clipper
        )
      )

      pos <- which(is.na(xxx[][checks]))

      if (length(pos) > 0) {

        xxx[][checks[pos]] <- 0
      }
    }


### Stack Layers by Projection    -------------------

    if (j == 1) {

      climate <- xxx

    } else {

      climate <- stack(climate, xxx)
    }
  }


### Export Stack                  -------------------

  save(
    climate,
    file = paste(
      root_path,
      "data",
      "ready",
      "climate",
       proj[i],
       sep = .Platform$file.sep
    )
  )

  cat("\n")
}
