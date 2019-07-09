################################################################################
###                                                                          ###
###                                 RUN BIOMOD                               ###
###                                                                          ###
###--------------------------------------------------------------------------###
###                                                                          ###
### AUTHORS : Nicolas Casajus                                                ###
### DATE    : 2019/07/04                                                     ###
###                                                                          ###
###--------------------------------------------------------------------------###
###                                                                          ###
### - This R script runs Biomod and extract TSS and SRC.                     ###
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
### [1] parallel stats graphics grDevices utils datasets methods base        ###
###                                                                          ###
### other attached packages:                                                 ###
### [1]  doParallel_1.0.14 iterators_1.0.10 foreach_1.4.4 biomod2_3.3-7.1    ###
### [5]  ggplot2_3.1.0 reshape_0.8.8 raster_2.9-5 rgeos_0.4-3 rgdal_1.4-4    ###
### [10] sp_1.3-1                                                            ###
###                                                                          ###
################################################################################


### Load Addings                  -------------------

library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(biomod2)
library(foreach)
library(doParallel)

source(paste0(root_path, "/R/create_buffer.R"))


### Save Current path             -------------------

saved_path <- getwd()


### Create Folders                -------------------

dir.create(
  paste(
    root_path,
    "data",
    "ready",
    "species",
    "_occ",
    sep = .Platform$file.sep
  ),
  showWarnings  = FALSE,
  recursive     = TRUE
)

dir.create(
  paste(
    root_path,
    "data",
    "ready",
    "species",
    "_mask",
    sep = .Platform$file.sep
  ),
  showWarnings  = FALSE,
  recursive     = TRUE
)

dir.create(
  paste(
    root_path,
    "outputs",
    "projs",
    sep = .Platform$file.sep
  ),
  showWarnings  = FALSE,
  recursive     = TRUE
)

dir.create(
  paste(
    root_path,
    "outputs",
    "metrics",
    sep = .Platform$file.sep
  ),
  showWarnings  = FALSE,
  recursive     = TRUE
)

dir.create(
  paste(
    root_path,
    "outputs",
    "archives",
    sep = .Platform$file.sep
  ),
  showWarnings  = FALSE,
  recursive     = TRUE
)


### Import Study Area Grid        -------------------

study <- raster(
  paste(
    root_path,
    "data",
    "grid",
    "grid_area_with_cell_ids.tif",
    sep = .Platform$file.sep
  )
)


### Import Species List           -------------------

splist <- get(
  load(
    paste(
      root_path,
      "data",
      "ready",
      "species",
      paste0("species_list_", user),
      sep = .Platform$file.sep
    )
  )
)

splist  <- splist[order(splist[ , "occurrence"]), ]
spnames <- as.character(splist[ , "species"])


### Import Current Climate        -------------------

vars <- get(
  load(
    paste(
      root_path,
      "data",
      "ready",
      "climate",
      "1979-2013",
      sep = .Platform$file.sep
    )
  )
)


### Loop on Species               -------------------

registerDoParallel(cores = n_cores)
foreach(spname = spnames) %dopar% {

  tars <- list.files(
    path = paste(
      root_path,
      "outputs",
      "archives",
      sep = .Platform$file.sep
    ),
    full.names = FALSE
  )
  tars <- gsub("\\.tar\\.gz2", "", tars)

  pos <- which(tars == spname)

  if (length(pos) == 0) {

    taxa <- as.character(splist[splist[ , "species"] == spname, "group"])


  ### Import Species Occurrence     -------------------

    spocc <- get(
      load(
        paste(
          root_path,
          "data",
          "raw",
          "species",
          paste0("occurrences_", taxa),
          sep = .Platform$file.sep
        )
      )
    )

    spocc <- spocc[[spname]]
    cells <- data.frame(
      cell_id     = spocc,
      occurrence  = 1
    )


  ### Convert Sp Occurrence to SPDF -------------------

    pos <- which(!is.na(study[]))
    cell_id <- study[][pos]

    xy <- xyFromCell(study, pos)

    grille <- data.frame(
      cell_id = cell_id,
      x       = xy[ , 1],
      y       = xy[ , 2]
    )

    pts <- merge(cells, grille, by = "cell_id", all = TRUE)

    pts$occurrence <- ifelse(is.na(pts$occurrence), 0, 1)

    pts <- SpatialPointsDataFrame(
      coords      = pts[ , c("x", "y")],
      data        = data.frame(occurrence = pts[ , "occurrence"]),
      proj4string = CRS(projection(study))
    )


  ### Select Buffer Absences        -------------------

    ptsBuf <- create_buffer(sp = pts, distMin = 0, distMax = get(paste0("dist_", taxa)))

    pts0 <- as.data.frame(pts@coords[ptsBuf, c("x", "y")])
    pts1 <- as.data.frame(pts@coords[pts@data == 1, c("x", "y")])

    pts0$occurrence <- rep(0, nrow(pts0))
    pts1$occurrence <- rep(1, nrow(pts1))

    spocc <- rbind(pts0, pts1)

    save(
      spocc,
      file = paste(
        root_path,
        "data",
        "ready",
        "species",
        "_occ",
         paste0("occ_", spname),
         sep = .Platform$file.sep
      )
    )


  ### Create Selection Mask         -------------------

    xyz <- data.frame(pts@coords, z = 0)
    xyz[c(ptsBuf, which(pts@data == 1)), "z"] <- 1
    spmask <- rasterFromXYZ(xyz, res = 50000, crs = CRS(projection(study)))

    save(
      spmask,
      file = paste(
        root_path,
        "data",
        "ready",
        "species",
        "_mask",
         paste0("mask_", spname),
         sep = .Platform$file.sep
      )
    )


  ### Sample Species 0 AND 1        -------------------

    n1  <- sum(spocc[ , "occurrence"])

    sp1 <- spocc[spocc[ , "occurrence"] == 1, ]
    sp0 <- spocc[spocc[ , "occurrence"] == 0, ]

    # Define rules

    if (n1 < 101) {

      prevalence  <- 0.2
      nb_presence <- 1.0

    } else {

      if (n1 < 1001) {

        prevalence  <- 0.3
        nb_presence <- 0.9

      } else {

        if (n1 < 10001) {

          prevalence  <- 0.4
          nb_presence <- 0.8

        } else {

          prevalence  <- 0.5
          nb_presence <- 0.7

        }
      }
    }

    # Sample Presences

    pos1 <- sample(
      x        = 1:n1,
      size     = n1 * nb_presence,
      replace  = FALSE
    )

    sp1 <- sp1[pos1, ]

    n1  <- sum(sp1[ , "occurrence"])


    # Sample Absences

    n0  <- (n1 / prevalence) - n1

    pos0 <- sample(
      x        = 1:nrow(sp0),
      size     = n0,
      replace  = FALSE
    )

    sp0 <- sp0[pos0, ]


  ### Create BIOMOD Dataset         -------------------

    spocc <- rbind(sp1, sp0)

    sp.occ <- spocc[ , "occurrence"]
    sp.xy  <- spocc[ , c("x", "y")]


  ### Format Data for BIOMOD        -------------------

    bm.form <- BIOMOD_FormatingData(
      resp.var   = sp.occ,
      expl.var   = vars,
      resp.xy    = sp.xy,
      resp.name  = spname
    )


  ### Change Working Directory      -------------------

    setwd(
      paste(
        root_path,
        "outputs",
        sep = .Platform$file.sep
      )
    )


  ### Build BIOMOD Single Models    -------------------

    bm.mod <- BIOMOD_Modeling(
      data              = bm.form,
      models            = mod.models,
      models.options    = bm.opt,
      NbRunEval         = mod.n.rep,
      DataSplit         = mod.data.split,
      Prevalence        = prevalence,
      VarImport         = mod.var.import,
      models.eval.meth  = mod.models.eval.meth,
      do.full.models    = FALSE,
      modeling.id       = "biomod"
    )


  ### Build BIOMOD Ensemble Models  -------------------

    bm.em.all <- BIOMOD_EnsembleModeling(
      modeling.output                = bm.mod,
      chosen.models                  = 'all',
      em.by                          = 'all',
      eval.metric                    = ens.eval.metric,
      eval.metric.quality.threshold  = ens.eval.metric.quality.threshold,
      models.eval.meth               = ens.models.eval.meth,
      prob.mean                      = FALSE,
      prob.mean.weight               = ens.prob.mean.weight,
      prob.mean.weight.decay         = ens.prob.mean.weight.decay,
      committee.averaging            = ens.committee.averaging
    )


  ### List Scenarios (Cur + Fut)    -------------------

    scenarios <- list.files(
      paste(
        root_path,
        "data",
        "ready",
        "climate",
        sep = .Platform$file.sep
      )
    )


    for (scenario in scenarios) {


  ### Open Future Scenario          -------------------

      preds <- get(
        load(
          paste(
            root_path,
            "data",
            "ready",
            "climate",
            scenario,
            sep = .Platform$file.sep
          )
        )
      )


  ### Ensemble Future Projection    -------------------

      bm.ef.all <- BIOMOD_EnsembleForecasting(
        EM.output       = bm.em.all,
        new.env         = preds,
        output.format   = ".grd",
        proj.name       = paste0(scenario, "_ENSEMBLE_all"),
        selected.models = "all",
        binary.meth     = ens.models.eval.meth
      )
    }


  ### Get Binary Cutoff             -------------------

    cutoff <- get(
      load(
        paste(
          root_path,
          "outputs",
          spname,
          "models",
          "biomod",
          paste0(
            spname,
            "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData"
          ),
          sep = .Platform$file.sep
        )
      )
    )

    cutoff <- cutoff@model_evaluation["TSS", "Cutoff"]

    save(
      cutoff,
      file = paste(
        root_path,
        "outputs",
        "metrics",
         paste("cutoff", spname, sep = "_"),
         sep = .Platform$file.sep
      )
    )


  ### Store Evaluation metric       -------------------

    tss <- bm.mod@models.evaluation@val[ , "Testing.data", , , ]

    save(
      tss,
      file = paste(
        root_path,
        "outputs",
        "metrics",
         paste("tss", spname, sep = "_"),
         sep = .Platform$file.sep
      )
    )


  ### Mega Ensemble Projection      -------------------

    for (horizon in c(baseline, horizons)) {

      projs <- list.files(
        path        = paste(
          root_path,
          "outputs",
          spname,
          sep = .Platform$file.sep
        ),
        full.names  = TRUE,
        recursive   = TRUE,
        pattern     = paste0("^proj_", horizon, ".+_ensemble.grd$")
      )

      projs  <- stack(projs)
      projs  <- subset(projs, grep("wmean", names(projs)))
      projs  <- mean(projs)
      bins   <- projs

      bins[] <- ifelse(projs[] < cutoff, 0, 1)

      names(projs) <- paste(spname, horizon, sep = "_")
      names(bins)  <- paste(spname, horizon, sep = "_")

      save(
        bins,
        file = paste(
          root_path,
          "outputs",
          "projs",
           paste("proj", spname, horizon, "bins", sep = "_"),
           sep = .Platform$file.sep
        )
      )

      save(
        projs,
        file = paste(
          root_path,
          "outputs",
          "projs",
           paste("proj", spname, horizon, "probs", sep = "_"),
           sep = .Platform$file.sep
        )
      )
    }


  ### Apply Mask on Current Proj    -------------------

  proj <- get(
    load(
      paste(
        root_path,
        "outputs",
        "projs",
         paste("proj", spname, baseline, "bins", sep = "_"),
         sep = .Platform$file.sep
      )
    )
  )

  mask <- get(
    load(
      paste(
        root_path,
        "data",
        "ready",
        "species",
        "_mask",
         paste("mask", spname, sep = "_"),
         sep = .Platform$file.sep
      )
    )
  )

  proj <- proj * mask

  save(
    proj,
    file = paste(
      root_path,
      "outputs",
      "projs",
       paste("proj", spname, baseline, "bins", sep = "_"),
       sep = .Platform$file.sep
    )
  )


  ### Define Mask for Future Proj   -------------------

  pos <- which(!is.na(proj[]))
  pts <- SpatialPointsDataFrame(
    coords  = xyFromCell(proj, pos),
    data    = data.frame(proj[pos])
  )
  proj4string(pts) <- projection(proj)

  ptsBuf <- create_buffer(sp = pts, distMin = 0, distMax = get(paste0("dist_", taxa)))

  mask   <- proj
  mask[] <- NA

  mask[cellFromXY(mask, pts@coords[ , 1:2])] <- 0
  mask[cellFromXY(mask, pts@coords[pts@data == 1, 1:2])] <- 1
  mask[cellFromXY(mask, pts@coords[ptsBuf, 1:2])] <- 1


  ### Constrain Future Projections  -------------------

  for (horizon in horizons) {

    proj <- get(
      load(
        paste(
          root_path,
          "outputs",
          "projs",
           paste("proj", spname, horizon, "bins", sep = "_"),
           sep = .Platform$file.sep
        )
      )
    )

    proj <- proj * mask

    save(
      proj,
      file = paste(
        root_path,
        "outputs",
        "projs",
         paste("proj", spname, horizon, "bins", sep = "_"),
         sep = .Platform$file.sep
      )
    )
  }


  ### Compute Species Range Change  -------------------

  cur <- get(
    load(
      paste(
        root_path,
        "outputs",
        "projs",
         paste("proj", spname, baseline, "bins", sep = "_"),
         sep = .Platform$file.sep
      )
    )
  )

  for (horizon in horizons) {

    fut <- get(
      load(
        paste(
          root_path,
          "outputs",
          "projs",
           paste("proj", spname, horizon, "bins", sep = "_"),
           sep = .Platform$file.sep
        )
      )
    )

    src <- BIOMOD_RangeSize(cur, fut)

    save(
      src,
      file = paste(
        root_path,
        "outputs",
        "metrics",
         paste("src", spname, horizon, "bins", sep = "_"),
         sep = .Platform$file.sep
      )
    )
  }


  ### Compress Biomod Results       -------------------

    system(
      paste0(
        "tar jcf ",
        "archives/",
        spname,
        ".tar.gz2 ",
        spname, "/"
      )
    )

    system(
      paste0(
        "rm -r ",
        spname, "/"
      )
    )
  }
}


### Reset Working Directory       -------------------

setwd(saved_path)
