### Horizons Names                -------------------

baseline <- "1979-2013"
horizons <- c("2041-2060", "2061-2080")


### Max. Dispersion (in meters)   -------------------

dist_mammals <- 3000000
dist_birds   <- 4000000


### Models Parameters             -------------------

mod.models           <- c("RF", "GLM", "GAM", "GBM")
mod.n.rep            <-  4
mod.data.split       <- 80
mod.var.import       <-  0
mod.models.eval.meth <- "TSS"


### Algo Options                  -------------------

bm.opt <- biomod2::BIOMOD_ModelingOptions(
  GLM = list(
    type              = "quadratic",
    interaction.level = 0,
    test              = "AIC"
  ),
  GBM = list(
    n.trees = 5000
  ),
  GAM = list(
    k = 3
  )
)


### Ensemble Parameters           -------------------

ens.eval.metric                   <- "TSS"
ens.eval.metric.quality.threshold <- 0.8
ens.models.eval.meth              <- "TSS"
ens.prob.mean.weight              <- TRUE
ens.prob.mean.weight.decay        <- "proportional"
ens.committee.averaging           <- TRUE
