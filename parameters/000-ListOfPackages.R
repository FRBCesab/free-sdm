pkgs <- c("sp", "rgdal", "rgeos", "raster", "biomod2", "foreach", "doParallel")

for (pkg in pkgs) {

  install.packages(pkg)  
}
