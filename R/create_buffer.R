#' @title Create buffer around species distribution
#'
#' @description
#' This function selects absences inside a buffer near the species distribution.
#'
#' @param sp a SpatialPointsDataFrame.
#' @param distMin a numeric indicating the minimum distance of the buffer (in meters).
#' @param distMax a numeric indicating the maximum distance of the buffer (in meters).
#'
#' @author Maya Guégan, \email{maya.gueguen@@univ-grenoble-alpes.fr}
#'
#' @export
#'
#' @return A vector with the rows index of the selected absences.
#'
#' @examples
#'
#' create_buffer(
#'   sp      = spdf,
#'   distMin = 0,
#'   distMin = 3000000
#' )


create_buffer <- function(sp, distMin, distMax) {

  # Determining selectable area

  coor    <- coordinates(sp)

  pres    <- which(sp@data[ , 1] == 1)
  abs     <- which(sp@data[ , 1] == 0)

  outside <- rep(0, length(abs))
  inside  <- rep(0, length(abs))

  for (i in 1:length(pres)) {

    pt_dist <- sqrt((coor[abs, 1] - coor[pres[i], 1]) ^ 2 + (coor[abs, 2] - coor[pres[i], 2]) ^ 2)

    # removing points too close from presences

    inside <- inside + (pt_dist > distMin)

    # keeping points not to far from presences

    if (!is.null(distMax)) {

      outside <- outside + (pt_dist < distMax)
    }
  }

  # No cells are too far

  if (is.null(distMax)) {

    outside <- outside + 1
  }

  selected.abs <- abs[(inside == length(pres)) & (outside > 0)]

  return(selected.abs)
}
