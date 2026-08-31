.embedded_sf_method = \(data, target, E = 3, tau = 1, style = 1,
                        stack = FALSE,detrend = FALSE,nb = NULL){
  vec = .uni_lattice(data,target,detrend)
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  if (!stack) {
    res = RcppGenLatticeEmbeddings(vec,nb,E,tau,style)
  } else {
    res = RcppGenLatticeEmbeddingsCom(vec,nb,E,tau,style)
  }
  return(res)
}

.embedded_spatraster_method = \(data, target, E = 3, tau = 1, style = 1, stack = FALSE,
                                detrend = FALSE, grid.coord = TRUE, embed.direction = 0){
  mat = .uni_grid(data,target,detrend,grid.coord)
  if (!stack) {
    res = RcppGenGridEmbeddings(mat,E,tau,style,embed.direction)
  } else {
    res = RcppGenGridEmbeddingsCom(mat,E,tau,style,embed.direction)
  }
  return(res)
}

#' spatial-lag embedding
#'
#' @param data observation data.
#' @param target name of target variable.
#' @param E (optional) embedding dimensions.
#' @param tau (optional) step of spatial lags.
#' @param style (optional) embedding style (`0` includes current state, `1` excludes it).
#' @param stack (optional) whether to stack embeddings.
#' @param detrend (optional) whether to remove the linear trend.
#' @param nb (optional) neighbours list.
#'
#' @return A matrix (when `stack` is `FALSE`) or list
#' @export
#' @name embedded
#' @aliases embedded,sf-method
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg",package="spEDM"))
#' v = spEDM::embedded(columbus,"crime")
#' v[1:5,]
#'
#' npp = terra::rast(system.file("case/npp.tif",package="spEDM"))
#' r = spEDM::embedded(npp,"npp")
#' r[which(!is.na(r),arr.ind = TRUE)[1:5],]
#'
methods::setMethod("embedded", "sf", .embedded_sf_method)

#' @rdname embedded
#' @param grid.coord (optional) whether to detrend using cell center coordinates (`TRUE`) or row/column numbers (`FALSE`).
#' @param embed.direction (optional) direction selector for embeddings (`0` returns all directions, `1-8` correspond to NW, N, NE, W, E, SW, S, SE).
methods::setMethod("embedded", "SpatRaster", .embedded_spatraster_method)
