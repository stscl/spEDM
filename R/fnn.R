.fnn_sf_method = \(data, target, E = 1:10, tau = 1, style = 1, stack = FALSE, lib = NULL, pred = NULL, 
                   dist.metric = "L1", rt = 10, eps = 2, threads = detectThreads(), detrend = TRUE, nb = NULL){
  vec = .uni_lattice(data,target,detrend)
  rt = .check_inputelementnum(rt,max(E))
  eps = .check_inputelementnum(eps,max(E))
  if (is.null(lib)) lib = which(!is.na(vec))
  if (is.null(pred)) pred = lib
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  return(RcppFNN4Lattice(vec, nb, rt, eps, lib, pred, E, tau, style,
                         stack, .check_distmetric(dist.metric), threads))
}

.fnn_spatraster_method = \(data, target, E = 1:10, tau = 1, style = 1, stack = FALSE, lib = NULL, pred = NULL, dist.metric = "L1", 
                           rt = 10, eps = 2, threads = detectThreads(), detrend = TRUE, grid.coord = TRUE, embed.direction = 0){
  mat = .uni_grid(data,target,detrend,grid.coord)
  rt = .check_inputelementnum(rt,max(E))
  eps = .check_inputelementnum(eps,max(E))
  if (is.null(lib)) lib = which(!is.na(mat), arr.ind = TRUE)
  if (is.null(pred)) pred = lib
  return(RcppFNN4Grid(mat, rt, eps, lib, pred, E, tau, style, stack,
                      .check_distmetric(dist.metric),threads,0,embed.direction))
}

#' false nearest neighbours
#'
#' @inheritParams embedded
#' @param lib (optional) libraries indices (input needed: `vector` - spatial vector, `matrix` - spatial raster).
#' @param pred (optional) predictions indices (input requirement same as `lib`).
#' @param dist.metric (optional) distance metric (`L1`: Manhattan, `L2`: Euclidean).
#' @param rt (optional) escape factor.
#' @param eps (optional) neighborhood diameter.
#' @param threads (optional) number of threads to use.
#'
#' @return A vector
#' @export
#' @name fnn
#' @aliases fnn,sf-method
#' @references
#' Kennel, M.B., Brown, R., Abarbanel, H.D.I., 1992. Determining embedding dimension for phase-space reconstruction using a geometrical construction. Physical Review A 45, 3403–3411.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg",package="spEDM"))
#' \donttest{
#' spEDM::fnn(columbus,"crime")
#' }
methods::setMethod("fnn", "sf", .fnn_sf_method)

#' @rdname fnn
methods::setMethod("fnn", "SpatRaster", .fnn_spatraster_method)
