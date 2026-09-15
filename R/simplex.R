.simplex_sf_method = \(data, column, target, E = 2:10, k = E+2, tau = 1, style = 1, stack = FALSE, lib = NULL, pred = NULL,
                       dist.metric = "L2", dist.average = TRUE, threads = detectThreads(), detrend = TRUE, nb = NULL){
  vx = .uni_lattice(data,column,detrend)
  vy = .uni_lattice(data,target,detrend)
  if (is.null(lib)) lib = .internal_library(cbind(vx,vy))
  if (is.null(pred)) pred = lib
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  res = RcppSimplex4Lattice(vx,vy,nb,lib,pred,E,k,tau,style,stack,
                            .check_distmetric(dist.metric),dist.average,threads)
  return(.bind_xmapself(res,target,"simplex"))
}

.simplex_spatraster_method = \(data, column, target, E = 2:10, k = E+2, tau = 1, style = 1, stack = FALSE, lib = NULL, pred = NULL, dist.metric = "L2",
                               dist.average = TRUE, threads = detectThreads(), detrend = TRUE, grid.coord = TRUE, embed.direction = 0){
  mx = .uni_grid(data,column,detrend,grid.coord)
  my = .uni_grid(data,target,detrend,grid.coord)
  if (is.null(lib)) lib = which(!(is.na(mx) | is.na(my)), arr.ind = TRUE)
  if (is.null(pred)) pred = lib
  res = RcppSimplex4Grid(mx, my, lib, pred, E, k, tau, style, stack,
                         .check_distmetric(dist.metric),dist.average,threads,embed.direction)
  return(.bind_xmapself(res,target,"simplex"))
}

#' optimal parameter search for simplex projection
#'
#' @inheritParams embedded
#' @param column name of library variable.
#' @param k (optional) number of nearest neighbors used.
#' @param lib (optional) libraries indices (input needed: `vector` - spatial vector, `matrix` - spatial raster).
#' @param pred (optional) predictions indices (input requirement same as `lib`).
#' @param dist.metric (optional) distance metric (`L1`: Manhattan, `L2`: Euclidean).
#' @param dist.average (optional) whether to average distance.
#' @param threads (optional) number of threads to use.
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{forecast performance}
#' \item{\code{varname}}{name of target variable}
#' \item{\code{method}}{method of cross mapping}
#' }
#' @export
#' @name simplex
#' @aliases simplex,sf-method
#' @references
#' Sugihara G. and May R. 1990. Nonlinear forecasting as a way of distinguishing chaos from measurement error in time series. Nature, 344:734-741.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg",package="spEDM"))
#' \donttest{
#' spEDM::simplex(columbus,"inc","crime")
#' }
methods::setMethod("simplex", "sf", .simplex_sf_method)

#' @rdname simplex
methods::setMethod("simplex", "SpatRaster", .simplex_spatraster_method)
