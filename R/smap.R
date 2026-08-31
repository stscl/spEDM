.smap_sf_method = \(data, column, target, E = 3, k = E+2, tau = 1, style = 1, stack = FALSE, lib = NULL, pred = NULL, dist.metric = "L2", dist.average = TRUE,
                    theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8), threads = detectThreads(), detrend = TRUE, nb = NULL){
  vx = .uni_lattice(data,column,detrend)
  vy = .uni_lattice(data,target,detrend)
  if (is.null(lib)) lib = .internal_library(cbind(vx,vy))
  if (is.null(pred)) pred = lib
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  res = RcppSMap4Lattice(vx,vy,nb,lib,pred,theta,E,tau,k,style,stack,
                         .check_distmetric(dist.metric),dist.average, threads)
  return(.bind_xmapself(res,target,"smap"))
}

.smap_spatraster_method = \(data, column, target, E = 3, k = E+2, tau = 1, style = 1, stack = FALSE, lib = NULL, pred = NULL, dist.metric = "L2",
                            dist.average = TRUE, theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8),
                            threads = detectThreads(), detrend = TRUE, grid.coord = TRUE, embed.direction = 0){
  mx = .uni_grid(data,column,detrend,grid.coord)
  my = .uni_grid(data,target,detrend,grid.coord)
  if (is.null(lib)) lib = which(!(is.na(mx) | is.na(my)), arr.ind = TRUE)
  if (is.null(pred)) pred = lib
  res = RcppSMap4Grid(mx, my, lib, pred, theta, E, tau, k, style, stack,
                      .check_distmetric(dist.metric),dist.average,threads,embed.direction)
  return(.bind_xmapself(res,target,"smap"))
}

#' optimal parameter search for s-mapping
#'
#' @inheritParams simplex
#' @param theta (optional) weighting parameter for distances.
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{forecast performance}
#' \item{\code{varname}}{name of target variable}
#' \item{\code{method}}{method of cross mapping}
#' }
#' @export
#' @name smap
#' @aliases smap,sf-method
#' @references
#' Sugihara G. 1994. Nonlinear forecasting for the classification of natural time series. Philosophical Transactions: Physical Sciences and Engineering, 348 (1688):477-495.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg",package="spEDM"))
#' \donttest{
#' spEDM::smap(columbus,"inc","crime",E = 5,k = 6)
#' }
methods::setMethod("smap", "sf", .smap_sf_method)

#' @rdname smap
methods::setMethod("smap", "SpatRaster", .smap_spatraster_method)
