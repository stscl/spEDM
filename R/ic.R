.ic_sf_method = \(data, column, target, E = 2:10, k = E+2, tau = 1, style = 1, lib = NULL, pred = NULL,
                  dist.metric = "L2", threads = detectThreads(), detrend = FALSE, nb = NULL){
  vx = .uni_lattice(data,column,detrend)
  vy = .uni_lattice(data,target,detrend)
  if (is.null(lib)) lib = .internal_library(cbind(vx,vy))
  if (is.null(pred)) pred = lib
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  res = RcppIC4Lattice(vx, vy, nb, lib, pred, E, k, tau, 0, style,
                       .check_distmetric(dist.metric),threads,0)
  return(.bind_xmapself(res,target,"ic"))
}

.ic_spatraster_method = \(data, column, target, E = 2:10, k = E+2, tau = 1, style = 1, lib = NULL, pred = NULL,
                          dist.metric = "L2", threads = detectThreads(), detrend = FALSE, grid.coord = TRUE, embed.direction = 0){
  mx = .uni_grid(data,column,detrend,grid.coord)
  my = .uni_grid(data,target,detrend,grid.coord)
  if (is.null(lib)) lib = which(!(is.na(mx) | is.na(my)), arr.ind = TRUE)
  if (is.null(pred)) pred = lib
  res = RcppIC4Grid(mx, my, lib, pred, E, k, tau, 0, style,
                    .check_distmetric(dist.metric),threads,0,embed.direction)
  return(.bind_xmapself(res,target,"ic"))
}

#' optimal parameter search for intersectional cardinality
#'
#' @inheritParams simplex
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{cross mapping performance}
#' \item{\code{varname}}{name of target variable}
#' \item{\code{method}}{method of cross mapping}
#' }
#' @export
#' @name ic
#' @aliases ic,sf-method
#' @references
#' Tao, P., Wang, Q., Shi, J., Hao, X., Liu, X., Min, B., Zhang, Y., Li, C., Cui, H., Chen, L., 2023. Detecting dynamical causality by intersection cardinal concavity. Fundamental Research.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg",package="spEDM"))
#' \donttest{
#' spEDM::ic(columbus,"hoval","crime",E = 7,k = 15:25)
#' }
methods::setMethod("ic", "sf", .ic_sf_method)

#' @rdname ic
methods::setMethod("ic", "SpatRaster", .ic_spatraster_method)
