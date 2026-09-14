.pc_sf_method = \(data, column, target, E = 2:10, k = E+2, tau = 1, style = 1, lib = NULL, pred = NULL, dist.metric = "L2", zero.tolerance = max(k),
                  relative = TRUE, weighted = TRUE, maximize = "dark", threads = detectThreads(), detrend = FALSE, nb = NULL){
  vx = .uni_lattice(data,column,detrend)
  vy = .uni_lattice(data,target,detrend)
  if (is.null(lib)) lib = .internal_library(cbind(vx,vy))
  if (is.null(pred)) pred = lib
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  res = RcppPC4Lattice(vx, vy, nb, lib, pred, E, k, tau, style, zero.tolerance,
                       .check_distmetric(dist.metric),relative,weighted,threads,1)
  return(.bind_xmapself(res,column,"pc",maximize = maximize))
}

.pc_spatraster_method = \(data, column, target, E = 2:10, k = E+2, tau = 1, style = 1, lib = NULL, pred = NULL, dist.metric = "L2", zero.tolerance = max(k),
                          relative = TRUE, weighted = TRUE, maximize = "dark", threads = detectThreads(), detrend = FALSE, grid.coord = TRUE, embed.direction = 0){
  mx = .uni_grid(data,column,detrend,grid.coord)
  my = .uni_grid(data,target,detrend,grid.coord)
  if (is.null(lib)) lib = which(!(is.na(mx) | is.na(my)), arr.ind = TRUE)
  if (is.null(pred)) pred = lib
  res = RcppPC4Grid(mx, my, lib, pred, E, k, tau, style, zero.tolerance,
                    .check_distmetric(dist.metric),relative,weighted,threads,1,embed.direction)
  return(.bind_xmapself(res,column,"pc",maximize = maximize))
}

#' optimal parameter search for pattern causality
#'
#' @inheritParams simplex
#' @param zero.tolerance (optional) maximum number of zeros tolerated in signature space.
#' @param relative (optional) whether to calculate relative changes in embeddings.
#' @param weighted (optional) whether to weight causal strength.
#' @param maximize (optional) causality metric to maximize: one of "positive", "negative", or "dark".
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{cross mapping performance}
#' \item{\code{varname}}{name of target variable}
#' \item{\code{method}}{method of cross mapping}
#' \item{\code{maximize}}{maximized causality metric}
#' }
#' @export
#' @name pc
#' @aliases pc,sf-method
#' @references
#' Stavroglou, S.K., Pantelous, A.A., Stanley, H.E., Zuev, K.M., 2020. Unveiling causal interactions in complex systems. Proceedings of the National Academy of Sciences 117, 7599–7605.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg",package="spEDM"))
#' \donttest{
#' spEDM::pc(columbus,"crime","hoval",E = 5:10,maximize = "negative")
#' }
methods::setMethod("pc", "sf", .pc_sf_method)

#' @rdname pc
methods::setMethod("pc", "SpatRaster", .pc_spatraster_method)
