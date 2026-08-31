.gpc_sf_method = \(data, cause, effect, libsizes = NULL, E = 3, k = E+2, tau = 1, style = 1, lib = NULL, pred = NULL, boot = 99, replace = FALSE, seed = 42L, dist.metric = "L2", zero.tolerance = max(k),
                   relative = TRUE, weighted = TRUE, threads = detectThreads(), detrend = FALSE, parallel.level = "low", bidirectional = TRUE, progressbar = TRUE, nb = NULL){
  varname = .check_character(cause, effect)
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  cause = .uni_lattice(data,cause,detrend)
  effect = .uni_lattice(data,effect,detrend)
  if (is.null(lib)) lib = which(!(is.na(cause) | is.na(effect)))
  if (is.null(pred)) pred = lib
  return(.run_gpc(cause, effect, E, k, tau, style, lib, pred, .check_distmetric(dist.metric),
                  zero.tolerance, relative, weighted, threads, bidirectional, varname, nb,
                  libsizes, boot, replace, seed, parallel.level, progressbar))
}

.gpc_spatraster_method = \(data, cause, effect, libsizes = NULL, E = 3, k = E+2, tau = 1, style = 1, lib = NULL, pred = NULL, boot = 99, replace = FALSE, seed = 42L, dist.metric = "L2", zero.tolerance = max(k),
                           relative = TRUE, weighted = TRUE, threads = detectThreads(), detrend = FALSE, parallel.level = "low", bidirectional = TRUE, progressbar = TRUE, grid.coord = TRUE, embed.direction = 0){
  varname = .check_character(cause, effect)
  cause = .uni_grid(data,cause,detrend,grid.coord)
  effect = .uni_grid(data,effect,detrend,grid.coord)
  if (is.null(lib)) lib = which(!(is.na(cause) | is.na(effect)), arr.ind = TRUE)
  if (is.null(pred)) pred = lib
  return(.run_gpc(cause, effect, E, k, tau, style, lib, pred, .check_distmetric(dist.metric),
                  zero.tolerance, relative, weighted, threads, bidirectional, varname, NULL,
                  libsizes, boot, replace, seed, parallel.level, progressbar, embed.direction))
}

#' geographical pattern causality
#'
#' @inheritParams gcmc
#' @param boot (optional) number of bootstraps to perform.
#' @param seed (optional) random seed.
#' @param replace (optional) should sampling be with replacement?
#' @param zero.tolerance (optional) maximum number of zeros tolerated in signature space.
#' @param relative (optional) whether to calculate relative changes in embedding.
#' @param weighted (optional) whether to weight causal strength.
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{cross mapping results (only present if `libsizes` is not `NULL`)}
#' \item{\code{causality}}{per-sample causality statistics (present if `libsizes` is `NULL`)}
#' \item{\code{summary}}{overall causal strength (present if `libsizes` is `NULL`)}
#' \item{\code{pattern}}{pairwise pattern relationships (present if `libsizes` is `NULL`)}
#' \item{\code{varname}}{names of causal and effect variables}
#' \item{\code{bidirectional}}{whether to examine bidirectional causality}
#' }
#' @export
#' @name gpc
#' @aliases gpc,sf-method
#' @references
#' Zhang, Z., Wang, J., 2025. A model to identify causality for geographic patterns. International Journal of Geographical Information Science 1–21.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg",package="spEDM"))
#' \donttest{
#' spEDM::gpc(columbus,"crime","hoval",E = 7,k = 9)
#'
#' # convergence diagnostics
#' g = spEDM::gpc(columbus,"crime","hoval",libsizes = seq(7,49,7),E = 7,k = 9)
#' plot(g)
#' }
methods::setMethod("gpc", "sf", .gpc_sf_method)

#' @rdname gpc
methods::setMethod("gpc", "SpatRaster", .gpc_spatraster_method)
