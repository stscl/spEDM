.sc_sf_method = \(data, cause, effect, k, block = 3, boot = 399, seed = 42L, base = 2, lib = NULL, pred = NULL,
                  threads = detectThreads(), detrend = TRUE, normalize = FALSE, progressbar = FALSE, nb = NULL){
  varname = .check_character(cause, effect)
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  block = RcppDivideLattice(nb,block)
  cause = .uni_lattice(data,cause,detrend)
  effect = .uni_lattice(data,effect,detrend)
  if (is.null(lib)) lib = which(!(is.na(cause) | is.na(effect)))
  if (is.null(pred)) pred = lib
  return(.bind_sc(RcppSGC4Lattice(cause,effect,nb,lib,pred,block,k,threads,boot,base,seed,TRUE,normalize,progressbar),varname))
}

.sc_spatraster_method = \(data, cause, effect, k, block = 3, boot = 399, seed = 42L, base = 2, lib = NULL, pred = NULL, threads = detectThreads(),
                          detrend = TRUE, normalize = FALSE, progressbar = FALSE, grid.coord = TRUE, embed.direction = 0){
  varname = .check_character(cause, effect)
  cause = .uni_grid(data,cause,detrend,grid.coord)
  effect = .uni_grid(data,effect,detrend,grid.coord)
  block = matrix(RcppDivideGrid(effect,block),ncol = 1)
  if (is.null(lib)) lib = which(!(is.na(cause) | is.na(effect)), arr.ind = TRUE)
  if (is.null(pred)) pred = lib
  return(.bind_sc(RcppSGC4Grid(cause,effect,embed.direction,lib,pred,block,k,threads,boot,base,seed,TRUE,normalize,progressbar),varname))
}

#' spatial causality test
#'
#' @param data observation data.
#' @param cause name of causal variable.
#' @param effect name of effect variable.
#' @param k (optional) number of nearest neighbors used in symbolization.
#' @param block (optional) number of blocks used in spatial block bootstrap.
#' @param boot (optional) number of bootstraps to perform.
#' @param seed (optional) random seed.
#' @param base (optional) logarithm base.
#' @param lib (optional) libraries indices (input needed: `vector` - spatial vector, `matrix` - spatial raster).
#' @param pred (optional) predictions indices (input requirement same as `lib`).
#' @param threads (optional) number of threads to use.
#' @param detrend (optional) whether to remove the linear trend.
#' @param normalize (optional) whether to normalize the result.
#' @param progressbar (optional) whether to show the progress bar.
#' @param nb (optional) neighbours list.
#'
#' @return A list
#' \describe{
#' \item{\code{sc}}{statistic for spatial causality}
#' \item{\code{varname}}{names of causal and effect variables}
#' }
#' @export
#' @name sc.test
#' @aliases sc.test,sf-method
#' @references
#' Herrera, M., Mur, J., Ruiz, M., 2016. Detecting causal relationships between spatial processes. Papers in Regional Science 95, 577–595.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg",package="spEDM"))
#' \donttest{
#' spEDM::sc.test(columbus,"hoval","crime",k = 15)
#' }
methods::setMethod("sc.test", "sf", .sc_sf_method)

#' @rdname sc.test
#' @param grid.coord (optional) whether to detrend using cell center coordinates (`TRUE`) or row/column numbers (`FALSE`).
methods::setMethod("sc.test", "SpatRaster", .sc_spatraster_method)
