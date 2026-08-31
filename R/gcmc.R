.gcmc_sf_method = \(data, cause, effect, libsizes = NULL, E = 3, k = min(E^2), tau = 1, style = 1, lib = NULL, pred = NULL, dist.metric = "L2",
                    threads = detectThreads(), detrend = FALSE, parallel.level = "low", bidirectional = TRUE, progressbar = TRUE, nb = NULL){
  varname = .check_character(cause, effect)
  E = .check_inputelementnum(E,4)
  tau = .check_inputelementnum(tau,4)
  k = .check_inputelementnum(k,2)
  pl = .check_parallellevel(parallel.level)
  .varname = .internal_varname()
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  coords = as.data.frame(sdsfun::sf_coordinates(data))
  data = sf::st_drop_geometry(data)
  data = data[,varname]
  names(data) = .varname

  if (detrend){
    data = .internal_detrend(data,.varname,coords)
  }
  cause = data[,"cause",drop = TRUE]
  effect = data[,"effect",drop = TRUE]

  if (is.null(lib)) lib = .internal_library(data)
  if (is.null(pred)) pred = lib
  if (is.null(libsizes)) libsizes = length(lib)

  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppGCMC4Lattice(cause,effect,nb,libsizes,lib,pred,E,tau,k[1],0,style,
                                .check_distmetric(dist.metric),threads,pl,progressbar)
  }
  y_xmap_x = RcppGCMC4Lattice(effect,cause,nb,libsizes,lib,pred,rev(E),rev(tau),k[2],0,
                              style,.check_distmetric(dist.metric),threads,pl,progressbar)

  return(.bind_intersectdf(varname,x_xmap_y,y_xmap_x,bidirectional))
}

.gcmc_spatraster_method = \(data, cause, effect, libsizes = NULL, E = 3, k = min(E^2), tau = 1, style = 1, lib = NULL, pred = NULL, dist.metric = "L2",
                            threads = detectThreads(), detrend = FALSE, parallel.level = "low", bidirectional = TRUE, progressbar = TRUE, grid.coord = TRUE, embed.direction = 0){
  varname = .check_character(cause, effect)
  E = .check_inputelementnum(E,4)
  tau = .check_inputelementnum(tau,4)
  k = .check_inputelementnum(k,2)
  pl = .check_parallellevel(parallel.level)
  .varname = .internal_varname()
  data = data[[varname]]
  names(data) = .varname

  dtf = .internal_grid2df(data,grid.coord)
  if (detrend){
    dtf = .internal_detrend(dtf,.varname)
  }
  causemat = matrix(dtf[,"cause"],nrow = terra::nrow(data),byrow = TRUE)
  effectmat = matrix(dtf[,"effect"],nrow = terra::nrow(data),byrow = TRUE)

  if (is.null(lib)) lib = .internal_library(dtf,TRUE)
  if (is.null(pred)) pred = lib
  if (is.null(libsizes)) libsizes = matrix(nrow(lib))

  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppGCMC4Grid(causemat,effectmat,libsizes,lib,pred,E,tau,k[1],0,style,
                             .check_distmetric(dist.metric),threads,pl,progressbar,embed.direction)
  }
  y_xmap_x = RcppGCMC4Grid(effectmat,causemat,libsizes,lib,pred,rev(E),rev(tau),k[2],0,style,
                           .check_distmetric(dist.metric),threads,pl,progressbar,embed.direction)

  return(.bind_intersectdf(varname,x_xmap_y,y_xmap_x,bidirectional))
}

#' geographical cross mapping cardinality
#'
#' @param data observation data.
#' @param cause name of causal variable.
#' @param effect name of effect variable.
#' @param libsizes (optional) number of spatial units used (input needed: `vector` - spatial vector, `matrix` - spatial raster).
#' @param E (optional) embedding dimensions.
#' @param k (optional) number of nearest neighbors.
#' @param tau (optional) step of spatial lags.
#' @param style (optional) embedding style (`0` includes current state, `1` excludes it).
#' @param lib (optional) libraries indices (input requirement same as `libsizes`).
#' @param pred (optional) predictions indices (input requirement same as `libsizes`).
#' @param dist.metric (optional) distance metric (`L1`: Manhattan, `L2`: Euclidean).
#' @param threads (optional) number of threads to use.
#' @param detrend (optional) whether to remove the linear trend.
#' @param parallel.level (optional) level of parallelism, `low` or `high`.
#' @param bidirectional (optional) whether to examine bidirectional causality.
#' @param progressbar (optional) whether to show the progress bar.
#' @param nb (optional) neighbours list.
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{cross mapping results}
#' \item{\code{cs}}{causal strength}
#' \item{\code{varname}}{names of causal and effect variables}
#' \item{\code{bidirectional}}{whether to examine bidirectional causality}
#' }
#' @export
#' @name gcmc
#' @aliases gcmc,sf-method
#' @references
#' Lyu, W., Dai, S., Song, Y., Zhao, W., Yi, W., Xiao, Y., Jia, N., 2026. Measuring causal strengths from spatial cross-sectional data with geographical cross mapping cardinality. International Journal of Geographical Information Science 1–23.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg",package="spEDM"))
#' \donttest{
#' g = spEDM::gcmc(columbus,"hoval","crime",E = 7,k = 19)
#' g
#' }
methods::setMethod("gcmc", "sf", .gcmc_sf_method)

#' @rdname gcmc
#' @param grid.coord (optional) whether to detrend using cell center coordinates (`TRUE`) or row/column numbers (`FALSE`).
#' @param embed.direction (optional) direction selector for embeddings (`0` returns all directions, `1-8` correspond to NW, N, NE, W, E, SW, S, SE).
methods::setMethod("gcmc", "SpatRaster", .gcmc_spatraster_method)
