.scpcm_sf_method = \(data, cause, effect, conds, libsizes = NULL, E = 3, k = E+2, tau = 1, style = 1, stack = FALSE, lib = NULL, pred = NULL, dist.metric = "L2", dist.average = TRUE,
                     theta = 1, algorithm = "simplex", threads = detectThreads(), detrend = TRUE, parallel.level = "low", bidirectional = TRUE, progressbar = TRUE, cumulate = FALSE, nb = NULL){
  varname = .check_character(c(cause, effect, rev(conds)))
  E = .check_inputelementnum(E,length(varname)+length(conds),2*length(conds))
  tau = .check_inputelementnum(tau,length(varname)+length(conds),2*length(conds))
  k = .check_inputelementnum(k,length(varname)+length(conds),2*length(conds))
  pl = .check_parallellevel(parallel.level)
  .varname = .internal_varname(conds)
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
  condmat = as.matrix(data[,.varname[-c(1,2)],drop = FALSE])

  if (is.null(lib)) lib = .internal_library(data)
  if (is.null(pred)) pred = lib
  if (is.null(libsizes)) libsizes = length(lib)

  simplex = ifelse(algorithm == "simplex", TRUE, FALSE)
  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppSCPCM4Lattice(cause,effect,condmat,nb,libsizes,lib,pred,E[-c(2,2+seq_along(conds))],tau[-c(2,2+seq_along(conds))],k[-c(2,2+seq_along(conds))],
                                 simplex, theta, threads, pl, cumulate, style, stack, .check_distmetric(dist.metric), dist.average, TRUE, progressbar)
  }
  y_xmap_x = RcppSCPCM4Lattice(effect,cause,condmat,nb,libsizes,lib,pred,E[-1],tau[-1],k[-1],simplex,theta,threads,pl,
                               cumulate, style, stack, .check_distmetric(dist.metric), dist.average, TRUE, progressbar)

  return(.bind_xmapdf2(varname,x_xmap_y,y_xmap_x,bidirectional))
}

.scpcm_spatraster_method = \(data, cause, effect, conds, libsizes = NULL, E = 3, k = E+2, tau = 1, style = 1, stack = FALSE, lib = NULL, pred = NULL, dist.metric = "L2", dist.average = TRUE, theta = 1, algorithm = "simplex",
                             threads = detectThreads(), detrend = TRUE, parallel.level = "low", bidirectional = TRUE, progressbar = TRUE, cumulate = FALSE, grid.coord = TRUE, embed.direction = 0, win.ratio = 0){
  varname = .check_character(cause, effect, rev(conds))
  E = .check_inputelementnum(E,length(varname)+length(conds),2*length(conds))
  tau = .check_inputelementnum(tau,length(varname)+length(conds),2*length(conds))
  k = .check_inputelementnum(k,length(varname)+length(conds),2*length(conds))
  win.ratio = .check_inputelementnum(win.ratio,2)
  pl = .check_parallellevel(parallel.level)
  .varname = .internal_varname(conds)
  data = data[[varname]]
  names(data) = .varname

  dtf = .internal_grid2df(data,grid.coord)
  if (detrend){
    dtf = .internal_detrend(dtf,.varname)
  }
  causemat = matrix(dtf[,"cause"],nrow = terra::nrow(data),byrow = TRUE)
  effectmat = matrix(dtf[,"effect"],nrow = terra::nrow(data),byrow = TRUE)
  condmat = as.matrix(dtf[,.varname[-c(1,2)],drop = FALSE])

  if (is.null(lib)) lib = .internal_library(dtf,TRUE)
  if (is.null(pred)) pred = lib
  if (is.null(libsizes)) libsizes = matrix(nrow(lib))

  simplex = ifelse(algorithm == "simplex", TRUE, FALSE)
  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppSCPCM4Grid(causemat,effectmat,condmat,libsizes,lib,pred,E[-c(2,2+seq_along(conds))],tau[-c(2,2+seq_along(conds))],k[-c(2,2+seq_along(conds))],simplex,
                              theta, threads, pl, cumulate, style, stack, .check_distmetric(dist.metric), dist.average, TRUE, progressbar, embed.direction, win.ratio)
  }
  y_xmap_x = RcppSCPCM4Grid(effectmat,causemat,condmat,libsizes,lib,pred,E[-1],tau[-1],k[-1],simplex,theta,threads,pl,cumulate,style,
                            stack, .check_distmetric(dist.metric), dist.average, TRUE, progressbar, embed.direction, win.ratio)

  return(.bind_xmapdf2(varname,x_xmap_y,y_xmap_x,bidirectional))
}

#' spatially convergent partial cross mapping
#'
#' @inheritParams gccm
#' @param conds name of conditioning variables.
#' @param cumulate (optional) serial or cumulative computation of partial cross mapping.
#'
#' @return A list
#' \describe{
#' \item{\code{pxmap}}{partial cross mapping results}
#' \item{\code{xmap}}{cross mapping results}
#' \item{\code{varname}}{names of causal, effect and conditioning variables}
#' \item{\code{bidirectional}}{whether to examine bidirectional causality}
#' }
#' @export
#' @name scpcm
#' @aliases scpcm,sf-method
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg",package="spEDM"))
#' \donttest{
#' g = spEDM::scpcm(columbus,"hoval","crime","inc",libsizes = seq(5,45,5),E = 6)
#' g
#' plot(g,ylimits = c(-0.1,0.4),ybreaks = seq(-0.1,0.4,0.1))
#' }
methods::setMethod("scpcm", "sf", .scpcm_sf_method)

#' @rdname scpcm
methods::setMethod("scpcm", "SpatRaster", .scpcm_spatraster_method)
