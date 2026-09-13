.multiview_sf_method = \(data, column, target, nvar, E = 3, k = E+2, tau = 1, style = 1, stack = FALSE, lib = NULL, pred = NULL,
                         dist.metric = "L2", dist.average = TRUE, top = NULL, threads = detectThreads(), detrend = TRUE, nb = NULL){
  xmat = .multivar_lattice(data,column,detrend)
  yvec = .uni_lattice(data,target,detrend)
  if (is.null(lib)) lib = .internal_library(cbind(xmat,yvec))
  if (is.null(pred)) pred = lib
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  if (is.null(top)) top = 0
  res = RcppMultiView4Lattice(xmat,yvec,nb,lib,pred,E,tau,k,top,nvar,style,stack,
                              .check_distmetric(dist.metric), dist.average, threads)
  return(res)
}

.multiview_spatraster_method = \(data, column, target, nvar, E = 3, k = E+2, tau = 1, style = 1, stack = FALSE, lib = NULL, pred = NULL, dist.metric = "L2", 
                                 dist.average = TRUE, top = NULL, threads = detectThreads(), detrend = TRUE, grid.coord = TRUE, embed.direction = 0){
  xmat = .multivar_grid(data,column,detrend,grid.coord)
  ymat = .multivar_grid(data,target,detrend,grid.coord)
  if (is.null(lib)) lib = .internal_library(cbind(xmat,ymat),TRUE)
  if (is.null(pred)) pred = lib
  if (is.null(top)) top = 0
  res = RcppMultiView4Grid(xmat, ymat, lib, pred, E, tau, k, top, nvar, style, stack,
                           .check_distmetric(dist.metric), dist.average, embed.direction, threads)
  return(res)
}

#' multiview embedding forecast
#'
#' @inheritParams simplex
#' @param nvar number of variable combinations.
#' @param top (optional) number of reconstructions used in MVE forecast.
#'
#' @return A vector (when input is sf object) or matrix
#' @export
#' @name multiview
#' @aliases multiview,sf-method
#' @references
#' Ye H., and G. Sugihara, 2016. Information leverage in interconnected ecosystems: Overcoming the curse of dimensionality. Science 353:922-925.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg",package="spEDM"))
#' \donttest{
#' spEDM::multiview(columbus,
#'                  column = c("inc","crime","open","plumb","discbd"),
#'                  target = "hoval", nvar = 3)
#' }
methods::setMethod("multiview", "sf", .multiview_sf_method)

#' @rdname multiview
methods::setMethod("multiview", "SpatRaster", .multiview_spatraster_method)
