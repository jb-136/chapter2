#' Get biogeographic data
#'
#' This will pull information from the specified service, clean it, and put it
#' in a uniform format. There is an option to clean the data to eliminate
#' potentially problematic points.
#'
#' @param species The species to get information about. Can be a vector.
#' @param service One of paleobiodb, gbif
#' @param clean Boolean; if TRUE, remove problematic points
#' @return A data.frame with a column for species, service, latitude, longitude, and altitude
#' @export
species_get_biogeographic_data <- function(species, service="gbif", clean=TRUE) {

}
