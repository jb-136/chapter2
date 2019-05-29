library(rgbif)
library(plyr)

#' Query gbif to return latitude and longitude for members of a clade of interest
#'
#' @param query name of the clade of interest
#' @param rank the clade's rank (e.g. "genus")
#' @return A data frame (tibble) of key, scientific name, decimal latitude, and decimal longitude
#' @export

gbif_taxon_query <- function (query, rank){
  key <- rgbif::name_suggest(q=query, rank=rank)$key[1]
  dat <- rgbif::occ_search(taxonKey=key, fields="minimal", limit=300)

  return(dat)
}

# e.g. gbif_clade_query("fagales", "class")


#' Query gbif to return latitude and longitude for a vector of species
#'
#' @param species a vector of species names
#' @return A data frame (tibble) of key, scientific name, decimal latitude, and decimal longitude
#' @export

gbif_species_query <- function (species){
  all.records <- data.frame()
  for (species_index in seq_along(species)) {
    all.records <- plyr::rbind.fill(all.records, rgbif::occ_search(scientificName = species[species_index], fields="minimal", limit = 300)$data)
  }
  return(all.records)
}


# e.g.  species <- c("Puma concolor", "Myrmecocystus mexicanus", "Ursus arctos")
# gbif_species_query(species=species)


#' Get information on specimens from the paleobiology data base (largely for extinct taxa)
#'
#' @param taxon taxon name
#' @return a data frame of name, paleo longitude, paleo latitude, earliest record in ma (minma), latest record in ma (maxma), earliest interval, and latest interval
#' @export
pbdb_taxon_query <- function(taxon){
  pbdb_data <- read.csv(paste0("http://paleobiodb.org/",
    "data1.2/occs/list.txt?base_name=",utils::URLencode(taxon),"&level=3&show=paleoloc"),
    stringsAsFactors = FALSE)
  lat_long <- data.frame(pbdb_data$accepted_name, pbdb_data$paleolng, pbdb_data$paleolat, pbdb_data$max_ma, pbdb_data$min_ma, pbdb_data$early_interval, pbdb_data$late_interval)
  lat_long$searched_taxon <- taxon
  return(lat_long)
}
