#' Continuous to discrete
#'
#' Convert a set of latitudes and longitudes to biogeographic regions
#'
#' @param latlon The data.frame of points, with species, latitude, and longitude columns
#' @return A data.frame as above but with a region column added
#' @export
convert_continuous_to_discrete <- function(latlon) {
}


#' Get a tree and data, resolve to same set of taxa, and make a shared object
#'
#' @param tree a phylo object
#' @param data A data.frame. Rownames are species names.
#' @return The function returns a list of two elements (phy and data) that are manipulated to include only those species found in both the tree and data supplied by the user. The class of the returned object is "chapter2".
#' @export
match_data <- function(tree, data){
  matchedData <- geiger::treedata(phy= tree, data = data, sort=TRUE) #Use geiger to merge trait with phylogeny and drop out unmatched taxa
  class(match_data) <- c("list","chapter2")
  return(match_data)
}

#' Convert the chapter2 format to one usable for Beaulieu-style
#'
#' Many of the packages authored by Jeremy Beaulieu with collaborators have the first column used to hold taxon names rather than store them as row names. Bless his heart. This converts a chapter2 object, which has its data with taxon names as row names, into having a data object internally with a column for taxon names.
#'
#' @param chapter2 The chapter2 object
#' @return a list of two objects, phy and data, with the data in the format preferred by Beaulieu
chapter2_convert_to_Beaulieu_data <- function(chapter2) {
  Beaulieu_Convert <- cbind.data.frame(rownames(match_data$data, match_data$data[,1:ncol(match_data$data)]))
}
