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
#' @param warnings Complain about mismatches if TRUE.
#' @return The function returns a list of two elements (phy and data) that are manipulated to include only those species found in both the tree and data supplied by the user. The class of the returned object is "chapter2".
#' @export
match_data <- function(tree, data, warnings = FALSE) {
  dm = length(dim(data))
  if (is.vector(data)) {
    data <- as.matrix(data)
  }
  if (is.factor(data)) {
    data <- as.matrix(data)
  }
  if (is.array(data) & length(dim(data)) == 1) {
    data <- as.matrix(data)
  }
  if (is.null(rownames(data))) {
    stop("names for 'data' must be supplied")
  }
  else {
    data.names <- rownames(data)
  }
  nc <- geiger::name.check(phy, data)
  if (is.na(nc[[1]][1]) | nc[[1]][1] != "OK") {
    if (length(nc[[1]] != 0)) {
      phy = ape::drop.tip(phy, as.character(nc[[1]]))
      if (warnings) {
        warning(paste("The following tips were not found in 'data' and were dropped from 'phy':\n\t",
                      paste(nc[[1]], collapse = "\n\t"), sep = ""))
      }
    }
    if (length(nc[[2]] != 0)) {
      m <- match(data.names, nc[[2]])
      data = as.data.frame(data[is.na(m), ],stringsAsFactors=FALSE)
      data.names <- data.names[is.na(m)]
      if (warnings) {
        warning(paste("The following tips were not found in 'phy' and were dropped from 'data':\n\t",
                      paste(nc[[2]], collapse = "\n\t"), sep = ""))
      }
    }
  }
  order <- match(data.names, phy$tip.label)
  rownames(data) <- phy$tip.label[order]

  index <- match(phy$tip.label, rownames(data))
  data <- as.data.frame(data[index, ], stringsAsFactors=FALSE)
  if (dm == 2) {
    data <- as.data.frame(data, stringsAsFactors=FALSE)
  }
  phy$node.label = NULL
  MatchReturn<-list(phy = phy, data = data)
  class(MatchReturn)<-c("list","chapter2")
  return(MatchReturn)
}


#' Convert the chapter2 format to one usable for Beaulieu-style
#'
#' Many of the packages authored by Jeremy Beaulieu with collaborators have the first column used to hold taxon names rather than store them as row names. Bless his heart. This converts a chapter2 object, which has its data with taxon names as row names, into having a data object internally with a column for taxon names.
#'
#' @param chapter2 The chapter2 object
#' @return a list of two objects, phy and data, with the data in the format preferred by Beaulieu
chapter2_convert_to_Beaulieu_data <- function(chapter2) {
  return(list(phy=chapter2$phy, data=cbind.data.frame(taxa=rownames(match_data$data), chapter2$data[,1:ncol(chapter2$data)])))
}

check_numeric_vector <- function(x) {
  x <- x[!is.na(x)]
  return(all(grepl("^[0-9\\.]{1,}$",x)))
}

check_continuous <- function(x) {
  continuous <- sapply(x, check_numeric)
  numeric_cols <- which(continuous)
  for (i in seq_along(numeric_cols)) {
    local_vector <- x[,numeric_cols[i]]
    local_vector <- local_vector[!is.na(local_vector)]
    local_vector <- as.numeric(local_vector)
    if(all(local_vector==round(local_vector))) {
      continuous[numeric_cols[i]] <- FALSE
    }
  }
  return(continuous)
}


chapter2_drop_type <- function(chapter2, keep=c("continuous", "discrete")) {
  keep_continuous <- TRUE
  if(keep=="discrete") {
    keep_continuous <- FALSE
  }

  column_checks <- check_continuous(chapter2$data)
  if(!keep_continuous) {
    column_check <- !column_check
  }
  chapter2$data  <- chapter2$data[, column_check, drop=FALSE]
  return(chapter2)
}

chapter2_fitContinuous <- function(chapter2, models=c("BM","OU","EB","trend","lambda","kappa","delta","drift","white"), ncores=NULL){
  ContDat <- chapter2_drop_type(chapter2, keep="continuous")
  fitContinuousResList <- list()
  for(model_index in seq_along(models)){
    for (char_index in sequence(ncol(chapter2$data))) {
      sliced_data <- chapter2
      sliced_data$data <- sliced_data$data[,char_index, drop=TRUE]
      sliced_data$data <- sliced_data$data[!is.na(sliced_data$data)]
      sliced_data <- geiger::treedata(sliced_data$phy, sliced_data$data, sort=TRUE, warnings=FALSE)
      fitContinuousResList[[models[model_index]]][[char_index]] <- geiger::fitContinuous(phy = sliced_data$phy, data = sliced_data$data, model = models[model_index], ncores=ncores)
    }
  }
  return(fitContinuousResList)
}
