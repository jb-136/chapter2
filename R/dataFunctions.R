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

remotes::install_github("ropensci/rfishbase")
library("rfishbase")
library("dplyr")

dat <- distribution(species_list(Genus='Labroides'))
dat <- dat[c("SpecCode", "Species", "NorthernLatitude", "NorthernLatitudeNS", "SouthernLatitude", "SouthernLatitudeNS", "WesternLongitude", "WesternLongitudeEW", "EasternLongitude", "EasternLongitudeEW")]


#' Get information on fish localities rfishbase
#'
#' @param genus fish genus name
#' @return a data frame of SpecCode, Species, NorthernLatitude, NorthernLatitudeNS, SouthernLatitude, SouthernLatitudeNS, WesternLongitude, WesternLongitudeEW, EasternLongitude, EasternLongitudeEW
#' @export

fishbase_genus_query <- function (genus) {
  dat <- rfishbase::distribution(species_list(Genus=genus))
  dat <- dat[c("SpecCode", "Species", "NorthernLatitude", "NorthernLatitudeNS", "SouthernLatitude", "SouthernLatitudeNS", "WesternLongitude", "WesternLongitudeEW", "EasternLongitude", "EasternLongitudeEW")]

  return(dat)
}

#' Get information on fish localities rfishbase
#'
#' @param species vector of fish species names
#' @return a data frame of SpecCode, Species, NorthernLatitude, NorthernLatitudeNS, SouthernLatitude, SouthernLatitudeNS, WesternLongitude, WesternLongitudeEW, EasternLongitude, EasternLongitudeEW
#' @export

fishbase_species_query <- function (species) {
  all.records <- data.frame()
  for (species_index in seq_along(species)) {
    all.records <- plyr::rbind.fill(all.records, rfishbase::distribution(species_list=species[species_index]))
  }
  return(all.records)
}

#' Get info from Open Tree of Life
#'
#' @param taxon The clade to investigate
#' @return data.frame of info about that taxon. Columns include year of the study, number of taxa in the tree from that study, the study citation, and the study DOI.
#' @export
#' @examples
#' info <- get_otol("Gallus")
#' histogram(info$year) # Years in which chicken papers in OpenTree were published
get_otol <- function(taxon) {
  clade.info <- rotl::tnrs_match_names(taxon)
  clade.name <- clade.info$unique_name[1]
  id <- clade.info$ott_id[1]
  node.info <- rotl::tol_node_info(id)
  relevant.studies <- rotl::studies_find_trees(property="ot:ottTaxonName", value=clade.name)
  tree.info <- data.frame()
  all.trees <- rotl::list_trees(relevant.studies)
  for (study.index in sequence(nrow(relevant.studies))) {
    study.info <- rotl::get_publication(rotl::get_study_meta(relevant.studies$study_ids[study.index]))
    for (tree.index in sequence(length(all.trees[study.index]))) {
      phy <- NULL
      try(phy <- rotl::get_study_tree(study_id=relevant.studies$study_ids[study.index], tree_id = gsub('tree/', '',all.trees[[study.index]][tree.index])))
      local.result <- NULL
      if(!is.null(phy)) {
        try(local.result <- data.frame(Year=relevant.studies$study_year[study.index], Ntax=ape::Ntip(phy), Pub=study.info[1], DOI=attr(study.info, "DOI")))
      } else {
        try(local.result <- data.frame(Year=relevant.studies$study_year[study.index], Ntax=NA, Pub=study.info[1], DOI=attr(study.info, "DOI")))
      }
      if(!is.null(local.result)) {
        if(nrow(tree.info)==0) {
          tree.info <- local.result
        } else {
          tree.info <- rbind(tree.info, local.result)
        }
      }

    }
  }
  return(tree.info)
}

#' Information about the number of species in the clade in genbank
#'
#' @param taxon The clade of interest
#' @return The count of species in genbank
#' @export
#' @examples
#' taxon <- "Myrmecocystus"
#' print(paste("There are", get_genbank_count(taxon), "species of", taxon, "in GenBank"))
get_genbank_count <- function(taxon) {
  clade.name<- rotl::tnrs_match_names(taxon)$unique_name[1]
  genbank.species.query <- paste0(clade.name, '[subtree] AND species[rank] AND specified[prop]')
  genbank.species.count <-  rentrez::entrez_search(db="taxonomy", genbank.species.query, use_history=TRUE)$count
}

#' Get the count of number of sequences for different genes for this taxon
#'
#' This will use a few popular genes by default, but you can pass your own instead.
#'
#' @param taxon The clade of interest
#' @param focal.genes Vector of gene names
#' @return vector of counts by gene, with label being the given gene name
#' @examples
#' get_genbank_count_by_gene("Myrmecocystus")
get_genbank_count_by_gene <- function(taxon, focal.genes=c("COI", "18S", "28S", "matk", "rbcl")) {
  clade.name<- rotl::tnrs_match_names(taxon)$unique_name[1]
  GetNucCount <- function(gene, taxon=clade.name) {
    gene.query <- paste0(taxon, '[organism] AND ',gene)
    Sys.sleep(3) #just to make sure we stay nice
    return(rentrez::entrez_search(db="nuccore", gene.query, use_history=TRUE)$count)
  }
  gene.count <- sapply(focal.genes, GetNucCount, taxon=clade.name) #make sure not to do mclapply or similar lest you
  return(gene.count)
}

#' Get information on the clade from pubmed
#'
#' This will give the total number of pubmed articles mentioning the taxon name and other search string and information on the most recent such papers.
#'
#' @param taxon The clade of interest
#' @param search.string Include spaces and AND and similar search elements
#' @param retmax how many papers to return
#' @return List with a count element (integer) and a recent.papers element (data.frame)
#' @export
#' @examples
#' taxon <- "Formicidae"
#' results <- get_pubmed(taxon)
#' print(paste("There are", results$count, "papers on", taxon, "and phylogeny"))
#' print(results$recent.papers[,c("Date", "FirstAuthor")])

get_pubmed <- function(taxon, search.string=' AND phylogeny',retmax=50) {
  clade.name<- rotl::tnrs_match_names(taxon)$unique_name[1]
  pubmed.query <- paste0(clade.name, search.string)
  pubmed.result <- rentrez::entrez_search(db="pubmed", pubmed.query, use_history=TRUE, retmax=retmax)
  pubmed.summaries <- rentrez::entrez_summary(db="pubmed", id=pubmed.result$id)
  pubmed.df <- data.frame(Date=rentrez::extract_from_esummary(pubmed.summaries, elements=c("sortpubdate")), FirstAuthor=rentrez::extract_from_esummary(pubmed.summaries, elements=c("sortfirstauthor")), Journal=rentrez::extract_from_esummary(pubmed.summaries, elements=c("fulljournalname")), Title=rentrez::extract_from_esummary(pubmed.summaries, elements=c("title")), row.names=NULL)
  return(list(count=pubmed.result$count,   recent.papers =   pubmed.df ))
}

#' Get biggest tree from datelife
#'
#' @param taxon Clade of interest
#' @return phylo object
#' @export
get_datelife_biggest <- function(taxon) {
  clade.name<- rotl::tnrs_match_names(taxon)$unique_name[1]
  datelife_biggest <- NULL
  try(datelife_biggest <- datelife::datelife_search(input=clade.name, get_spp_from_taxon=TRUE, summary_format="phylo_biggest"))
  return(datelife_biggest)
}

#' Get Wikipedia summary
#'
#' @param taxon Clade of interest
#' @return text string of summary of page
#' @export
get_wikipedia_summary <- function(taxon) {
  URL <- paste0('https://en.wikipedia.org/w/api.php?format=json&action=query&prop=extracts&exintro&explaintext&redirects=1&titles=', utils::URLencode(taxon))
  return(jsonlite::fromJSON(URL)$query$pages[[1]]$extract)
}
