#' Do chapter2 analysis
#'
#' @param taxon Clade of interest
#' @return A list with all output
run_chapter2 <- function(taxon) {
  datelife_biggest <- get_datelife_biggest(taxon)
  pubmed <- get_pubmed(taxon)
  genbank_count_by_gene <- get_genbank_count_by_gene(taxon)
  genbank_count <- get_genbank_count(taxon)
  otol <- get_otol(taxon)
  return(list(datelife_biggest=datelife_biggest, pubmed=pubmed, genbank_count_by_gene=genbank_count_by_gene, genbank_count=genbank_count , otol=otol))

}

#' Run biogeobears analyses
#'
#' Run a set of biogeobears analyses to estimate parameters of biogeographic models
#'
#' @param phy The phylo object
#' @param locations The states for biogeographic regions
#' @return A list with fits for each model
#' @export
run_biogeobears <- function(phy, locations) {

}
