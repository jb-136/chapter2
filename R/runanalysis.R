#' Do chapter2 analysis
#'
#' @param taxon Clade of interest
#' @return A list with all output
#' @export
run_chapter2 <- function(taxon) {
  wikipedia_summary <- get_wikipedia_summary(taxon)
  datelife_biggest <- get_datelife_biggest(taxon)
  pubmed <- get_pubmed(taxon)
  genbank_count_by_gene <- get_genbank_count_by_gene(taxon)
  genbank_count <- get_genbank_count(taxon)
  otol <- get_otol(taxon)
  return(list(wikipedia_summary=wikipedia_summary, datelife_biggest=datelife_biggest, pubmed=pubmed, genbank_count_by_gene=genbank_count_by_gene, genbank_count=genbank_count , otol=otol))

}

#' Create a file of results
#'
#' @param taxon Clade of interest
#' @param format Format: pdf or html
#' @param output_dir Where to put the output; by default, current working directory
#' @return Nothing, though a file is created in the current working directory
#' @export
#' @examples
#' render_chapter2("Tyto")
render_chapter2 <- function(taxon, format="pdf", output_dir=getwd()) {
  rmarkdown::render(system.file("rmd", "summary.Rmd", package="chapter2"), params=list(taxon=taxon),output_file=paste0("Report_",gsub(" ", "_",taxon), ".", format), output_dir=output_dir)
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
