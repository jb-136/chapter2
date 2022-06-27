#' Query Encyclopedia of Life with a URL to return all data for a given taxon
#'
#' @param species Species name to query
#' @return A data frame of trait category, value, source, and definitions
#' @export
 
 
 #https://eol.org/api/search/1.0.json?q=Camelis%20dromedarius&exact=1&page=1
 #https://eol.org/api/search/1.0.json?q=Camelus%20dromedarius&exact=1&page=1
eol_data <- function(species) {
	#library(rvest)
	#url <- "https://eol.org/pages/919224/data"
	#url <- "https://eol.org/pages/491995/data"
	searchurl <- paste('http://eol.org/api/search/1.0.json?q=', URLencode(species), '&exact=1&page=1&key=', sep="")
	url <- NA
	url <- paste0(jsonlite::fromJSON(searchurl)$results$link[1], "/data")

	input <-  rvest::read_html(url)
	all_ul <-  rvest::html_elements(input,'ul')
	trait_ul <- all_ul[[5]]
	trait_list_text <- rvest::html_text2(rvest::html_nodes(trait_ul, "div"))
  trait_list_text <- gsub("([0-9]*) records hidden", " \\1 records hidden", trait_list_text)
  trait_list_text <- gsub("([0-9]*) record hidden", " \\1 record hidden", trait_list_text)
	trait_list_text <- gsub('\\d* record hidden \\— show all', "", trait_list_text)
	trait_list_text <- gsub('\\d* records hidden \\— show all', "", trait_list_text)

	trait_list_text <- gsub('\nshow all records', "", trait_list_text)
	trait_list_raw <- as.character(rvest::html_nodes(trait_ul, "div"))

	empty <- which(nchar(trait_list_text)==0)
	trait_list_text <- trait_list_text[-empty]
	trait_list_raw <- trait_list_raw[-empty]
	data_heads <- which(grepl("h3", trait_list_raw))
	trait_df <- data.frame(matrix(nrow=0, ncol=6))
	colnames(trait_df) <- c("species", "trait", "value", "source", "URI", "definition")
	data_head_plus_end <- c(-1+data_heads[-1], length(trait_list_text))
	for (i in seq_along(data_heads)) {
		relevant_rows <- trait_list_text[(data_heads[i]+1):(data_head_plus_end[i])]
		
		relevant_rows <- relevant_rows[grepl(".+\\\n.+\\\nURI", relevant_rows)]
		
		for(j in seq_along(relevant_rows)) {
			trait_info <- strsplit(relevant_rows[j], "\n")[[1]][2]
			source_info <- strsplit(relevant_rows[j], "\n")[[1]][1]
			URI_info <- strsplit(relevant_rows[j], "\n")[[1]][3]
			definition_info <- NA
			try(definition_info <- strsplit(relevant_rows[j], "\n")[[1]][4])

			trait_df <- rbind(trait_df, data.frame(species=species, trait=gsub("\n", "", trait_list_text[data_heads[i]]), value=trait_info, source=source_info, URI=URI_info, definition=definition_info))
		}
	}
	return(trait_df)
}
