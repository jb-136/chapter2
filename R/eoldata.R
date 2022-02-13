#' Query Encyclopedia of Life with a URL to return all data for a given taxon
#'
#' @param url URL to query, include /data 
#' @return A data frame of trait category, value, source, and definitions
#' @export
 
eol_data <- function(url="https://eol.org/pages/491995/data") {
	#library(rvest)
	#url <- "https://eol.org/pages/919224/data"
	#url <- "https://eol.org/pages/491995/data"
	input <-  rvest::read_html(url)
	all_ul <- input %>% html_elements('ul')
	trait_ul <- all_ul[[5]]
	trait_list_text <- trait_ul %>% html_nodes("div") %>% html_text2()
	trait_list_text <- gsub('\\d* record hidden \\— show all', "", trait_list_text)
	trait_list_text <- gsub('\nshow all records', "", trait_list_text)
	trait_list_raw <- trait_ul %>% html_nodes("div") %>% as.character()

	empty <- which(nchar(trait_list_text)==0)
	trait_list_text <- trait_list_text[-empty]
	trait_list_raw <- trait_list_raw[-empty]
	data_heads <- which(grepl("h3", trait_list_raw))
	trait_df <- data.frame(matrix(nrow=0, ncol=2))
	colnames(trait_df) <- c("trait", "value", "source", "URI", "definition")
	data_head_plus_end <- c(-1+data_heads[-1], length(trait_list_text))
	for (i in seq_along(data_heads)) {
	relevant_rows <- trait_list_text[(data_heads[i]+1):(data_head_plus_end[i])]
	
	relevant_rows <- relevant_rows[grepl(".+\\\n.+\\\nURI", relevant_rows)]
	
	for(j in seq_along(relevant_rows)) {
		trait_info <- strsplit(relevant_rows[j], "\n")[[1]][2]
		source_info <- strsplit(relevant_rows[j], "\n")[[1]][1]
		URI_info <- strsplit(relevant_rows[j], "\n")[[1]][3]
		definition_info <- NA
		try(definition_info <- strsplit(relevant_rows[j], "\n")[[1]][4]
	)

		trait_df <- rbind(trait_df, data.frame(trait=gsub("\n", "", trait_list_text[data_heads[i]]), value=trait_info, source=source_info, URI=URI_info, definition=definition_info))
	}
	}
	return(trait_df)
}