# mergingTreesWithMRP analysis

devtools::install_github("phylotastic/datelife")
devtools::install_github("dwbapst/paleotree", ref="developmentBranch")
# not sure its necessary to load these packages really
# library(datelife)
# library(paleotree)
# probably necessary to load ape though
library(ape)

# to test with otol Fagales tree
dq <- datelife::make_datelife_query("fagales", 
	get_spp_from_taxon = TRUE)
tree_backbone<- datelife::get_otol_synthetic_tree(dq$cleaned_names,
	ott_ids = dq$ott_ids)
# and with Fagales PBDB taxon tree
faData <- paleotree::getCladeTaxaPBDB("Fagales")
# make the taxon tree
tree_secondary <- paleotree::makePBDBtaxonTree(
     taxaDataPBDB = faData,
     rankTaxon = "species",
     method = "parentChild"
     )

source("merging_trees_with_MRP.R")
mergedTree <- merging_trees_with_MRP(tree_backbone, tree_secondary)

# then get the strict consensus?