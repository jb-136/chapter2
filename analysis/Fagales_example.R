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
	ott_ids = dq$ott_ids, resolve = FALSE)
is.binary(tree_backbone)
# and with Fagales PBDB taxon tree
faData <- paleotree::getCladeTaxaPBDB("Fagales")
# make the taxon tree
tree_secondary <- paleotree::makePBDBtaxonTree(
     taxaDataPBDB = faData,
     rankTaxon = "species",
     method = "parentChild"
     )
# look at stats of the trees
Ntip(tree_backbone)
Nnode(tree_backbone)
Ntip(tree_secondary)
Nnode(tree_secondary)

source("merging_trees_with_MRP.R")
source("~/chapter2/R/merging_trees_with_MRP.R")
mergedTrees <- merging_trees_with_MRP(
	tree_backbone, tree_secondary, trace=1)

# then get the strict consensus?

strictMerged <- consensus(mergedTrees, p = 1, check.labels = TRUE)
plot(strictMerged)

pdf("fagales_mrp_fan.pdf", height = 10)
mergedTreesx <- ape::compute.brlen(mergedTrees)
plot(mergedTreesx, cex = 0.2, type = "fan")
dev.off()

pdf("fagales_mrp.pdf", height = 50)
plot(mergedTrees, cex = 0.2)
dev.off()
