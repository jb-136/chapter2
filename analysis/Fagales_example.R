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


source("~//chapter2//R//merging_trees_with_MRP.R")
mergedTree <- merging_trees_with_MRP(
	tree_backbone, tree_secondary, 
	reduce_collapse = TRUE,
	trace=1)
mergedTree


pdf_path <- paste0("~//chapter2//analysis/",
		"merged_Fagales_tree_",
		format(Sys.time(), "%m-%d-%y"),
		".pdf")


pdf(file = pdf_path,
	height = 100,
	width = 4)
plot(mergedTree,
	cex=0.3,
	show.tip.label = TRUE,
	no.margin = TRUE)
dev.off()
# need to use path expand to convert ~
shell.exec(path.expand(pdf_path))


########################################################


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
