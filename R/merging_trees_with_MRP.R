# merging_trees_with_MRP 05-29-19
# GPL v3
# Authors: David Bapst and Luna Luisa Sanchez Reyes

merging_trees_with_MRP <- function(
		tree_backbone, tree_secondary,
		backbone_reweighting = 1,
		reduce_collapse = TRUE,
		trace = 0){
	###############################
	# add an artificial outgroup to both trees
		# this also removes any edge lengths
	tree_backbone <- add_single_taxon_to_tree(tree = tree_backbone, 
		new_tip_label = "placeholder_artificial_outgroup")
	tree_secondary <- add_single_taxon_to_tree(tree = tree_secondary, 
		new_tip_label = "placeholder_artificial_outgroup") 	
	#############
	# Make sure the trees have properly structured clade labels!! No missing!
	#
	# append a vector of NAs if either tree lacks node labels
	# cannot have trees lacking $node.label
	# we can use datelife::tre_add_nodelabels
		# it only adds names to unnamed nodes
		# (using a prefix and a consecutive number)
		# and generates a vector of names if there is no $node.label
		#
		# tree_backbone<- datelife::tree_add_nodelabels(tree_backbone, node_prefix = "unnamed_mrp_backbone")
		# tree_secondary<- datelife::tree_add_nodelabels(tree_secondary, node_prefix = "unnamed_mrp_secondary")
		#
		# but what you do here works exactly the same
	if(is.null(tree_backbone$node.label)) {
		tree_backbone$node.label <- rep(NA, Nnode(tree_backbone))
		}
	if(is.null(tree_secondary$node.label)) {
		tree_secondary$node.label <- rep(NA, Nnode(tree_secondary))
		}
	#
	# rename all unnamed clades in both
	# otherwise we're gonna run into major issues
	# count unnamed clades in tree_backbone
	unnamed_backbone <- is.na(tree_backbone$node.label) | tree_backbone$node.label == ""
	unnamed_secondary <-  is.na(tree_secondary$node.label) | tree_secondary$node.label == ""
	total_unnamed_nodes <- sum(unnamed_backbone) + sum(unnamed_secondary)
	#
	# give a number to each unnamed node
	# labeled with 'unnamed_mrp_' just in case numbers somehow ended up as clade labels
	tree_backbone$node.label[unnamed_backbone] <- paste0(
		"unnamed_mrp_", 1:sum(unnamed_backbone)
		)
	tree_secondary$node.label[unnamed_secondary] <- paste0(
		"unnamed_mrp_", (sum(unnamed_backbone)+1):total_unnamed_nodes
		)
	################
	# first, tree_backbone
			# for MRP we want a table
			# and so we can make the matrix from tree A (modern)
			# each node is a column
			# values should be 0/1
	#
	# easiest way to get parent-child info is always prop.part
	#cool
	tree_bb_proppart <- ape::prop.part(tree_backbone)
	mrp_backbone <- sapply(tree_bb_proppart, function(x){
		   member <- rep(0, Ntip(tree_backbone))
		   member[x] <- 1
		   return(member)
	  }
		)
	# relabel each column with the node label, if present
	colnames(mrp_backbone) <- tree_backbone$node.label
	# and label rows with taxon labels
	rownames(mrp_backbone) <- tree_backbone$tip.label
	#
	##########################
	# then take tree B, and code them as unknown if there is no matching label
	# code them 0/1 for matching labels
	# and code them as 0/1 for any additional nodes that don't have matching labels
	# while coding all taxa on tree A as ?
	# then we would end up with a matrix of 0,1,? values
	# with the number of rows equal to the number of unique tip labels
	# and the number of columns equal to the number of unique nodes
	#
	# use prop.part again
	tree_sec_proppart <- ape::prop.part(tree_secondary)
	mrp_sec <- sapply(tree_sec_proppart, function(x){
			member <- rep(0, Ntip(tree_secondary))
			member[x] <- 1
			return(member)
			}
		)
	# relabel each column with the node label, if present
	colnames(mrp_sec) <- tree_secondary$node.label
	rownames(mrp_sec) <- tree_secondary$tip.label
	#
	# First, identify all matching columns (nodes)
	matchingNodes <- match(colnames(mrp_sec), colnames(mrp_backbone))
	# or identify all matching OTUs first (rows)
	matchingOTUs <- match(rownames(mrp_sec), rownames(mrp_backbone))
	#
	# this will match nodes in sec that are in backbone
		# just wondering if for whatever reason tree_sec
		# has more named nodes, will this be affected?)
	# answer: no we should be fine
		# we will always get back a response as long as the sec vector
	#
	# can see this as three matrices that need to be constructed
		# in addition to the pre-existing backbone mrp matrix (top left)
			# new nodes, old OTUs (the top right quater)
			# new nodes , new OTUs (the 'diagonal', bottom right quarter)
			# old nodes, new OTUs (the bottom left quater)
	# the backbone mrp matrix won't be altered at all
		# nothing gets added or removed from that
	# get dims for all of these
	n_old_OTUs <-  nrow(mrp_backbone)
	n_old_nodes <- ncol(mrp_backbone)
	n_new_OTUs <- sum(is.na(matchingOTUs))
	n_new_nodes <- sum(is.na(matchingNodes))
	#
	# get the name vectors for rows and columns before proceeding
	# need to make sure name vectors are right
	new_mrp_rownames <- c(rownames(mrp_backbone),
	rownames(mrp_sec)[is.na(matchingOTUs)])
	new_mrp_colnames <- c(colnames(mrp_backbone),
	colnames(mrp_sec)[is.na(matchingNodes)])
	# make these new matrices
	#
	# first new nodes, old OTUs
	# make sure there are new nodes
	if(n_new_nodes > 0){
		# make fake matrices full of '?' for the old/new and new/old matrices
		newNodes_oldOTUs <- matrix('?', n_old_OTUs, n_new_nodes)
		# replace where there are matches
		newNodes_oldOTUs [matchingOTUs[!is.na(matchingOTUs)],] <- mrp_sec[
			!is.na(matchingOTUs), is.na(matchingNodes)
			]
		# combine with mrp_backbone (left and right)
		mrp_backbone <- cbind(mrp_backbone, newNodes_oldOTUs)
		}
	#
	# now old nodes, new OTUs
		# make sure there are new OTUs
	if(n_new_OTUs > 0){
		# make fake matrices full of '?' for the old/new and new/old matrices
		oldNodes_newOTUs <- matrix('?', n_new_OTUs, n_old_nodes)
		#
		# replace where there are matches
		oldNodes_newOTUs [, matchingNodes[!is.na(matchingNodes)]] <- mrp_sec[
			is.na(matchingOTUs), !is.na(matchingNodes)
			]
		# now do the new/new matrix INSIDE this if statement
		# need to make sure there are new OTUs and new nodes
		if(n_new_OTUs > 0 & n_new_nodes > 0){
			# now let's make the bottom-right quarter
				# its the matrix with no 'matches'
			newNodes_newOTUs <- mrp_sec[is.na(matchingOTUs),
				is.na(matchingNodes)]
			# combine left and right
			oldNodes_newOTUs <- cbind(oldNodes_newOTUs, newNodes_newOTUs)
			}
		# combine, top and bottom
		mrp_backbone <-rbind(mrp_backbone, oldNodes_newOTUs)
		}
	#
	mrp_full <- mrp_backbone
	#
	# replace row and col names
	rownames(mrp_full) <- new_mrp_rownames
	colnames(mrp_full) <- new_mrp_colnames
	#print(dim(mrp_full))
	#
	################################
	# reweight the representation in the backbone (so its higher than the secondary?
	# nice!
	if(backbone_reweighting > 1){
		# repeat the columns in the backbone mrp X times, so that the same node is being
			# considered as multiple characters (thus weighting those groups 'higher'
		newColumns <- rep(1:n_old_nodes,backbone_reweighting)
		newColumns <- c(newColumns, (1+n_old_nodes):ncol(mrp_full))
		mrp_full <- mrp_full[,newColumns]
		# NOTE 
		# this will repeat the old nodes many times, possibly including the membership
			# of new OTUs in those old nodes, which wasn't in the backbone
		# but if we weight the inclusion of new OTUs less
				# then they will look like they belong to a sister group 
				# or will appear paraphyletic to the taxa in the backbone tree
		}
	#####################################
	# and then we do MRP
	# or really just a basic parsimony search with phangorn
	#
	supertrees_out <- parsimony_search_clade_collapse(
		char_matrix = mrp_full,
		levels = 0:1,  ambiguity = "?", trace = trace,
		reduce_collapse = reduce_collapse
		)
	#
	# is it one tree or more?
		# root the trees based on artificial outgroup
			#and remove artificial outgroup
		# code modeled on phangorn's supertree functions
	if (inherits(supertrees_out, "multiPhylo")) {
		supertrees_out  <- lapply(supertrees_out, root,
			"placeholder_artificial_outgroup")
		supertrees_out <- lapply(supertrees_out, drop.tip,
			"placeholder_artificial_outgroup")
        class(supertrees_out) <- "multiPhylo"
      }else{
        supertrees_out <- root(supertrees_out, 
			"placeholder_artificial_outgroup")
        supertrees_out <- drop.tip(supertrees_out, 
			"placeholder_artificial_outgroup")
      }
	#
	# and voilla, you'd get a tree sample you can do
		# a strict consensus on, or whatever
	#
	return(supertrees_out)
	}
	
	
add_single_taxon_to_tree <- function(tree, 
		new_tip_label,
		# default location to add tip is the root
		nodeID = Ntip(tree) + 1
		){
	###############################
	# currently only handles trees without branch lengths
		# in fact the branch lengths will be removed from the input tree
	# remove branch lengths
	tree$edge.length <- NULL
	#
	# make an artificial 1 tip tree
	one_tip_tree <- list(
		edge = matrix(c(2,1),1,2),
		tip.label = new_tip_label,
		edge.length = NULL,
		Nnode = 1)
	# make it class phylo
	class(one_tip_tree)<-"phylo"
	#
	tree <- bind.tree(
		x = tree,
		y = one_tip_tree,
		where = nodeID
		)
	return(tree)
	}


parsimony_search_clade_collapse <- function(
		char_matrix, 
		levels,  ambiguity, trace,
		reduce_collapse = TRUE
		){
	####################
	#
	if(reduce_collapse){
		# identify sets of taxa that share the same coding
		taxonSets <- apply(char_matrix, 1,function(x)
			which(apply(char_matrix, 1,function(y)
				identical(x,y)
				))[1]
			)
		# count how many are repeated
		nRepeats <- table(taxonSets)
		if(any(nRepeats > 2)){
			# remove all but two, rename with placeholders
			collapse_these_sets <- names(nRepeats[nRepeats > 2])
			# save information on the removed rows
			saved_sets <- list()
			matrix_modified <- char_matrix
			for(i in 1:length(collapse_these_sets)){
				rows_in_set <- which(taxonSets == collapse_these_sets[i])
				set_OTU_labels <- rownames(char_matrix)[rows_in_set]
				saved_sets[[i]] <- set_OTU_labels
				# find the rows in the modified matrix
					# use which so we can modify as vector
				mod_rows_in_set <- which(sapply(rownames(matrix_modified),
					function(x) any(x == set_OTU_labels)))
				# remove all but two of the set
					# dropping all but first two will not modify location of first two
				matrix_modified <- matrix_modified[-mod_rows_in_set[-(1:2)],]
				# new tip names
				new_OTU_names <- paste0("placeholder_set_", collapse_these_sets[i])
				new_OTU_names <- paste0(new_OTU_names, "_", c("a", "b"))
				rownames(matrix_modified)[mod_rows_in_set[1:2]] <- new_OTU_names				
				}
			names(saved_sets) <- paste0("placeholder_set_",collapse_these_sets)
			matrix_final <- matrix_modified 
		}else{
			matrix_final <- char_matrix
			reduce_collapse <- FALSE
			}
	}else{
		matrix_final <- char_matrix
		}
	##############
	# now do the pratchet
	#
	# phangorn requires everything to be phyDat format
	mrp_phyDat <- phangorn::phyDat(matrix_final, type="USER",
		levels = levels,  ambiguity = ambiguity)
	#
	# and now we can do parsimony
	outTree <- phangorn::pratchet(mrp_phyDat, trace = trace)
	#	
	if(reduce_collapse){
		# handle properly depending on if multi phylo or not
		if (inherits(outTree, "multiPhylo")) {
			outTree <- lapply(outTree,
				expand_collapsed_clades_post_pratchet,
				saved_sets = saved_sets, 
				expected_num_OTUs = nrow(char_matrix)
				)		
			class(outTree) <- "multiPhylo"
		}else{
			outTree <- expand_collapsed_clades_post_pratchet(
				tree = outTree, 
				saved_sets = saved_sets, 
				expected_num_OTUs = nrow(char_matrix)
				)
			}
	}else{
		res <- outTree
		}
	return(res)
	}


expand_collapsed_clades_post_pratchet<-function(
		tree, saved_sets, expected_num_OTUs
		){
	############################################
	# now replace each set in turn
	for( i in 1:length(saved_sets)){
		# get the set name
		set_name <- names(saved_sets)[i]
		# get the set name with addendums for the tips
		tip_names <- paste0(set_name, "_", c("a", "b"))
		# find the two tips with the corresponding name
		whichReplace <- sapply(tree$tip.label,function(x)
			any(x == tip_names))
		whichReplace <- which(whichReplace)
		# 
		# if these tips are *somehow* not direct children of same mother node
			# collapse all edges that are in their way
		#
		# do they have the same mom node?
		mom_nodes <- tree$edge[match(whichReplace, tree$edge[,2]),1]
		#
		if(mom_nodes[1] != mom_nodes[2]){
			tree <- collapse_all_nodes_between(
				tree = tree, tip_labels = tip_names)
			}
		###########
		# now let's add in all labels we removed as new descendants of the mom node
		tips_to_add <- saved_sets[i]
		#
		for(j in 1:length(tips_to_add)){
			# every cycle, need to reidentify which tips are to be replaced
			whichReplace <- sapply(tree$tip.label,function(x)
				any(x == tip_names))
			whichReplace <- which(whichReplace)		
			# get the mom node
			mom_node <- (tree$edge[match(whichReplace, tree$edge[,2]),1])[1]
			#
			tree <- add_single_taxon_to_tree(tree= tree, 
				new_tip_label = tips_to_add[j],
				# default location to add tip is the root
				nodeID = Ntip(tree) + 1
				)
			}
		#
		# drop the tips we don't want anymore
		for(j in tip_names){
			tree <- drop.tip(phy = tree, tip = j)
			}
		}
	# check that it has the correct number of taxa
	if(Ntip(tree) != expected_num_OTUs){
		stop(
			"Iterative collapsing & expansion of OTUs failed to produce a tree of the right number of OTUs"
			)
		}
	#################
	return(tree)
	}


collapse_all_nodes_between <- function(tree, tip_labels){
	# identify tips based on labels
	which_tips <- sapply(tree$tip.label,function(x)
		any(x == tip_labels))
	which_tips <- which(which_tips)
	#
	# first, find mom nodes
	mom_nodes <- tree$edge[match(which_tips, tree$edge[,2]),1]
	#
	if(mom_nodes[1] == mom_nodes[2]){
		stop("these tips already share the same direct mother node??")
		}
	########
	# find all nodes that aren't shared
	unshared_nodes <- find_unshared_nodes(
		tree = tree, tip_labels = tip_labels)
	#
	while(length(unshared_nodes) > 1 ){
		# pick one at random
		collapse_this_node <- unshared_nodes[1]
		#
		tree <- paleotree::collapseNodes(
			nodeID = collapse_this_node,
			tree = tree,
			collapseType = "backward") 
		#
		# remake unshared_nodes
		unshared_nodes_new <- find_unshared_nodes(
			tree = tree, tip_labels = tip_labels)
		# make sure the length of unshared_nodes changed
		if(length(unshared_nodes_new) >= length(unshared_nodes)){
			stop("somehow removing unshared nodes made the number of unshared nodes not decrease")
			}		
		# 
		unshared_nodes <- unshared_nodes_new
		}
	###########
	# check tree
	# identify tips based on labels
	which_tips <- sapply(tree$tip.label,function(x)
		any(x == tip_labels))
	which_tips <- which(which_tips)
	# first, find mom nodes
	mom_nodes <- tree$edge[match(which_tips, tree$edge[,2]),1]
	#######
	if(mom_nodes[1] != mom_nodes[2]){
		stop("tips of interest are still not have common ancestor even after collapsing!")
		}
	########
	return(tree)	
	}


find_unshared_nodes <- function(tree, tip_labels){
	# identify tips based on labels
	which_tips <- sapply(tree$tip.label,function(x)
		any(x == tip_labels))
	which_tips <- which(which_tips)
	#
	# first, find mom nodes
	mom_nodes <- tree$edge[match(which_tips, tree$edge[,2]),1]
	#
	# find all nodes leading up to each mom node
	# make them into a single vector
	node_lineages <- c(
		get_node_lineage(node = mom_nodes[1], tree = tree),
		get_node_lineage(node = mom_nodes[2], tree = tree)
		)		
	# find all nodes that aren't shared
	unshared_nodes <- names(table(node_lineages))[table(node_lineages) == 1]
	return(unshared_nodes)
	}


get_node_lineage <- function(tree, node){
	# find all nodes leading up to each mom node
	lineage <- node
	while(lineage[1] != (Ntip(tree) + 1)){
		mom_node <- tree$edge[tree$edge[,2] == lineage[1],1]
		lineage <- c(mom_node, lineage)
		}
	return(lineage)
	}	
