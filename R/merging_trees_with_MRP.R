# merging_trees_with_MRP 05-29-19
# GPL v3
# Authors: David Bapst and Luna Luisa Sanchez Reyes

merging_trees_with_MRP <- function(
		tree_backbone, tree_secondary, 
		backbone_reweighting = 1,
		trace = 0){
	###############################
	# remove branch lengths
	tree_backbone$edge.length <- NULL
	tree_secondary$edge.length <- NULL
	#
	##########################
	# add an artificial outgroup to both trees
	# make an artificial 1 tip tree
	outgroup <- list(
		edge = matrix(c(2,1),1,2),
		tip.label = "placeholder_artificial_outgroup",
		edge.length = NULL,
		Nnode = 1)
	class(outgroup)<-"phylo"
	#
	tree_backbone <- bind.tree(tree_backbone, outgroup)
	tree_secondary <- bind.tree(tree_secondary, outgroup) 
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
	################################
	# reweight the representation in the backbone (so its higher than the secondary?
	# nice!
	if(backbone_reweighting > 1){
		# repeat the columns in the backbone mrp X times, so that the same node is being
	# considered as multiple characters (thus weighting those groups 'higher'
	   mrp_backbone <- mrp_backbone[,]
		}
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
	#####################################
	# and then we do MRP
	# or really just a basic parsimony search with phangorn
	#
	# phangorn requires everything to be phyDat format
	mrp_phyDat <- phangorn::phyDat(mrp_full, type="USER", 
		levels = 0:1,  ambiguity = "?")
	#
	# and now we can do parsimony
	supertrees_out <- phangorn::pratchet(mrp_phyDat, trace = trace)
	#
	# is it one tree or more?
		# root the trees based on artificial outgroup
			#and remove artificial outgroup
		# code modeled on phangorn's supertree functions
	if (inherits(supertrees_out, "multiPhylo")) {
		supertrees_out  <- lapply(supertrees_out , root, 
			"placeholder_artificial_outgroup")
		supertrees_out <- lapply(supertrees_out , drop.tip, 
			"placeholder_artificial_outgroup")	
        class(supertrees_out) <- "multiPhylo"
      }else{
        supertrees_out <- root(supertrees_out, "placeholder_artificial_outgroup")
        supertrees_out <- drop.tip(supertrees_out, "placeholder_artificial_outgroup")
      }
	#
	# and voilla, you'd get a tree sample you can do
		# a strict consensus on, or whatever
	#
	return(supertrees_out)
	}

