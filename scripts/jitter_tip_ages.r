#####Jitter Nodes######################################
#####This script is a modified version of the##########
#####script by Thompson et al., see:#########################
#####https://github.com/ammonthompson/phylogeo_epi_cnn/tree/main/real_data_analysis/scripts/jitter_tip_ages.r


input_tree = "nadeau_combined_chains_mcc.typed.node.tree"
outfile <- "nadeau.jittered.nexus"

tree <- ape::read.nexus(input_tree)

#root-to-node distance for every node, indexed by node id
node_depth  <- ape::node.depth.edgelength(tree)
tree_height <- max(node_depth)

key       <-  node_depth
dup_nodes <- which(duplicated(key) | duplicated(key, fromLast = TRUE))

if (length(dup_nodes) == 0) {
  message("No duplicate node depths found")
} else {

  subtend_edges <- which(tree$edge[, 2] %in% dup_nodes)
  #additive jitter instead of multiplied
  n <- length(subtend_edges)
  tree$edge.length[subtend_edges] <-
    tree$edge.length[subtend_edges] +
    runif(n, 1e-5 + tree_height, 1e-4 + tree_height)
  
  
}

ape::write.nexus(tree, file = outfile, translate = TRUE)
