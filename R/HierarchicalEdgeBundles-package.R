#' Visualize graphs using Hierarchical Edge Bundles
#'
#' Visualize graphs using a crude version of the hierarchical edge bundling
#' method.
#' Hierarchical edge bundling visualizes graphs by guiding edges along along 
#' a hierarchical tree of the nodes. A bundling parameter controls how tightly
#' the edges follow the tree.
#' For details confer the reference below.
#'
#' @name HierarchicalEdgeBundles-package
#' @aliases HierarchicalEdgeBundles-package HierarchicalEdgeBundles HEB
#' @author
#'   Anders Ellern Bilgrau \cr
#'   Maintainer: Anders Ellern Bilgrau <anders.ellern.bilgrau@@gmail.com>
#' @docType package
#' @references
#'   Holten, Danny. "Hierarchical edge bundles: Visualization of adjacency
#'   relations in hierarchical data." Visualization and Computer Graphics, IEEE
#'   Transactions on 12.5 (2006): 741-748.
#' @seealso
#'   Core user function: \code{\link{plotHEB}}
#' @examples
#' library("igraph")
#' library("ape")
#'
#' graph <- watts.strogatz.game(1, size = 10, nei = 2, p = 0.5)
#' adj.mat <- get.adjacency(graph)
#' phylo <- as.phylo(hclust(as.dist(1-adj.mat)))
#'
#' plotHEB(graph, phylo, type = "fan")
#' @import igraph
#' @importFrom ape plot.phylo
NULL
