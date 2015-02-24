
Hierarchical Edge Bundles in R
------------------------------

[![Build Status](https://api.travis-ci.org/AEBilgrau/HierarchicalEdgeBundles.svg?branch=master)](https://travis-ci.org/AEBilgrau/HierarchicalEdgeBundles)

This packages implements a crude version of the excellent visualization method of hierarchical edge bundling (Holten, 2006) in the statistical software R (R Core Team, 2014).

## Installation
To install the latest version of **HierarchicalEdgeBundles** directly from the master branch at GitHub, run 

```R
#install.packages("devtools")  # Uncomment if devtools is not installed
devtools::install_github("AEBilgrau/HierarchicalEdgeBundles")
```

The package is still under development and may be considered unstable. Be sure that you have the [package development prerequisites](http://www.rstudio.com/ide/docs/packages/prerequisites) if you wish to install the package from the source.

## Usage

The following example should run without problems.

```R
library("HierarchicalEdgeBundles")
library("igraph")
library("ape")

# Create a graph and a tree
graph <- watts.strogatz.game(1, size = 10, nei = 2, p = 0.5)
adj.mat <- get.adjacency(graph)
phylo <- as.phylo(hclust(as.dist(1-adj.mat)))

# Plot
plotHEB(graph, phylo, beta = 0.8, type = "fan")
```

## References
 * Holten, Danny. "Hierarchical Edge Bundles: Visualization of Adjacency 
   Relations in Hierarchical Data." Visualization and Computer Graphics, 
   IEEE Transactions on 12.5 (2006): 741-748.
   
 * R Core Team (2014). R: A language and environment for statistical
   computing. R Foundation for Statistical Computing, Vienna, Austria.
   http://www.R-project.org/.
