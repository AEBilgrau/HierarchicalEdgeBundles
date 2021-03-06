% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/plotHEB.R
\name{plotHEB}
\alias{plotHEB}
\title{Hierarchical Edge Bundling}
\usage{
plotHEB(graph, phylo, beta = 0.8, include.mrca = FALSE, simplify = FALSE,
  ..., args.lines = list(), args.points = list(pch = 16, cex = 0.1),
  debug = FALSE, e.cols, v.use.only, e.use.only)
}
\arguments{
\item{graph}{A \code{igraph} object to be drawn.}

\item{phylo}{A \code{phylo} object defining the hierarchy.}

\item{beta}{A number between 0 and 1 controlling the bundling strength.}

\item{include.mrca}{Should the only the most recent common ancestor be used
in the shortest path used for the splines?}

\item{simplify}{Simplify the paths by taking the convex hull of the control
points. Can sometimes yield better results.}

\item{...}{Arguments passed to \code{\link[ape]{plot.phylo}}.}

\item{args.lines}{A list of arguments passed to \code{lines}.}

\item{args.points}{A list of arguments passed to \code{points}.}

\item{debug}{Plot some extra debug info.}

\item{e.cols}{A character vector giving the colors of the edges.
Overrides any use of \code{args.lines$col}.}

\item{v.use.only}{An integer vector giving the nodes from which edges are
to be drawn. E.g. \code{use.only = 1} will only plot the edges from
vertex 1.}

\item{e.use.only}{An integer vector giving the edges to be drawn.
E.g. \code{use.only = 1} will only draw the first edge.}
}
\value{
Plots to a new device.
}
\description{
Visualization of networks using hierarchical edge bundles referenced below.
Plots a graph using a hierarchical tree as guides for the edges.
}
\examples{
library("igraph")
library("ape")
n <- 10

# Create a graph
corr <- cor(matrix(rnorm(n^2), n, n))
rownames(corr) <- colnames(corr) <- LETTERS[1:n]
adj.mat <- abs(corr)
graph <- graph.adjacency(adj.mat, mode = "un", weighted = TRUE, diag = FALSE)
E(graph)$color <- ifelse(corr[lower.tri(corr)] < 0, "blue", "red")
E(graph)$width <- 10*adj.mat[lower.tri(adj.mat)]

# Generate dendrogram
phylo <- as.phylo(hclust(as.dist(1 - adj.mat), method = "complete"))

# Plot 1 - the graph and the tree
par(mfrow = c(1,2))
plot(graph, layout = layout.circle)
plot(phylo, type = "fan")

# Plot 2 - Combining the two
plotHEB(graph, phylo, type = "fan", e.cols = E(graph)$color)


par(mfrow = c(1, 2), mar = c(0, 0, 2, 0))
plotHEB(graph, phylo, type = "fan", beta = 0.95,
                            args.lines = list(col = "#4682B470"),
                            main = "High bundling strength")
plotHEB(graph, phylo, type = "fan", beta = 0.55,
                            args.lines = list(col = "#4682B470"),
                            main = "Low bundling strength")

# Plot 3 - Extra control of plotting and debugging
par(mfrow = c(1,2))
plot(phylo, type = "unrooted")
plotHEB(graph, phylo, type = "unrooted", beta = 0.8,
                            v.use.only = 1, debug = FALSE,
                            args.lines = list(col = "red", lwd = 2))
}
\author{
Anders Ellern Bilgrau <anders.ellern.bilgrau (at) gmail.com>
}
\references{
Holten, Danny. "Hierarchical edge bundles: Visualization of adjacency
  relations in hierarchical data." Visualization and Computer Graphics, IEEE
  Transactions on 12.5 (2006): 741-748.
}
\seealso{
See \code{\link[ape]{plot.phylo}}.
}

