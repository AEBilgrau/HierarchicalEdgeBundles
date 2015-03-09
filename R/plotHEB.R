#' Hierarchical Edge Bundling
#'
#' Visualization of networks using hierarchical edge bundles referenced below.
#' Plots a graph using a hierarchical tree as guides for the edges.
#'
#' @param graph A \code{igraph} object to be drawn.
#' @param phylo A \code{phylo} object defining the hierarchy.
#' @param beta A number between 0 and 1 controlling the bundling strength.
#' @param include.mrca Should the only the most recent common ancestor be used
#'   in the shortest path used for the splines?
#' @param simplify Simplify the paths by taking the convex hull of the control
#'   points. Can sometimes yield better results.
#' @param ... Arguments passed to \code{\link[ape]{plot.phylo}}.
#' @param args.lines  A list of arguments passed to \code{lines}.
#' @param args.points A list of arguments passed to \code{points}.
#' @param debug Plot some extra debug info.
#' @param e.cols A character vector giving the colors of the edges.
#'   Overrides any use of \code{args.lines$col}.
#' @param v.use.only An integer vector giving the nodes from which edges are
#'   to be drawn. E.g. \code{use.only = 1} will only plot the edges from
#'   vertex 1.
#' @param e.use.only An integer vector giving the edges to be drawn.
#'   E.g. \code{use.only = 1} will only draw the first edge.
#' @return Plots to a new device.
#' @seealso See \code{\link[ape]{plot.phylo}}.
#' @author
#'   Anders Ellern Bilgrau <abilgrau (at) math.aau.dk>
#' @references
#'   Holten, Danny. "Hierarchical edge bundles: Visualization of adjacency
#'   relations in hierarchical data." Visualization and Computer Graphics, IEEE
#'   Transactions on 12.5 (2006): 741-748.
#' @examples
#' library("igraph")
#' library("ape")
#' n <- 10
#'
#' # Create a graph
#' corr <- cor(matrix(rnorm(n^2), n, n))
#' rownames(corr) <- colnames(corr) <- LETTERS[1:n]
#' adj.mat <- abs(corr)
#' graph <- graph.adjacency(adj.mat, mode = "un", weighted = TRUE, diag = FALSE)
#' E(graph)$color <- ifelse(corr[lower.tri(corr)] < 0, "blue", "red")
#' E(graph)$width <- 10*adj.mat[lower.tri(adj.mat)]
#'
#' # Generate dendrogram
#' phylo <- as.phylo(hclust(as.dist(1 - adj.mat), method = "complete"))
#'
#' # Plot 1 - the graph and the tree
#' par(mfrow = c(1,2))
#' plot(graph, layout = layout.circle)
#' plot(phylo, type = "fan")
#'
#' # Plot 2 - Combining the two
#' plotHEB(graph, phylo, type = "fan", e.cols = E(graph)$color)
#'
#' 
#' par(mfrow = c(1, 2), mar = c(0, 0, 2, 0))
#' plotHEB(graph, phylo, type = "fan", beta = 0.95,
#'                             args.lines = list(col = "#4682B470"),
#'                             main = "High bundling strength")
#' plotHEB(graph, phylo, type = "fan", beta = 0.55,
#'                             args.lines = list(col = "#4682B470"),
#'                             main = "Low bundling strength")
#'
#' # Plot 3 - Extra control of plotting and debugging
#' par(mfrow = c(1,2))
#' plot(phylo, type = "unrooted")
#' plotHEB(graph, phylo, type = "unrooted", beta = 0.8,
#'                             v.use.only = 1, debug = FALSE,
#'                             args.lines = list(col = "red", lwd = 2))
#' @import adephylo
#' @export
plotHEB <- function(graph,
                    phylo,
                    beta = 0.80,  # Bundling strength
                    include.mrca = FALSE,
                    simplify = FALSE,
                    ...,
                    args.lines = list(),
                    args.points = list(pch = 16, cex = 0.1),
                    debug = FALSE,
                    e.cols,
                    v.use.only,
                    e.use.only) {

  plot(phylo, edge.color = ifelse(debug,"grey","#00000000"), ...)
  last.plot.phylo <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  pos <- with(last.plot.phylo, data.frame(i = seq_along(xx), x = xx, y = yy))
  # Add points
  if (debug) {
    points(pos$x, pos$y, col = "black", pch = 16, cex = 0.5)
    text(pos$x, pos$y, col = "black", pos$i, cex = 1.5)
  }

  # Get edgelist (and convert to numeric if names are present)
  es <- get.edgelist(graph)
  if (!is.null(V(graph)$name)) {
    es <- structure(match(es, V(graph)$name), dim = dim(es))
  }

  if (!missing(v.use.only)) es <- es[es[, 1] %in% v.use.only, , drop = FALSE]
  if (!missing(e.use.only)) es <- es[e.use.only, , drop = FALSE]

  # Shortest paths (with start and end) or path through Most Recent Common
  # Ancestor
  sp2 <- sp.tips2(phylo, es[, 1], es[, 2], include.mrca = include.mrca,
                  useTipNames = TRUE)
  stopifnot(nrow(es) == length(sp2))
  # Add start and end
  sp <- lapply(seq_along(sp2), function(i) unname(c(es[i,1],sp2[[i]],es[i,2])))
  names(sp) <- names(sp2)

  # Get plot order if weights are present
  if (!is.null(E(graph)$weight)) {
    o <- order(abs(E(graph)$weight))
    indicies <- seq_along(sp)[o]
  } else {
    indicies <- seq_along(sp)
  } 

  # Plot spline curve for each path, 
  for (i in indicies) {
    path <- sp[[i]]
    d <- pos[path, ]
    if (simplify) {
      ch <- chull(d$x, d$y)
      d <- d[match(intersect(d$i, d$i[ch]), d$i), ]
    }
    ord <- ifelse(length(d$x) >= 4, 4, length(d$x))
    if (debug) lines(d$x, d$y)
    tmp <- straightenedBSpline(d$x, d$y, order = ord, beta = beta,
                               n.evals = 250)
    if (!missing(e.cols)) args.lines$col <- e.cols[i]
    do.call(lines, c(tmp, args.lines))
  }
  do.call(points, c(pos[seq_len(vcount(graph)), -1], args.points))
}





# MODIFIED sp.tips from adephylo!
# Will be updated in a later version of adephlyo
# This overwrites the current function
sp.tips2 <- function(x, tip1, tip2, useTipNames = FALSE,
                     quiet = FALSE, include.mrca = TRUE) {
  x <- as(x, "phylo4")
  if (is.character(checkval <- checkPhylo4(x)))
    stop(checkval)
  t1 <- getNode(x, tip1)
  t2 <- getNode(x, tip2)
  if (any(is.na(c(t1, t2))))
    stop("wrong tip specified")
  if (any(c(t1, t2) > nTips(x)))
    stop("specified nodes are internal nodes")
  if (length(t1) != length(t2)) {
    maxLength <- max(length(t1), length(t2))
    t1 <- rep(t1, length.out = maxLength)
    t2 <- rep(t2, length.out = maxLength)
  }
  toRemove <- (t1 == t2)
  if (sum(toRemove) > 0) {
    t1 <- t1[!toRemove]
    t2 <- t2[!toRemove]
    if (length(t1) == 0)
      stop("tip1 and tip2 are the same vectors")
    if (!quiet)
      warning("tip1 and tip2 are sometimes the same; erasing these cases")
  }
  N <- nTips(x)
  root <- getNode(x, N + 1)
  E <- x@edge
  allTips <- unique(c(t1, t2))
  pathTwoTips <- function(path1, path2) {
    cpath <- c(path1, rev(path2)) # <- CHANGED HERE, added rev()
    temp <- factor(cpath, levels = unique(cpath))
    CA <- temp[table(temp) == 2][1]
    CA <- as.integer(as.character(CA))
    path1 <- path1[1:(which(path1 == CA))]
    temp <- which(path2 == CA)
    if (temp == 1)
      return(path1)
    path2 <- path2[1:(temp - 1)]
    return(c(path1, path2))
  }
  pathTwoTips.no.mrca <- function(path1, path2) {
    cpath <- c(path1, rev(path2)) # <- CHANGED HERE, added rev()
    temp <- intersect(path1, path2)
    res <- setdiff(cpath, temp)
    return(res)
  }
  allPathToRoot <- lapply(allTips, function(i)
    .tipToRoot(x, i, root, include.root = TRUE))
  names(allPathToRoot) <- allTips
  allPath1 <- allPathToRoot[as.character(t1)]
  allPath2 <- allPathToRoot[as.character(t2)]
  if (include.mrca) {
    res <- lapply(1:length(allPath1), function(i)
      pathTwoTips(allPath1[[i]], allPath2[[i]]))
  }
  else {
    res <- lapply(1:length(allPath1), function(i)
      pathTwoTips.no.mrca(allPath1[[i]], allPath2[[i]]))
    temp.names <- names(res)
    temp <- sapply(res, function(vec) length(vec) > 0)
    res[temp] <- lapply(res[temp], function(vec) getNode(x, vec))
    names(res) <- temp.names
  }
  if (useTipNames) {
    names(res) <- paste(names(t1), names(t2), sep = "-")
  }
  else {
    names(res) <- paste(t1, t2, sep = "-")
  }
  return(res)
}
environment(sp.tips2) <- asNamespace("adephylo")

