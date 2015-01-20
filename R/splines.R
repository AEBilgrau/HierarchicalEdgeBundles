#' Plot piecewise B-splines
#' 
#' Add piecewise B-splines of supplied order to a plot. The function
#' uses the \code{lines} function to plot to the current device.
#' 
#' @param x A vector giving the x-coordinates of the control points. 
#'   Alternatively, \code{x} can be a \code{list} or \code{data.frame} where
#'   the two first elements (or columns) are used as \code{x} and \code{y}.
#' @param y An optional vector giving the y-coordinates of the control points. 
#' @param order the order of the B-splines. Default is cubic (\code{order = 4}).
#' @param n.evals integer giving the number of points the spline should be 
#'   evaluated at. Default is 100.
#' @return 
#'   Returns a \code{list} of length 2 with entries \code{x} and \code{y} of 
#'   x and y coordinates corresponding to the evaluations of the B-spline.
#' @seealso See also \code{\link{lines}} and \code{\link[Hmisc]{bezier}}.
#' @author Anders Ellern Bilgrau
#' @examples
#' n <- 20
#' x <- list(x = cumsum(rnorm(2*n)), y = cumsum(rnorm(2*n)))
#' plot(x, pch = 16, col = "grey")
#' i <- seq_along(x[[1]] - 1)
#' with(x, arrows(x[i], y[i], x[i+1], y[i+1], col = "grey"))
#' lines(bSpline(x, order = 4), col = "red")
#' lines(bSpline(x, order = 5), col = "blue")
#' lines(bSpline(x, order = 6), col = "green")
#' @export
bSpline <- function(x, y, order = 4, n.evals = 100) {
  if (is.numeric(x)) stopifnot(length(x) >= order)
  if (missing(y)) {
    if (!is.list(x)) {
      stop("If y is not supplied, then x must be a list or data.frame.")
    }
    y <- x[[2]]
    x <- x[[1]]
  }
  mat <- cbind(x, y)
  n <- nrow(mat)
  t <- seq(0, 1, length.out = n - order + 2)
  y <- seq(0, 1, length.out = n.evals)
  ans <- deboor(mat, t, y, order)
  return(list(x = ans[,1], y = ans[,2]))
}


#' De Boor's algorihmn for evaluating B-splines
#' 
#' Evaluate piecewise B-splines a various parameter values. Internal function.
#' 
#' @param x A n by 2 matrix of control points.
#' @param t A vector of support points
#' @param y A vector of evaluation points. Should be sorted!
#' @param order The spline order. Cubic splines is order 4.
#' @return A 2 by \code{length(y)} matrix of the evaluated points on the spline.
#' @examples
#' order <- 3
#' x <- matrix(rnorm(12), 6, 2)
#' t <- seq(0, 1, l = nrow(x) - order + 2)
#' y <- seq(0, 1, l = 50)
#' print(res <- Bmisc:::deboor(x, t, y, order))
#' plot(x, type = "b", col = "grey", pch = 16)
#' points(res, col = "red", type = "l", pch = 16, lwd = 3)
#' @keywords internal
deboor <- function(x, t, y, order) {
  m <- nrow(x)
  n <- length(y)
  X <- Y <- matrix(0, order, order)
  t <- c(rep(min(t), order - 1), t, rep(max(t), order - 1))
  ans <- matrix(NA, n, ncol(x))
  ans[n, ] <- x[m, ]
  for (l in seq_len(n)) {
    k  <- max(which(y[l] >= t))
    if (k > m) {
      break
    }
    
    X[, 1] <- x[(k - order + 1):k, 1]
    Y[, 1] <- x[(k - order + 1):k, 2]
    for (i in 2:order) {
      for (j in i:order) {
        num <- y[l] - t[k - order + j]
        s <- t[k + j - i + 1] - t[k - order + j]
        alpha <- (num != 0)*(num/s)
        X[j, i] = (1 - alpha)*X[j - 1, i - 1] + alpha*X[j, i - 1]
        Y[j, i] = (1 - alpha)*Y[j - 1, i - 1] + alpha*Y[j, i - 1]
      }
    }
    
    ans[l, 1] <- X[order, order];
    ans[l, 2] <- Y[order, order];
  }
  return(ans)
}




#' Straightened Spline curves
#' 
#' Function for straightning Bezier or B-spline curves.
#' 
#' @param x A numeric vector giving the x-coordinates of the control points. 
#'   Alternatively, it can be a \code{list} (or \code{data.frame}) containing
#'   two vectors of control points.
#' @param y A numeric vector giving the y-coordinates of the control points. 
#' @param beta A numeric values that control amount of straightening.
#'   A \code{beta} of \eqn{1} is a normal spline curve and where as a value of
#'   \eqn{0} is a straight line.
#' @param order The order of the B-spline. See \code{\link{bSpline}}.
#' @param n.evals The number of evaluations of the Bezier curve.
#' @return
#'   Returns a \code{list} of length 2 with entries \code{x} and \code{y} of 
#'   x and y coordinates corresponding to the evaluations of the B-spline.
#' @author Anders Ellern Bilgrau
#' @note The straightening is done by modifying the control points as described
#'   in the reference.
#' @references
#'   Holten, Danny. "Hierarchical edge bundles: Visualization of adjacency 
#'   relations in hierarchical data." Visualization and Computer Graphics, 
#'   IEEE Transactions on 12.5 (2006): 741-748.
#' @examples
#' # Bezier curves
#' x <- list(x = cumsum(rnorm(7)), 
#'           y = cumsum(rnorm(7)))
#' plot(x, type = "b", col = "grey", pch = 16)
#' for (b in seq(0.1, 1, l = 9)) 
#'   lines(Bmisc:::straightenedBezier(x, beta = b))
#' @keywords internal
straightenedBezier <- function(x, y, beta = 0.5, n.evals = 100) {
  stopifnot(require("Hmisc"))
  if (missing(y)) {
    if (!is.list(x)) {
      stop("If y is not supplied, then x must be a list or data.frame.")
    }
    y <- x[[2]]
    x <- x[[1]]
  }
  n <- length(x)
  
  # Straigten by control points
  sequ <- (seq_len(n) - 1)/(n - 1)
  x <- beta*x + (1 - beta)*( x[1] + sequ*(x[n] - x[1]) )
  y <- beta*y + (1 - beta)*( y[1] + sequ*(y[n] - y[1]) )
  
  # Construct and return Bezier
  return(Hmisc::bezier(x, y, evaluation = n.evals))
}

#' @rdname straightenedBezier
#' @examples
#' # B-splines
#' x <- cumsum(rnorm(7))
#' y <- cumsum(rnorm(7))
#' plot(x, y, type = "b", col = "grey", pch = 16)
#' for (b in seq(0.1, 0.9, l = 9)) 
#'   lines(Bmisc:::straightenedBSpline(x, y, beta = b))
#' @keywords internal
straightenedBSpline <- function(x, y, order = 4, beta = 0.5, n.evals = 100) {
  if (missing(y)) {
    y <- x[[2]]
    x <- x[[1]]
  }
  n <- length(x)
  
  # Straigten by control points
  sequ <- (seq_len(n) - 1)/(n - 1)
  x <- beta*x + (1 - beta)*( x[1] + sequ*(x[n] - x[1]) )
  y <- beta*y + (1 - beta)*( y[1] + sequ*(y[n] - y[1]) )
  
  # Construct and return B-spline
  return(bSpline(x, y, order = order, n.evals = n.evals))
}




