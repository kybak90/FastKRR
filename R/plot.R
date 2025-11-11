#' Plot method for fitted Kernel Ridge Regression (KRR) models
#'
#' @description
#' Visualizes fitted results from a Kernel Ridge Regression (KRR) model.
#' Automatically generates predictions on a regular grid
#' (120\% of training sample size) and overlays them with training data.
#'
#' @details
#' For multivariate inputs (\eqn{d \ge 2}), visualization requires fixing
#' all but one variable. For example, in 2D, one can plot
#' \eqn{f(x_1, x_2 = \bar{x}_2)} to examine the effect of \eqn{x_1}
#' while holding \eqn{x_2} at its mean.
#'
#' @param x A fitted KRR model (class \code{"krr"}) returned by \code{\link{fastkrr}}.
#' @param show_points Logical; if \code{TRUE}, displays the training data points.
#'   Default = \code{TRUE}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A ggplot object showing the fitted regression curve.
#'
#' @seealso \code{\link{fastkrr}}, \code{\link{predict.krr}}
#'
#' @examples
#' set.seed(1)
#' n = 1000
#' X = runif(n, 0, 1)
#' y = sin(2*pi*X^3) + rnorm(n, 0, 0.1)
#'
#' model_exact = fastkrr(X, y, kernel = "gaussian", rho = rho, opt = "exact", verbose = FALSE)
#' plot(model_exact)
#'
#' @import ggplot2
#' @export

plot.krr = function(x, show_points = TRUE, ...) {
  library(ggplot2)

  if (!inherits(x, "krr"))
    stop("Input must be a fitted model of class 'krr'.")

  X = attr(x, "x")
  y = as.vector(attr(x, "y"))
  d = if (is.matrix(X)) ncol(X) else 1
  if (d != 1)
    stop("plot.krr() currently supports only 1D input.")

  n_grid = max(50, round(1.2 * nrow(as.matrix(X))))
  newdata = seq(min(X), max(X), length.out = n_grid)

  pred = predict(x, as.matrix(newdata))

  ord = order(as.numeric(newdata))
  df_pred = data.frame(x = as.numeric(newdata[ord]),
                       yhat = as.numeric(pred[ord]))
  df_train = data.frame(x = as.numeric(X), y = as.numeric(y))

  p = ggplot() +
    {if (show_points)
      geom_point(data = df_train, aes(x = x, y = y),
                 alpha = 0.25, color = "gray")} +
    geom_line(data = df_pred, aes(x = x, y = yhat),
              color = "#0072B2", size = 1) +
    labs(
      title = paste0("FastKRR Fit (opt = '", attr(x, "opt"), "')"),
      x = "x", y = "Predicted f(x)"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5))

  return(p)
}
