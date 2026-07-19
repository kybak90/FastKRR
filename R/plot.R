utils::globalVariables(c("yhat", "x", "y"))

#' Plot method for fitted Kernel Ridge Regression (KRR) models
#'
#' @description
#' Visualizes the fitted regression curve from a Kernel Ridge Regression (KRR) model.
#' Automatically generates predictions on a regular grid (1.2 times the training
#' sample size) and overlays them with the training data.
#'
#' @details
#' Currently, \code{plot.krr} supports only uni-variate inputs (\eqn{d = 1}).
#' For multivariate settings (\eqn{d \ge 2}), the plot method will return an error,
#' and users are encouraged to manually slice their data to visualize conditional
#' main effects.
#'
#' @param x A fitted KRR model (class \code{"krr"}) returned by \code{\link{fastkrr}}.
#' @param show_points Logical; if \code{TRUE}, displays the original training data points
#'   as a background layer. Default is \code{TRUE}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A \code{ggplot} object showing the fitted regression line and training data.
#'
#' @seealso \code{\link{fastkrr}}, \code{\link{predict.krr}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n = 1000
#' rho = 1
#' d = 1
#' X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d); colnames(X) = paste0("X", seq_len(d))
#' y = sin(2 * pi * rowMeans(X)^3) + rnorm(n, mean = 0, sd = 0.1)
#'
#' data = data.frame(X, y = y)
#'
#' model_exact = fastkrr(data = data, response = "y",
#'                        kernel = "gaussian", rho = rho, opt = "exact", verbose = FALSE)
#' plot(model_exact)
#' }
#'
#' @importFrom ggplot2 ggplot geom_point geom_line labs theme_minimal theme element_text aes
#' @export
plot.krr = function(x, show_points = TRUE, ...) {

  if (!inherits(x, "krr"))
    stop("Input must be a fitted model of class 'krr'.")

  model = x
  X = model$x
  y = as.vector(model$y)
  d = if (is.matrix(X)) ncol(X) else 1
  if (d != 1)
    stop("plot.krr() currently supports only 1D input.")

  n_grid = max(50, round(1.2 * nrow(as.matrix(X))))
  newdata = seq(min(X), max(X), length.out = n_grid)

  pred = predict(model, as.matrix(newdata))

  ord = order(as.numeric(newdata))
  df_pred = data.frame(x = as.numeric(newdata[ord]),
                       yhat = as.numeric(pred[ord]))
  df_train = data.frame(x = as.numeric(X), y = as.numeric(y))

  p = ggplot() +
    {if (show_points)
      geom_point(data = df_train, aes(x = x, y = y),
                 alpha = 0.25, color = "gray")} +
    geom_line(data = df_pred, aes(x = x, y = yhat),
              color = "#0072B2", linewidth = 1) +
    labs(
      title = paste0("FastKRR Fit (opt = '", model$opt, "')"),
      x = "x", y = "Predicted f(x)"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5))

  return(p)
}


