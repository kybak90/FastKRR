set.seed(1)
n = 1000 ; d = 1

# generating data
X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
y = as.vector(sin(2*pi*X^3) + rnorm(n, 0, 0.1))

X_new = matrix(seq(0, 1, length.out = 1200), ncol = 1)
y_new = as.vector(sin(2*pi*X_new^3) + rnorm(1200, 0, 0.1))



fit_exact   = fastkrr(X, y, kernel = "gaussian", opt = "exact",   fastcv = TRUE)
fit_nystrom = fastkrr(X, y, kernel = "gaussian", opt = "nystrom", m = 100, fastcv = TRUE)
fit_pivoted = fastkrr(X, y, kernel = "gaussian", opt = "pivoted", m = 100, fastcv = TRUE)
fit_rff     = fastkrr(X, y, kernel = "gaussian", opt = "rff",     m = 200, fastcv = TRUE)


pred_exact   = predict(fit_exact,   X_new)
pred_nystrom = predict(fit_nystrom, X_new)
pred_pivoted = predict(fit_pivoted, X_new)
pred_rff     = predict(fit_rff,     X_new)


# Plots
library(dplyr)
library(tidyr)
library(ggplot2)
df_preds = data.frame(
  X = X_new[,1],
  Exact   = pred_exact,
  Nystrom = pred_nystrom,
  Pivoted = pred_pivoted,
  RFF     = pred_rff) %>%
  pivot_longer(-X, names_to = "Method", values_to = "y_hat")

df_preds$Method = factor(
  df_preds$Method,
  levels = c("Exact", "Nystrom", "Pivoted", "RFF"),
  labels = c("(a) Exact", "(b) Nystrom", "(c) Pivoted Cholesky", "(d) RFF"))

df_train = data.frame(X = X[,1], y = y)

ggplot() +
  geom_point(data = df_train, aes(x = X, y = y),
             color = "grey70", alpha = 0.35, size = 0.6) +
  geom_line(data = df_preds, aes(x = X, y = y_hat, color = Method),
            linewidth = 0.6) +
  facet_wrap(~ Method, ncol = 2) +
  scale_color_manual(values = c(
    "(a) Exact"   = "#1f77b4",
    "(b) Nystrom" = "#d62728",
    "(c) Pivoted Cholesky" = "#2ca02c",
    "(d) RFF"     = "#9467bd"
  )) +
  labs(x = "X", y = "y") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")





# Approx_kernel()
set.seed(1)
d = 1
n = 1000
m = 50
X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
K = make_kernel(X, kernel = "gaussian", rho = 1)

# Example: RFF approximation
#rff_pars = rff_random(m = m, d = d, rho = 1, kernel = "gaussian")
K_rff = approx_kernel(X = X, opt = "rff", kernel = "gaussian",
                      m = m, d = d, rho = 1,
                      #W = rff_pars$W, b = rff_pars$b,
                      n_threads = 1)

# Exapmle: Nystrom approximation
K_nystrom = approx_kernel(K = K, opt = "nystrom",
                          m = m, d = d, rho = 1,
                          n_threads = 1)

# Example: Pivoted Cholesky approximation
K_pivoted = approx_kernel(K = K, opt = "pivoted",
                          m = m, d = d, rho = 1)




library(tidymodels)
library(parsnip)
library(stats)
library(modeldata)

# Data analysis
data(ames)
ames = ames %>% mutate(Sale_Price = log10(Sale_Price))

set.seed(502)
ames_split = initial_split(ames, prop = 0.80, strata = Sale_Price)
ames_train = training(ames_split) # dim (2342, 74)
ames_test  = testing(ames_split) # dim (588, 74)

# Model spec
krr_spec = krr_reg(kernel = "gaussian", opt = "exact",
                   m = 50, eps = 1e-6, n_threads = 4,
                   rho = 1, penalty = tune()) %>%
  set_engine("fastkrr") %>%
  set_mode("regression")

# Define rec
rec = recipe(Sale_Price ~ Longitude + Latitude, data = ames_train)

# workflow
wf = workflow() %>%
  add_recipe(rec) %>%
  add_model(krr_spec)

# Define hyper-parameter grid
param_grid = grid_regular(
  dials::penalty(range = c(-10, -3)),
  levels = 5
)

# CV setting
set.seed(123)
cv_folds = vfold_cv(ames_train, v = 5, strata = Sale_Price)

# Tuning
tune_results = tune_grid(
  wf,
  resamples = cv_folds,
  grid = param_grid,
  metrics = metric_set(rmse),
  control = control_grid(verbose = TRUE, save_pred = TRUE)
)

# Result check
collect_metrics(tune_results)

# Select best parameter
best_params = select_best(tune_results, metric = "rmse")

# Finalized model spec using best parameter
final_spec = finalize_model(krr_spec, best_params)
final_wf = workflow() %>%
  add_recipe(rec) %>%
  add_model(final_spec)

# Finalized fitting using best parameter
final_fit = final_wf %>% fit(data = ames_train)

# Prediction
predict(final_fit, new_data = ames_test)
print(best_params)

plot(final_fit)

