gen_test_data <- function(g_formula = formula(~ I(exp(x1) / (1 + exp(x1))) + I(x3 / 4)), # formula that defines how x affects y indepdent of treatment
                          m_formula = formula(~ x1 + I(exp(x3) / (1 + exp(x3)) / 4)), # formula that defines how x affects treatment (d)
                          te_formula = formula(~ I(0.5 * d)), # formula that defines theta(x) * t
                          n_obs = 500,
                          n_vars = 20,
                          mu_x = 0,
                          vcov_x = NULL,
                          sigma = 1 # sd of the error term in the y equation
) {
  if (is.null(vcov_x)) {
    vcov_x <- matrix(rep(0, n_vars^2), nrow = n_vars)
    for (i in seq_len(n_vars)) {
      vcov_x[i, ] <- 0.7^abs(i - seq_len(n_vars))
    }
  }

  # === draw from multivariate normal ===#
  data <-
    mvrnorm(n_obs, mu = rep(0, n_vars), Sigma = vcov_x) %>%
    data.table() %>%
    setnames(names(.), paste0("x", 1:n_vars))

  # === generate d ===#
  if (m_formula == "independent") {
    data[, d := rnorm(n_obs)]
  } else {
    data[, d := model.frame(m_formula, data = data) %>% rowSums() + rnorm(n_obs)]
  }

  # === generate g ===#
  data[, g := model.frame(g_formula, data = data) %>% rowSums()]

  # === generate treatment effect ===#
  data[, te := model.frame(te_formula, data = data) %>% rowSums()]

  # === generate y ===#
  data[, y := te + g + rnorm(n_obs, sd = sigma)]

  return(data[])
}
