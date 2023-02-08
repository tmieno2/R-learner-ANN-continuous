run_CF_c <- function(Y,
                     T,
                     X,
                     W,
                     X_test,
                     model_y,
                     model_t,
                     t_gam_formula,
                     cv = 5,
                     criterion = 5,
                     n_estimators = 2000,
                     min_samples_leaf = 10,
                     min_impurity_decrease = 0.001,
                     random_state = 3747823) {
  gam_setup <- gam(t_gam_formula, data = data)

  T_sp <-
    predict(gam_setup, data = data, type = "lpmatrix") %>%
    #* get rid of the intercept
    .[, -1] %>%
    data.table() %>%
    setnames(names(.), paste0("T_", 1:ncol(.)))

  te_hat_cf <-
    run_CF_c(
      Y = Y,
      T = T_sp,
      X = X,
      W = W,
      X_test = X,
      model_y = model_y,
      model_t = model_t,
      cv = cv,
      criterion = criterion,
      n_estimators = n_estimators,
      min_samples_leaf = min_samples_leaf,
      min_impurity_decrease = min_impurity_decrease,
      random_state = random_state
    ) %>%
    data.table()
}
