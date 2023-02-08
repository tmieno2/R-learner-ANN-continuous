# /*===========================================================
#' # CF
# /*===========================================================

run_analysis_CF <- function(train_data, test_data, id_var, x_vars, t_formula_fs, report = "eonr") {
  # /*+++++++++++++++++++++++++++++++++++
  #' ## Preparation
  # /*+++++++++++++++++++++++++++++++++++
  #--- create spline basis on the treatment variable ---#
  T_info <- prepare_T_mat(t_formula_fs, data = train_data)
  #--- dependent vairable name ---#
  y_var <- all.vars(t_formula_fs)[1]

  Y <-
    train_data[, ..y_var] %>%
    as.matrix() %>%
    .[, 1]
  X <- train_data[, ..x_vars] %>% as.matrix()
  W <- X

  #--- train CF ---#
  cf_trained <- run_CF_c_py(Y, T_info$T_sp, X, W, n_estimators = 1000)
  
  #--- find treatment effect ---#
  te_info <- get_te_CF(cf_trained, test_data, x_vars, id_var)

  T_seq <-
    copy(train_data) %>%
    .[, get_min_max_seq(eval(parse(text = T_info$T_var_name)))]

  response_data <- find_response_semi(T_seq, T_info, te_info)

  eonr_hat <-
    response_data %>%
    .[, profit := est * p_crop - p_input * T] %>%
    .[, .SD[which.max(profit), ], by = id_var] %>%
    .[, c(id_var, "T"), with = FALSE] %>%
    setnames("T", "opt_N_hat") %>%
    .[, method := "CF"]

  if (report == "eonr") {
    return(eonr_hat)
  } else if (report == "all") {
    return(list(eonr_hat = eonr_hat, response_data = response_data, trained_model = cf_trained, te_info = te_info))
  }
}

get_te_CF <- function(trained_model, test_data, x_vars, id_var, alpha) {
  X_test <- test_data[, ..x_vars] %>% as.matrix()

  te_hat <-
    trained_model$const_marginal_effect(X_test)

  te_hat_se <-
    trained_model$const_marginal_effect_interval(X_test, alpha = 0.05) %>%
    data.table(se = .) %>%
    rowwise() %>%
    mutate(se = list(
      data.table(se)
    )) %>%
    ungroup() %>%
    mutate(type = c("lower bound", "higher bound"))

  return(list(te_hat = data.table(te_hat), te_hat_se = te_hat_se, id_data = test_data[, ..id_var]))
}

find_response_semi <- function(T_seq, T_info, te_info) {
  eval_T <-
    data.table(T = T_seq) %>%
    setnames("T", T_info$T_var_name) %>%
    predict(T_info$gam_setup, newdata = ., type = "lpmatrix") %>%
    #* get rid of the intercept
    .[, -1] %>%
    data.table()

  curv_data <-
    as.matrix(te_info$te_hat) %*% t(as.matrix(eval_T)) %>%
    data.table() %>%
    setnames(names(.), c(as.character(T_seq))) %>%
    .[, id := 1:.N] %>%
    melt(id.var = "id") %>%
    setnames(c("variable", "value"), c("T", "est")) %>%
    .[, T := as.numeric(as.character(T))]

  id_data_to_merge <-
    te_info$id_data %>%
    .[, id := 1:.N]

  final_data <-
    id_data_to_merge[curv_data, on = "id"] %>%
    .[, id := NULL]

  return(final_data)
}

# /*===========================================================
#' # RF
# /*===========================================================
run_analysis_RF_ranger <- function(y_var, t_var, x_vars, id_var, train_data, test_data, p_crop, p_input, report = "eonr") {
  num_cov <- length(x_vars) + 1
  ranger_formula <-
    paste0(
      y_var, " ~ ", t_var, "+",
      paste0(x_vars, collapse = "+")
    )

  # defualt mtry is sqrt(num_cov)
  x_weights <-
    (sqrt(num_cov) - 1) / length(x_vars) %>%
      rep(., length(x_vars))

  trained_RF <-
    ranger(
      ranger_formula,
      data = train_data,
      num.trees = 1000,
      split.select.weights = c(1, x_weights)
    )

  T_seq <-
    copy(train_data) %>%
    .[, get_min_max_seq(eval(parse(text = t_var)))]

  N_data <- data.table(N = T_seq)

  X_test <-
    test_data[, c(x_vars, id_var), with = FALSE] %>%
    expand_grid_df(N_data, .) %>%
    data.table()

  response_data <-
    copy(X_test) %>%
    .[, y_hat := predict(trained_RF, data = .)$predictions]

  eonr_hat <-
    response_data %>%
    .[, pi_hat := (p_crop * y_hat) - (p_input * N)] %>%
    .[, .SD[which.max(pi_hat), ], by = id_var] %>%
    setnames("N", "opt_N_hat") %>%
    .[, c(id_var, "opt_N_hat"), with = FALSE] %>%
    .[, method := "RF"]

  # ggplot(st_as_sf(opt_EONR)) +
  #   geom_sf(aes(fill = N))

  # ggplot(opt_EONR[temp_id == 5, ]) +
  #   geom_line(aes(y = y_hat, x = N))

  # dup_units <- test_data[, .N, by = aunit_id] %>%
  #   .[N == 2, aunit_id]

  if (report == "eonr") {
    return(eonr_hat)
  } else if (report == "all") {
    return(list(eonr_hat = eonr_hat, response_data = response_data, trained_model = trained_RF))
  }
}

# /*===========================================================
#' # BRF
# /*===========================================================

run_analysis_BRF <- function(y_var, t_var, x_vars, id_var, train_data, test_data, p_crop, p_input, report = "eonr") {
  Y <- train_data[, ..y_var] %>% as.matrix()
  cov_ls <- c(t_var, x_vars)
  X <- train_data[, ..cov_ls] %>% as.matrix()
  X_test <- test_data[, ..cov_ls] %>% as.matrix()

  trained_brf <- boosted_regression_forest(X, Y, num.trees = 1000, tune.parameters = "all")

  T_seq <-
    copy(train_data) %>%
    .[, get_min_max_seq(eval(parse(text = t_var)))]

  N_data <- data.table(N = T_seq)

  X_test <-
    test_data[, c(x_vars, id_var), with = FALSE] %>%
    expand_grid_df(N_data, .) %>%
    data.table()

  response_data <-
    copy(X_test) %>%
    .[, y_hat := predict(trained_brf, newdata = .[, cov_ls, with = FALSE])]

  eonr_hat <-
    response_data %>%
    .[, pi_hat := (p_crop * y_hat) - (p_input * N)] %>%
    .[, .SD[which.max(pi_hat), ], by = id_var] %>%
    setnames("N", "opt_N_hat") %>%
    .[, c(id_var, "opt_N_hat"), with = FALSE] %>%
    .[, method := "BRF"]

  # ggplot(st_as_sf(opt_EONR)) +
  #   geom_sf(aes(fill = opt_N_hat)) +
  #   scale_fill_viridis_c()

  # dup_units <- test_data[, .N, by = aunit_id] %>%
  #   .[N == 2, aunit_id]
  if (report == "eonr") {
    return(eonr_hat)
  } else if (report == "all") {
    return(list(eonr_hat = eonr_hat, response_data = response_data, trained_model = trained_brf))
  }
}
