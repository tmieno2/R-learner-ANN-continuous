# /*===========================================================
#' # First stage models
# /*===========================================================
predict_ranger_c <- function(training_data, eval_data, y_var, w_vars) {
  rf_formula_y <- formula(paste0(y_var, " ~ ", paste0(w_vars, collapse = "+")))

  # === train ===#
  rf_trained_y <-
    ranger(
      rf_formula_y,
      data = training_data
    )

  #--- fit ---#
  eval_data[, y_hat := predict(rf_trained_y, data = eval_data)$predictions]

  return(eval_data[, .(id, y_hat)])
}

# /*+++++++++++++++++++++++++++++++++++
#' ## Regression Forest from the grf package
# /*+++++++++++++++++++++++++++++++++++
predict_rf_c <- function(training_data, eval_data, y_var, xw_var) {
  Y <- training_data[, ..y_var] %>% as.matrix()
  X <- training_data[, ..xw_var] %>% as.matrix()
  X_eval <- eval_data[, ..xw_var] %>% as.matrix()

  trained_rf <-
    grf::regression_forest(
      X,
      Y,
      num.trees = 1000,
      tune.parameters = "all",
      # tune.num.trees = 100,
      # tune.num.reps = 500
      num.thread = 1
    )

  #--- fit ---#
  eval_data[, y_hat := predict(trained_rf, newdata = X_eval)$predictions]

  return(eval_data)
}

# /*+++++++++++++++++++++++++++++++++++
#' ## Boosted Regression Forest from the grf package
# /*+++++++++++++++++++++++++++++++++++
predict_brf_c <- function(training_data, eval_data, y_var, xw_var) {
  Y <- training_data[, ..y_var] %>% as.matrix()
  X <- training_data[, ..xw_var] %>% as.matrix()
  X_eval <- eval_data[, ..xw_var] %>% as.matrix()

  trained_rf <-
    grf::boosted_regression_forest(
      X,
      Y,
      num.trees = 1000,
      tune.parameters = "all",
      # tune.num.trees = 100,
      # tune.num.reps = 500,
      num.threads = 1
    )

  #--- fit ---#
  y_hat <- predict(trained_rf, newdata = X_eval)$predictions
  y_tilde <- eval_data[, ..y_var] - y_hat

  eval_data[, y_hat := y_hat]
  eval_data[, y_tilde := y_tilde]

  return(eval_data)
}

rf_t_ranger <- function(training_data, eval_data, t_var, w_vars) {
  rf_formula_t <- formula(paste0(t_var, " ~ ", paste0(w_vars, collapse = "+")))

  # === train ===#
  rf_trained_t <-
    ranger(
      rf_formula_t,
      data = training_data
    )

  #--- fit ---#
  eval_data[, t_hat := predict(rf_trained_t, data = eval_data)$predictions]
  setnames(eval_data, "t_hat", paste0(t_var, "_hat"))

  return(eval_data[, c("id", paste0(t_var, "_hat")), with = FALSE])
}

# /*===========================================================
#' # Second stage models
# /*===========================================================
final_stage_brf <- function(data_2nd_stage, x_vars) {
  X <- data_2nd_stage[, ..x_vars] %>% as.matrix()
  Y <- data_2nd_stage[, y_to_t] %>% as.matrix()
  weights <- data_2nd_stage[, weight]

  trained_BRF <- boosted_regression_forest(X, Y, sample.weights = weights, tune.parameter = "all")

  return_data <- data_2nd_stage[, theta_hat := predict(trained_BRF, newdata = X)][, .(id, theta_hat)]

  return(return_data)
}

final_stage_xgb <- function(data_2nd_stage, x_vars) {
  data_2nd_xgb <-
    xgb.DMatrix(
      data = data_2nd_stage[, ..x_vars] %>% as.matrix(),
      label = data_2nd_stage[, y_to_t],
      weight = data_2nd_stage[, weight]
    )

  # === train ===#
  xgb_trained_2nd <-
    xgboost(
      data = data_2nd_xgb,
      nrounds = 200,
      objective = "reg:squarederror"
    )

  data_2nd_stage[, theta_hat := predict(xgb_trained_2nd, newdata = data_2nd_xgb)]

  return(data_2nd_stage)
}


# /*===========================================================
#' # Internal
# /*===========================================================
#' Implement cross-fitting for the first stage of an R-leaner training
#'
#' @description An extension of expand.grid() to expand on two `data.frame`s
#' @param n: dependent variable
#' @param data_folds: dependent variable
#' @param y_var: dependent variable
#' @param t_var: treatment variable
#' @param w_var: all the covariates (including non-heterogeneity drivers)
#' @import data.table
#' @export

# y_var <- "Y"
# t_var <- "T"
# w_vars <- c("x1", "x2", "x3")
# model_y <- rf_y_ranger
# model_t <- rf_t_ranger

# cross_fit <- function(n, model_y, model_t, data_folds, y_var, t_var, w_vars) {
#   training_data <- analysis(data_folds[n, ]$splits[[1]])
#   training_data <- analysis(data_folds[n, ]$splits[[1]])
#   eval_data <- assessment(data_folds[n, ]$splits[[1]])

#   #--- E[Y|X] ---#
#   y_hat_data <- model_y(training_data, eval_data, y_var, w_vars)

#   #--- E[T|X] ---#
#   t_hat_data <- model_t(training_data, eval_data, t_var, w_vars)

#   #--- combine ---#
#   return_data <- t_hat_data[y_hat_data, on = "id"]

#   return(return_data)
# }

cross_fit <- function(n, model_y, model_t, data_folds, y_var, t_var, xw_var, id_var) {
  print(paste0("working on ", n, "th fold"))
  training_data <- analysis(data_folds[n, ]$splits[[1]])
  eval_data <- assessment(data_folds[n, ]$splits[[1]])

  #--- E[Y|X] ---#
  y_hat_data <- model_y(training_data, eval_data, y_var, xw_var)

  #--- E[T|X] ---#
  t_hat_data <-
    data.table(t_var = t_var) %>%
    rowwise() %>%
    mutate(t_hat = list(
      model_t(training_data, eval_data, y_var = t_var, xw_var) %>%
        .[, c(id_var, "y_hat", "y_tilde"), with = FALSE] %>%
        setnames("y_hat", paste0(t_var, "_hat")) %>%
        setnames("y_tilde", paste0(t_var, "_tilde"))
    )) %>%
    .$t_hat %>%
    reduce(left_join, by = id_var)

  all_vars <- c(id_var, y_var, "y_hat", "y_tilde", t_var, paste0(t_var, "_hat"), paste0(t_var, "_tilde"))

  return_data <-
    cbind(y_hat_data, t_hat_data) %>%
    .[, ..all_vars]

  return(return_data)
}

run_first_stage <- function(data_folds, model_y, model_t, y_var, t_var, xw_var, id_var) {
  data_2nd_stage <-
    #--- cross-fitting ---#
    lapply(
      seq_len(nrow(data_folds)),
      function(x) cross_fit(x, model_y, model_t, data_folds, y_var, t_var, xw_var, id_var)
    ) %>%
    rbindlist() %>%
    collapse::fgroup_by(id_var) %>%
    collapse::fmean()

  return(data_2nd_stage)
}


prepare_2nd_stage_data <- function(t_formula_fs, model_y, model_t, x_var, id_var, data, num_fold = 5, num_repeats = 1) {
  #--- create spline basis on the treatment variable ---#
  T_info <- prepare_T_mat(t_formula_fs, data = data)

  #--- dependent vairable name ---#
  y_var <- all.vars(t_formula_fs)[1]

  #--- treatment variable names ---#
  t_var <- names(T_info$T_sp)

  #--- add the data spline basis data ---#
  data_with_T <- cbind(copy(data), T_info$T_sp)

  #--- data for cross-fitting ---#
  data_folds <- rsample::vfold_cv(data_with_T, v = 5, repeats = num_repeats)

  #--- first stage models ---#
  data_ss <-
    run_first_stage(
      data_folds,
      model_y = model_y,
      model_t = model_t,
      y_var = y_var,
      t_var = t_var,
      xw_var = x_var,
      id_var = id_var
    )

  data_reg <-
    data_ss[data, on = id_var] %>%
    .[, c(id_var, "y_tilde", paste0(t_var, "_tilde"), x_var, t_var, y_var, "Y_det"), with = FALSE]

  return(data_reg)
}

prepare_eval_data <- function(t_formula_fs, data) {
  t_seq <- data[, seq(min(T), max(T), length = 100)]

  eval_data <-
    prepare_T_mat(
      t_formula_fs,
      data = data.table(
        T = t_seq,
        Y = 1
      )
    ) %>%
    .$T_sp %>%
    .[, T := t_seq]
}

# /*===========================================================
#' # Misc
# /*===========================================================
expand_grid_df <- function(data_1, data_2) {
  data_1_ex <-
    data_1[rep(1:nrow(data_1), each = nrow(data_2)), ] %>%
    data.table() %>%
    .[, rowid := 1:nrow(.)]
  data_2_ex <-
    data_2[rep(1:nrow(data_2), nrow(data_1)), ] %>%
    data.table() %>%
    .[, rowid := 1:nrow(.)]
  expanded_data <-
    data_1_ex[data_2_ex, on = "rowid"] %>%
    .[, rowid := NULL]
  if ("tbl" %in% class(data_1)) {
    expanded_data <- as_tibble(expanded_data)
  }
  if ("rowwise_df" %in% class(data_1)) {
    expanded_data <- rowwise(expanded_data)
  }
  return(expanded_data)
}

prepare_T_mat <- function(gam_formula, data) {
  gam_setup <- gam(gam_formula, data = data)

  T_var_name <- all.vars(gam_formula)[2]

  T_sp <-
    predict(gam_setup, data = data, type = "lpmatrix") %>%
    #* get rid of the intercept
    .[, -1] %>%
    data.table() %>%
    setnames(names(.), paste0("T_", 1:ncol(.)))

  return(list(gam_setup = gam_setup, T_sp = T_sp, T_var_name = T_var_name))
}

get_min_max_seq <- function(x, length) {
  seq(min(x), max(x), length = 200)
}



# /*===========================================================
#' #
# /*===========================================================
# === Quadratic-Plateau response
gen_yield_QP <- function(b0, b1, b2, Nk, N) {
  yield <- (N < Nk) * (b0 + b1 * N + b2 * N^2) + (N >= Nk) * (b0 + b1 * Nk + b2 * Nk^2)
  return(yield)
}

# === Quadratic response
gen_yield_QD <- function(b0, b1, b2, N) {
  yield <- b0 + b1 * N + b2 * N^2
  return(yield)
}

unnest_dt <- function(data, col, by) {
  eval(parse(text = paste0("return_data <- data[, ", col, "[[1]], by = .(", paste0(by, collapse = ","), ")]")))

  return(return_data)
}
