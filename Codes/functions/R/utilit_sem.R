#' @title Run a series of analysis: regression and yield prediction
#'
#' @description Run a series of analysis: regression and yield prediction
#' @param data: data used for model training
#' @param id_var: name of the id variable
#' @param x_var: name of heterogeneity drivers
#' @param dml (logical): whether to do DML or not
#' @param t_formula_fs: gam formula to create spline bases from the treatment variable
#' @param model_y: prediction procedure used to orthogonalize y
#' @param model_t: prediction procedure used to orthogonalize t
#' @param cv_k: number of folds in the first-stage cross-fitting
#' @param cv_repeat: number of repeats in the first-stage cross-fitting
#' @param w: weight matrix
#' @param eval_data: data for which yields are predicted
#' @param t_seq: sequence of values of the treatment variable at which yields are predicted
#' @examples
se_semi <- function(data,
                    id_var,
                    x_var,
                    dml,
                    # /*+++++++++++++++++++++++++++++++++++
                    #' ## first-stage
                    # /*+++++++++++++++++++++++++++++++++++
                    t_formula_fs,
                    model_y = predict_brf_c,
                    model_t = predict_brf_c,
                    cv_k = 5,
                    cv_repeat = 3,
                    # /*+++++++++++++++++++++++++++++++++++
                    #' ## second-stage
                    # /*+++++++++++++++++++++++++++++++++++
                    w,
                    # /*+++++++++++++++++++++++++++++++++++
                    #' ## post-estimation analysis
                    # /*+++++++++++++++++++++++++++++++++++
                    eval_data,
                    t_seq,
                    p_crop,
                    p_N) {

  # /*+++++++++++++++++++++++++++++++++++
  #' ## Preparation
  # /*+++++++++++++++++++++++++++++++++++
  #--- create spline basis on the treatment variable ---#
  T_info <- prepare_T_mat(t_formula_fs, data = data)
  #--- dependent vairable name ---#
  y_var <- all.vars(t_formula_fs)[1]
  #--- treatment variable names ---#
  t_var <- names(T_info$T_sp)
  #--- add the data spline basis data ---#
  data_with_T <- cbind(copy(data), T_info$T_sp)

  if (dml == TRUE) {
    # /*+++++++++++++++++++++++++++++++++++
    #' ## First stage
    # /*+++++++++++++++++++++++++++++++++++
    data_folds <- rsample::vfold_cv(data_with_T, v = cv_k, repeats = cv_repeat)

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

    data_reg <- data_ss[data, on = id_var]
  } else {
    data_reg <- data_with_T
  }

  # /*+++++++++++++++++++++++++++++++++++
  #' ## Second stage
  # /*+++++++++++++++++++++++++++++++++++
  if (dml == TRUE) {
    y_reg_var <- "y_tilde"
    t_reg_var <- paste0(t_var, "_tilde")
  } else {
    y_reg_var <- y_var
    t_reg_var <- t_var
  }

  se_trained <- run_sem(y_reg_var, t_reg_var, x_var, w, data = data_reg, method = "GMM")

  # /*+++++++++++++++++++++++++++++++++++
  #' ## yield prediction
  # /*+++++++++++++++++++++++++++++++++++
  yield_eval <- predict_yield_profit_sem(se_trained, t_seq, T_info, eval_data, x_var, id_var)

  # /*+++++++++++++++++++++++++++++++++++
  #' ## find EONR
  # /*+++++++++++++++++++++++++++++++++++
  eonr <-
    yield_eval %>%
    .[, pi_hat := p_crop * y_hat - p_N * N] %>%
    .[, .SD[which.max(pi_hat), ], by = id_var] %>%
    .[, .(aunit_id, N)] %>%
    setnames("N", "opt_N_hat")

  # /*+++++++++++++++++++++++++++++++++++
  #' ## Assimilate
  # /*+++++++++++++++++++++++++++++++++++
  return_ls <-
    list(
      t_formula = t_formula_fs,
      se_trained = se_trained,
      yield_eval = yield_eval,
      eonr = eonr,
      T_info = T_info,
      dml = dml
    )

  return(return_ls)
}

# se_dml_semi(
#   data = data,
#   id_var = "aunit_id",
#   x_var = c("b1", "plateau", "Nk"),
#   t_formula_fs = formula(yield ~ s(N, k = 4), m = 2),
#   model_y = predict_brf_c,
#   model_t = predict_brf_c,
#   cv_k = 5,
#   cv_repeat = 3,
#   w = w,
#   eval_data = data,
#   t_seq = data[, seq(min(N), max(N), length = 100)]
# )


#' @title Predict yield at various levels of the treatment variable
#'
#' @description Predict yield at various levels of the treatment variable
#' @param se_trained: coefficient estimates from run_sem()
#' @param T_seq: sequence of treatment levels
#' @param T_info: treatment spline information
#' @param eval_data: data for which yields are estimated
#' @param x_var: name of heterogeneity drivers
#' @param id_var: name of the id variable
#' @examples
predict_yield_profit_sem <- function(se_trained, t_seq, T_info, eval_data, x_var, id_var) {
  eval_data <- eval_data[, c("aunit_id", x_var), with = FALSE]

  T_data <-
    data.table(T = t_seq) %>%
    setnames("T", T_info$T_var_name)

  T_mat_data <-
    T_data %>%
    predict(T_info$gam_setup, newdata = ., type = "lpmatrix") %>%
    #* get rid of the intercept
    .[, -1] %>%
    data.table() %>%
    setnames(names(.), se_trained$t_reg_var) %>%
    cbind(., T_data)

  yhat_data <-
    expand_grid_df(eval_data, T_mat_data) %>%
    .[, temp_y := 1] %>%
    setnames("temp_y", se_trained$y_reg_var)

  X_mat <-
    yhat_data %>%
    model.matrix(se_trained$reg_formula, data = .) %>%
    .[, -1]

  y_hat_ls <- X_mat %*% se_trained$beta[!(var %in% c("Intercept", "lambda")), coef]

  yhat_data[, y_hat := y_hat_ls]

  return(yhat_data)
}

#' @title Create a formula where treatment spline variables and heterogeneity drivers (x) are fully interacted
#'
#' @description Create a formula where treatment spline variables and heterogeneity drivers (x) are fully interacted
#' @param t_var: name of treatment spline variables
#' @param x_var: name of heterogeneity drivers
#' @examples
make_formula_se <- function(y_var, t_var, x_var) {
  formula_se <-
    str_c(
      "I(",
      str_c(x_var, rep(t_var, each = length(x_var)), sep = "*"),
      ")",
      collapse = "+"
    ) %>%
    str_c(y_var, " ~ ", str_c(t_var, collapse = "+"), "+", .) %>%
    formula()
  return(formula_se)
}

#' @title Run spatial error model (GMM or ML) using the spreg package in Python
#'
#' @description Run spatial error model (GMM or ML) using the spreg package in Python
#' @param se_formula: regression formula
#' @param w: weight matrix
#' @param data: data.frame
#' @param method: GMM (default) or ML
#' @export
#' @examples
run_sem <- function(y_reg_var, t_reg_var, x_var, w, data, method = "GMM") { # nolint
  reg_formula <- make_formula_se(y_reg_var, t_reg_var, x_var)
  data_mat <- model.frame(formula = reg_formula, data = data)
  x_var_names <- labels(terms(reg_formula))

  y <- data_mat[, 1] %>% as.matrix()
  x <- data_mat[, -1] %>% as.matrix()

  if (method == "GMM") {
    se_trained <- run_se_GMM_py(y, x, w)
  } else {
    se_trained <- run_se_ML_py(y, x, w)
  }

  beta <-
    data.frame(
      var = c("Intercept", x_var_names, "lambda"),
      coef = se_trained$betas,
      se = se_trained$std_err
    ) %>%
    data.table()

  return(list(beta = beta, reg_formula = reg_formula, y_reg_var = y_reg_var, t_reg_var = t_reg_var))
}

# library(tidyverse)
# library(reticulate)
# library(data.table)
# library(spdep)

# source_python(here::here("Codes/functions/python/run_se_ML.py"))
# source_python(here::here("Codes/functions/python/run_se_GMM.py"))

# sim_data <- readRDS(here::here("Data/sim_data.rds"))
# num_fields <- 5
# data <- sim_data$reg_data[[1]]$data[1:num_fields] %>% rbindlist()
# w <-
#   sim_data$weights_matrix[[1]]$Wls_50 %>%
#   listw2mat() %>%
#   kronecker(diag(num_fields), .)

# run_sem(yield ~ N + I(N^2), w, data, method = "GMM")

find_eonr <- function(sim_results, p_crop, p_N) {
  eonr <-
    sim_results$yield_eval %>%
    .[, pi_hat := p_crop * y_hat - p_N * N] %>%
    .[, .SD[which.max(pi_hat), ], by = aunit_id] %>%
    .[, .(aunit_id, N)] %>%
    setnames("N", "opt_N_hat")
}
