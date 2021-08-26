# helpers for method evaluation

# Compute the relative L1 difference between fitted and true
# p tensors
evaluate_fit <- function(fit_decomp, p, X_array, 
                         R = length(fit_decomp$lambda)) {
  Yhat <- get_Yhat(fit_decomp, R = R)
  phat <- apply(Yhat, c(2, 3), function(x) {
    x <- x - max(x)
    return(exp(x) / sum(exp(x)))
  })
  
  # Relative L1 difference
  mat_L1_relative <- matrix(NA, nrow = dim(phat)[2], ncol = dim(phat)[3])
  for(j in seq(1, nrow(mat_L1_relative)))
    for(k in seq(1, ncol(mat_L1_relative)))
      mat_L1_relative[j, k] <- mean(abs(phat[, j, k] - p[, j, k]) / p[, j, k])
  
  results <- list(mat_L1_relative = mat_L1_relative)
  
  return(results)
}

# estimate overinflation from the fitted model
estimate_overInflation <- function(fit_decomp, p, X_array, R, 
                                   weighted = FALSE) {
  Yhat <- get_Yhat(fit_decomp, R = R)
  phat <- apply(Yhat, c(2, 3), function(x) {
    x <- x - max(x)
    return(exp(x) / sum(exp(x)))
  })
  
  Xsum_m <- apply(X_array, c(2, 3), sum)
  
  Xhat <- vapply(seq_len(dim(X_array)[1]),
                 function(i) phat[i, , ] * Xsum_m,
                 matrix(0.0, nrow = dim(X_array)[2], ncol = dim(X_array)[3]))
  Xhat <- aperm(Xhat, perm = c(3, 1, 2))
  var_obs <- apply((X_array - Xhat)^2, c(2, 3), sum)
  var_exp <- apply(Xhat * (1 - phat), c(2, 3), sum)
  if(weighted)
    var_exp <- var_exp * (1 + (Xsum_m - 1) * fit_decomp$wt_estimate$phi)
  over_inflation <- var_obs / var_exp
  
  return(over_inflation)
}

# create feature/subject/time/sample loadings for
# microTensor/CTF/PCA fit
create_loading <- function(fit_decomp,
                           feature_names,
                           subject_names,
                           time_names,
                           class = "feature") {
  lambdas <- fit_decomp$lambda
  R <- length(lambdas)
  PCA <- !("t" %in% names(fit_decomp))
  
  if(is.null(feature_names))
    subject_names <- paste0("Feature", seq(1, nrow(fit_decomp$m)))
  if(is.null(subject_names))
    subject_names <- paste0("Subject", seq(1, nrow(fit_decomp$s)))
  if(is.null(time_names))
    time_names <- paste0("Time", seq(1, nrow(fit_decomp$t)))
  
  if(class == "feature") {
    tb_loading <- cbind(data.frame(Feature = feature_names),
                        fit_decomp$m)
    colnames(tb_loading) <- c("Feature", paste0("Axis ", seq(1, R)))
  } else {
    if(PCA) {
      mat_s <- fit_decomp$s_full
      tb_loading <- seq(1, R) %>% 
        purrr::map(function(r) {
          r_tb <- tidyr::expand_grid(time_names,
                                     subject_names) %>% 
            dplyr::mutate(Loading = mat_s[, r] * lambdas[r])
          colnames(r_tb) <- c("Time", "Subject", paste0("Axis ", r))
          return(r_tb)
        }) %>% 
        purrr::reduce(dplyr::left_join, by = c("Subject", "Time"))
      
      tb_loading_nolambda <- seq(1, R) %>% 
        purrr::map(function(r) {
          r_tb <- tidyr::expand_grid(time_names,
                                     subject_names) %>% 
            dplyr::mutate(Loading = mat_s[, r])
          colnames(r_tb) <- c("Time", "Subject", paste0("Axis ", r))
          return(r_tb)
        }) %>% 
        purrr::reduce(dplyr::left_join, by = c("Subject", "Time"))
    } 
    
    if(class == "sample") {
      if(!PCA) {
        mat_s <- fit_decomp$s
        mat_t <- fit_decomp$t
        
        tb_loading <- seq(1, R) %>% 
          purrr::map(function(r) {
            r_tb <- tidyr::expand_grid(subject_names,
                                       time_names) %>% 
              dplyr::mutate(Loading = as.vector(outer(mat_t[, r],
                                                      mat_s[, r]) * lambdas[r]))
            colnames(r_tb) <- c("Subject", "Time", paste0("Axis ", r))
            return(r_tb)
          }) %>% 
          purrr::reduce(dplyr::left_join, by = c("Subject", "Time"))
      }
    }
    
    if(class == "subject") {
      if(!PCA) {
        mat_s <- fit_decomp$s
        tb_loading <- cbind(subject_names, as.data.frame(t(t(mat_s) * lambdas)))
      } else {
        tb_loading <- tb_loading %>% 
          tidyr::pivot_longer(cols = paste0("Axis ", seq(1, R)),
                              values_to = "Loading",
                              names_to = "Axis") %>% 
          dplyr::group_by(Subject, Axis) %>% 
          dplyr::summarize(Loading = mean(Loading, na.rm = TRUE)) %>% 
          dplyr::ungroup() %>% 
          tidyr::pivot_wider(names_from = Axis, values_from = Loading)
      }
      
      colnames(tb_loading) <- c("Subject", paste0("Axis ", seq(1, R)))
    }
    
    if(class == "time") {
      if(!PCA) {
        mat_t <- fit_decomp$t
        tb_loading <- cbind(time_names, as.data.frame(mat_t))
      } else {
        tb_loading <- tb_loading_nolambda %>% 
          tidyr::pivot_longer(cols = paste0("Axis ", seq(1, R)),
                              values_to = "Loading",
                              names_to = "Axis") %>% 
          dplyr::group_by(Time, Axis) %>% 
          dplyr::summarize(Loading = mean(Loading, na.rm = TRUE)) %>% 
          dplyr::ungroup() %>% 
          tidyr::pivot_wider(names_from = Axis, values_from = Loading)
      }
      
      colnames(tb_loading) <- c("Time", paste0("Axis ", seq(1, R)))
    }
  }
  return(tb_loading)
}

# generate ggplot like color hues
gg_clor_hue <- function(values = NULL, n = NULL) 
{
  if (!is.null(values)) {
    if (anyDuplicated(values)) 
      stop("Values should be unique if provided!")
    if (!is.character(values)) 
      stop("Values should be of character class!")
    n <- length(values)
    hues <- seq(15, 375, length = n + 1)
    colors <- hcl(h = hues, l = 65, c = 100)[1:n]
    names(colors) <- values
    return(colors)
  }
  hues <- seq(15, 375, length = n + 1)
  colors <- hcl(h = hues, l = 65, c = 100)[1:n]
}
