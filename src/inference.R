
################################################################################
#
# Purpose: New HMM method for behavioural state identification from high frequency movement data
#
# Author:  ..
#
# Date created: December 02 2025
#
################################################################################



# Install missing packages
for (pkg in  c("Rcpp", "dplyr", "accelerometry")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
  
}

sourceCpp("src/forward_algorithm.cpp")



#####################################################
#            Log Likelohood Calculation
#####################################################


# log likelihood
log_emission <- function(data, params, j,  dis,  variable) {
  
  

  # von Mises with dynamic variable name
  log_vonmises <- function(x, params, var, j) {
    kappa_name <- paste0("kappas_", var)
    dvm_cpp(x, params[[kappa_name]][j])
  }
  
  # Soft truncated von Mises
  log_tvonmises <- function(x, params, var, j) {
    kappa_name <- paste0("kappas_", var)
    w_name <- paste0("w_", var)
    k =  min(params[[kappa_name]][j],100) # 
    w =  min(params[[w_name]][j],100)
    dsoft_tvm_cpp(x, params[[kappa_name]][j], params[[w_name]][j])
  }
  
  
  # Gamma
  log_gamma <- function(x, params, var, j) {
    shape_name <- paste0("shapes_", var)
    rate_name <- paste0("rates_", var)
    dgamma_cpp(x, shape = params[[shape_name]][j], rate = params[[rate_name]][j])
    #dgamma(x , shape = params[[shape_name]][j], rate = params[[rate_name]][j] , log = TRUE)
  }
  
  # Exponential
  log_exp <- function(x, params, var, j) {
    lambda_name <- paste0("lambdas_", var)
    dexp_cpp(x, lambda = params[[lambda_name]][j])
  }
  
  
  
  lpdf_dict <- list(
    "tvonmises" = log_tvonmises,
    "gamma"     = log_gamma,
    "exp"       = log_exp,
    "vonmises"  = log_vonmises
  )
  
  
  total <- 0

  for (i in seq_along(variable)) {
    varname <- variable[i]
    distfun <-   lpdf_dict[[dis[[i]]]]
    total <- total + distfun(data[[varname]], params, varname , j)
  }
  return(total)
  
}



#####################################################
#                 Simulation
#####################################################

simulate_trajectory <- function( n , m , params , variables , dis_names ,  covariate = NULL , formula = NULL, seed = NULL ){
  


  simulate_gamma <- function(n, params, var, j, seed = NULL) {
    shape_name <- paste0("shapes_", var)
    rate_name  <- paste0("rates_", var)
    rgamma(n, shape = params[[shape_name]][j], rate = params[[rate_name]][j])
  }
  
  
  simulate_exp <- function(n, params, var, j, seed = NULL) {
    lambda_name <- paste0("lambdas_", var)
    rexp(n, rate = params[[lambda_name]][j])
  }
  
  
  simulate_vm <- function(n, params, var, j, seed = NULL) {
    kappa_name <- paste0("kappas_", var)
    (Rfast::rvonmises(n, m = 0, k = params[[kappa_name]][j]) + pi) %% (2 * pi) - pi 
  }
  
  

  simulate_tvm <- function(n, params, var, j, seed = NULL) {
    
    cdf_tvonmises <- function(k, w, n = 100000) {
      x_vals <- seq(-pi, pi, length.out = n)
      dx <- x_vals[2] - x_vals[1]
      pdf_vals <-  exp(dsoft_tvm_cpp(x_vals, k, w))
      cdf_vals <- cumsum(pdf_vals) * dx
      cdf_vals <- cdf_vals / max(cdf_vals)
      list(x = x_vals, cdf = cdf_vals)
    }
    
    
    kappa_name <- paste0("kappas_", var)
    w_name     <- paste0("w_", var)
    cdf_obj <- cdf_tvonmises(params[[kappa_name]][j], params[[w_name]][j])
    inv_cdf <- approxfun(cdf_obj$cdf, cdf_obj$x, yleft = -pi, yright = pi)
    seed <- ifelse(!is.null(seed), seed, sample(1:1e6, 1))
    set.seed(seed)
    u_vals <- runif(n)
    inv_cdf(u_vals)
  }

  
  sim_dict <- list(
    "tvonmises" = simulate_tvm ,
    "gamma"     = simulate_gamma,
    "exp"       = simulate_exp,
    "vonmises"  = simulate_vm
  )
  
  
  m_matrix <- cov.matrix(n , covariate , formula )
  trans_p <- exp(log_transition_probabilities( m_matrix , Filter( is.matrix , params) ))
  seed = ifelse( !is.null(seed)  , seed , sample(1:1e5,1) )
  states <- numeric(n)
  set.seed(seed)
  states[1] <- sample(1:m ,1, prob = exp(params$ini_p)/sum(exp(params$ini_p)) )
  for (i in 2:n) {
    seed = seed+1
    set.seed(seed)
    # Simulate next state based on the current state
    current_state <- states[i - 1]
    states[i] <- sample(1:m, size = 1, prob = trans_p[,,i][current_state, ])
  }

  out <- as.data.frame(matrix(NA, nrow = N, ncol = length(variables)))
  colnames(out) <- variables
  out = cbind(out,  data.frame(Step = 1:N, states = states))
  
  for (i in  unique(states) ) {
    state_indices <- which(states == i)
    for(v in 1:length(variables)){
      out[[variables[v]]][state_indices] = sim_dict[dis_names][[v]](length(state_indices) , params , variables[v], i , seed)
    }
  }
  
  out
  

}



#################################################################
#                        Initialisation
#################################################################


mle_emission_params <- function(data, states, dis_names, variables , m ){
  
  
  mle_gamma <- function(data, var) {
    x <- data[[var]]
    vals <- c(mean(x)^2 / var(x), mean(x) / var(x))
    names(vals) <- c(paste0("shapes_", var), paste0("rates_", var))
    as.list(vals)
  }
  
  mle_exp <- function(data, var) {
    x <- data[[var]]
    val <- 1 / mean(x)
    setNames(list(val), paste0("lambdas_", var))
  }
  
  
  mle_vonmises <- function(data, var) {
    x <- data[[var]]
    est <- Rfast::vm.mle(x - mean(x))$param[2]
    setNames(list(est), paste0("kappas_", var))
  }
  
  
  
  mle_tvonmises <- function(data, var) {
    x <- na.omit(data[[var]])
    neg_loglik <- function(par, x) {
      k <- pmin(exp(par[1]),100)
      w <- pmin(exp(par[2]),100)
      -sum(dsoft_tvm_cpp(x, k, w))
    }
    
    best_fit <- NULL
    best_value <- Inf
    ini_vals_k <- exp(seq( -0.1 , 4 , 1))
    ini_vals_w <- exp(seq( -0.1 , 4 , 1))
    
    
    # try many starting values 
    for (ini_k in ini_vals_k) {
      for (ini_w in ini_vals_w) {
        init_params <- log(c(ini_k, ini_w))
        
        fit <- tryCatch(
          optim(par = init_params, fn = neg_loglik, x = x, control = list(maxit = 100000000), ),
          error = function(e) NULL
        )
        
        if (!is.null(fit) && fit$value < best_value) {
          best_fit <- fit
          best_value <- fit$value
        }
      }
    }
    
    
    k <- exp(best_fit$par[1])
    w <- exp(best_fit$par[2])
    
    result <- list()
    result[[paste0("kappas_", var)]] <- k
    result[[paste0("w_", var)]] <- w
    return(result)
  }
  
  
  estimation_dict <- list(
    "tvonmises" = mle_tvonmises,
    "gamma"     = mle_gamma,
    "exp"       = mle_exp,
    "vonmises"  = mle_vonmises
  )
  
  
  
  all_state <- sort(unique(states))
  params = build_param_template_from_dict(dis_names , variables , m , "random" )
  for (state in all_state) {
    sub <- na.omit(data[states == state, ])
    
    if(nrow(sub)!=0){
      for (i in seq_along(dis_names)) {
        dist_name <- dis_names[i]
        var_name <- variables[i]
        est_fun <- estimation_dict[[dist_name]]
        
        est <- est_fun(sub, var_name)
        
        for (pname in names(est)) {
          params[[pname]][state] <- est[[pname]]
        }
      }
    }
  }
  
  
  return(params)
}







initialize_based_on_kmean <- function(data, m, w, dis_names ,  variables ,
                                      covariate = NULL , formula = NULL ) {
  
  
  
  
  
  # Helpers
  rolling_stats <- function(data, window = 5) {
    out <- list()
    for (col in names(data)) {
      x <- data[[col]]
      if (is.numeric(x)) {
        m_mean_acc  <-  forecast::ma(x,window)
        out[[paste0(col, "_mean")]] <- m_mean_acc
        out[[paste0(col, "_var")]]  <- sqrt(forecast::ma(x^2,window) - m_mean_acc^2)
      }
    }
    return(  as.data.frame(out)  )
  }
  
  
  calculate_transition_matrix <- function(cluster, covariate , formula , m ){
    
    m_matrix <- cov.matrix( length(cluster), covariate , formula)
    res <- setNames(
      lapply(seq_len(ncol(m_matrix)), function(x) matrix(0, m, m)),
      paste0("t", seq_len(ncol(m_matrix)), "p")
    )
    
    trans.p <- matrix(0.00001 , m , m)
    for (i in 1:(length(cluster)-1)) {
      current_state <- cluster[i]
      next_state <- cluster[i + 1]
      trans.p[current_state, next_state] <- trans.p[current_state, next_state] + 1
    }
    trans.p <- trans.p/rowSums(trans.p)
    
    X <-  log(trans.p/( 1- trans.p ) )
    X[X ==  Inf]  <- 50  ;  X[X == -Inf] <- -50
    res[[1]] <- X
    
    return(res) 
  }
  
  
  
  calculate_regression_params <- function(covariate, cluster, formula, m) {
    
    m_matrix <- model.matrix(formula, covariate)
    res <- setNames(lapply(seq_len(ncol(m_matrix)), function(x) matrix(0, m, m)), paste0("t", seq_len(ncol(m_matrix)), "p"))
    States <- cbind( cluster , next_clust = c(cluster[-1], NA) )
    covariate <- na.omit(cbind(covariate, States))
    new_formula <- as.formula(paste("next_clust ~", as.character(formula)[2]))
    all_possible <- expand.grid(cluster = seq_len(m), next_clust = seq_len(m))
    weights = c(rep(1, nrow(covariate)), rep(0.1, 2*nrow(all_possible)))
    covariate <- bind_rows(covariate, all_possible , all_possible)
    covariate$weights <- weights
    for (i in seq_len(m)) {
      subsets <- covariate %>% filter(cluster == i)  %>% 
        mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))  
      fit <- nnet::multinom(new_formula, data = subsets, trace = FALSE, weights = subsets$weights)
      cof <- as.matrix(coef(fit))
      cof <- if (m == 2) t(cof) else cof
      res[seq_len(ncol(cof))] <- lapply(seq_len(ncol(cof)), function(z) { res[[z]][i, ] <- c(0, cof[, z]); res[[z]] })
    }
    for (i in seq_along(res)) {
      res[[i]] <- res[[i]] - res[[i]][, ncol(res[[i]])]
    }
    return(res)
  }
  #####################
  
  
  
  m_matrix <- cov.matrix(nrow(data) , covariate,formula)
  fdata <- data
  data <- data[, variables, drop = FALSE]
  if ("log_vonmises" %in% dis_names | "log_tvonmises" %in% dis_names) {
    data$cosx <- cos(data$x)
    data$sinx <- sin(data$x)
    data <- subset(data, select = -x)
  }
  
  
  data_rolling = rolling_stats(data, w)
  na_cols <- !complete.cases(data_rolling)
  cluster <- rep(NA, nrow(data))
  kmeans_result <- kmeans(scale(data_rolling[!na_cols, ]), centers = m, nstart =  round(0.25*nrow(data)) )
  cluster[!na_cols] <- as.numeric(kmeans_result$cluster)
  
  
  emission <-  mle_emission_params(fdata[!na_cols,] , cluster[!na_cols] , dis_names , variables , m )
  
  if(is.null(formula)){
    trans_p <- calculate_transition_matrix(cluster , covariate , formula , m)
  }else{
    trans_p <- calculate_regression_params(covariate , cluster , formula , m)
  }

  
  trans <- exp(log_transition_probabilities(as.matrix(m_matrix[1,]) , trans_p ))
  eig <- eigen(t(trans[,,1]))
  init_dist <- Re(eig$vectors[, 1])
  init_dist <- init_dist / sum(init_dist)
  epsilon <- 1e-10
  init_dist[init_dist < epsilon] <- epsilon
  init_dist[init_dist > 1 - epsilon] <- 1 - epsilon
  use_logit <- TRUE
  
  ini_p <- if (use_logit) {
    (log(init_dist) - log(1 - init_dist))
  } else {
    (init_dist)
  }
  
  
  res <- c(emission , trans_p)
  res$ini_p = ini_p
  
  return(res)
}











#############################################################################
#                      template
#############################################################################


build_param_template_from_dict <- function(dis_names, variable, N, fill = NA) {
  stopifnot(length(dis_names) == length(variable))
  
  param_dict <- list(
    "tvonmises" = c("kappas", "w"),
    "gamma"    = c("shapes", "rates"),
    "exp"       = c("lambdas"),
    "vonmises"  = c("kappas")
  )
  
  params <- list()
  
  for (i in seq_along(dis_names)) {
    dist <- dis_names[i]
    var  <- variable[i]
    param_names <- param_dict[[dist]]
    
    for (p in param_names) {
      full_name <- paste0(p, "_", var)
      
      if (is.character(fill) && fill == "random") {
        params[[full_name]] <- runif(N, 0.0001, 20)  # uniform(0,1)
      } else {
        params[[full_name]] <- rep(fill, N)
      }
    }
  }
  
  return(params)
}




##########################################################################
#  Fit HMM
##########################################################################






# Prepare Covariate Model Matrix

cov.matrix <- function(n , covariate = NULL,formula = NULL){
  if( is.null(covariate) ){
    m_matrix <- model.matrix(~1, data.frame(x = rep(1,n) ))
  }else{
    m_matrix <- model.matrix(formula, covariate)
  }
  m_matrix
}


# Check Data
check_data <- function(data , variables ,  covariate = NULL , formula = NULL){
  # Check movement matrix
  true_names <- variables
  if(!is.data.frame(data)) stop("data must be dataframe")
  provided_names <- names(data)
  if (any(!true_names %in% provided_names)) {
    missing_columns <- setdiff(true_names , provided_names)
    stop(sprintf(
      "The given dataframe is not suitable. The following required columns are missing: %s",
      paste(missing_columns, collapse = ", ")
    ))
  }
  
  # Check covariate matrix
  if(!is.null(covariate)){
    if(nrow(covariate) != nrow(data)) stop("Covariate and data must have the same length")
    if(!is.data.frame(covariate)) stop("covariate must be dataframe")
    if( is.null(formula) ){
      stop("formula is missed")
    }else{
      formula_vars <- all.vars(formula)
      data_vars <- colnames(covariate)
      missing_vars <- setdiff(formula_vars, data_vars)
      if (length(missing_vars) != 0) {
        stop(sprintf(
          "The following variables are missing from the dataset: %s",
          paste(missing_vars, collapse = ", ")
        ))
      }
    }
  }
}




# Viterbi Algorithm
viterbi <- function(params , data , covariate , formula , dis,  variable , m ) {
  
  m_matrix <- cov.matrix(nrow(data) , covariate ,formula )
  
  n <- nrow(m_matrix)
  
  logA <- log_transition_probabilities(m_matrix, Filter(is.matrix, params))
  logpi <- params$ini_p - log(sum(exp(params$ini_p)))
  logB <- sapply(1:m, function(i) {
    log_emission(data ,  params , i,   dis,  variable )
  })
  
  states_decoded <- viterbi_cpp(logA, logpi , logB)
}


# Optimise likelihood 
fit <- function( data , covariate , formula ,   dis,  variable  , initial_params , m, maxit = 1e50 ) {
  
  
  transform_to_realline <- function(params,m) {
    name <- names(params)
    for (i in seq_along(params)) {
      vec <- params[[i]]
      nm <- name[i]
      if (length(vec) == m && !nm %in% c("ini_p")) {
        params[[i]] <- log(vec)
      } else if (length(vec) > m) {
        params[[i]] <- vec - vec[, m]
        params[[i]] <- params[[i]][, 1:(m - 1)]
      }
    }
    params
  }
  
  
  # A function to take transform optim output
  params_unflatten <- function(data, flat ,group_names, m) {
    reversed <- split(flat, group_names)
    for (nm in names(reversed)) {
      if (grepl("^t[0-9]+p$", nm)) {
        reversed[[nm]] <- cbind(matrix(reversed[[nm]], ncol = m - 1), 0)
      } else if (length(reversed[[nm]]) == m && !nm %in% c("ini_p")) {
        reversed[[nm]] <- exp(reversed[[nm]])
      }
    }
    reversed
  }
  
  
  
  # objective Function
  objective_function <- function(flat , data , m_matrix ,  dis,  variable , group_names, m) {

    params <- params_unflatten(data,flat , group_names,m)
    B <- sapply(seq_along(params$ini_p), function(i) {
      log_emission(data , params , i , dis , variable)})
    tr_prob <- log_transition_probabilities( m_matrix, Filter(is.matrix, params))
    ini_p <- params$ini_p - log( sum(exp(params$ini_p)) )
    likelihood <- forward_algorithm(tr_prob, ini_p, B)
    -likelihood$log_likelihood
  }
  
  
  

  params <- transform_to_realline(initial_params, m)
  flat <- unlist(params)
  group_names <- gsub("\\d+$", "", names(flat))
  m_matrix <- cov.matrix(nrow(data),covariate ,formula )
  
  

  
  
  res <- optim(par = flat, fn = objective_function ,
               data = data ,
               m_matrix = m_matrix ,
               dis = dis,
               variable=variable,
               group_names = group_names,
               control = list(maxit = maxit) ,m =m)
  
  params <- params_unflatten(data,res$par, group_names,m)
  params$neg_lik_value <- res$value
  params$convergence <- res$convergence
  params$counts <- res$counts
  params$n_params <- length(flat)
  return(params)
}




fit_HMM_with_kmean_initialsation <- function(data, m, covariate, formula, dis_names, variables, w_values, maxit = 1e10) {
  check_data(data, variables, covariate, formula)
  
  if("y" %in% c(variables) ){
    meanY <- mean(data$y)
    data$y <- data$y/meanY
  }
  
  if("z" %in% c(variables) ){
    meanZ <- mean(data$z)
    data$z <- data$z/meanZ
  }
  
  # Step 1: Fit models with different w initializations
  starting_values <- lapply(w_values, function(w) {
    tryCatch({
      initial_params <- initialize_based_on_kmean(data, m, w, dis_names, variables, covariate, formula)
      est_params <- fit(data, covariate, formula, dis_names, variables, initial_params, m, maxit)
      return(est_params)
    }, error = function(e) {
      message("Error encountered: ", e$message)
      return(NULL)
    })
  })
  
  
  indx <- which.min(sapply(starting_values, function(x) x$neg_lik_value))
  model <- starting_values[[indx]]
  #model$states <- viterbi(model, data, covariate, formula, dis_names, variables, m)
  #model <- match_model_params(1:m , model)
  
  if("y" %in% c(variables) ){
    ydis <- dis_names[which(variables=="y")]
    if(ydis == "gamma" || ydis == "exp"){
      model$rates_y <- model$rates_y/meanY
    }
  }
  
  if("z" %in% c(variables) ){
    zdis <- dis_names[which(variables=="z")]
    if(zdis == "gamma" || zdis == "exp" ){
      model$rates_z <- model$rates_z/meanZ
    }
  }
  
  model <- lapply(model, function(x) { names(x) <- NULL; x })
  model$history <- starting_values
  return(model)
}



