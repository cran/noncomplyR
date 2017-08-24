# #' Column Names
# #'
# #' \code{name_columns} gives the appropriate column names for the output from \code{compliance_chain}.
# #'
# #' @param outcome_model character string indicating the outcome model used
# #' @param exclusion_restriction logical value indicating whether the Exclusion Restriction assumption was made.
# #' @param strong_access logical value indicating whether the Strong Access Monotonicity assumption was made.
# #'
# #' @keywords internal
# #' @examples
# #' name_columns(outcome_model = "binary", exclusion_restriction = TRUE, strong_access = TRUE)
# #' name_columns(outcome_model = "normal", exclusion_restriction = FALSE, strong_access = FALSE)
# #' @return a character string containing the names of the model parameters
name_columns <- function(outcome_model, exclusion_restriction, strong_access){

  compliance_params <- c("omega_c", "omega_n", "omega_a")

  if(outcome_model == "binary"){
    binary_params <- c("p_c0", "p_c1", "p_n0", "p_n1", "p_a0", "p_a1")

    if(exclusion_restriction & strong_access) param_names <- c(compliance_params[1:2], binary_params[1:2], "p_n")
    if(exclusion_restriction & !strong_access) param_names <- c(compliance_params, binary_params[1:2], "p_n", "p_a")
    if(!exclusion_restriction & strong_access) param_names <- c(compliance_params[1:2], binary_params[1:4])
    if(!exclusion_restriction & !strong_access) param_names <- c(compliance_params, binary_params)


  }
  if(outcome_model == "normal"){
    normal_params <- c("mu_c0", "sigma2_c0", "mu_c1", "sigma2_c1", "mu_n0", "sigma2_n0", "mu_n1", "sigma2_n1", "mu_a0",
                     "sigma2_a0", "mu_a1", "sigma2_a1")

    if(exclusion_restriction & strong_access) param_names <- c(compliance_params[1:2], normal_params[1:4], "mu_n", "sigma2_n")
    if(exclusion_restriction & !strong_access) param_names <- c(compliance_params, normal_params[1:4], "mu_n", "sigma2_n", "mu_a", "sigma2_a")
    if(!exclusion_restriction & strong_access) param_names <- c(compliance_params[1:2], normal_params[1:8])
    if(!exclusion_restriction & !strong_access) param_names <- c(compliance_params, normal_params)


  }

  param_names


}

# #' Normal Outcome Model
# #'
# #' \code{g_normal} is a wrapper function for computing the complete-data outcome
# #' likelihood for the Normal model
# #'
# #' @param y a numeric vector representing the outcome values at which the likelihood should be computed
# #' @param param_vals a numeric vector containing the parameter values at which the likelihood should be computed
# #' @examples g_normal(5, c(4, 2)); g_normal(rnorm(10), c(0, 1))
# #' @keywords internal
# #'
g_normal <- function(y, param_vals){

  stats::dnorm(y, mean = param_vals[1], sd = sqrt(param_vals[2]))

}

# #' Binary Outcome Model
# #'
# #' \code{g_binary} is a wrapper function for computing the complete-data outcome likelihood for
# #' the Binary outcome model
# #'
# #' @param y a numeric vector representing the outcome values at which the likelihood should be computed
# #' @param param_vals a numeric vector containing the parameter value at which the likelihood should be computed
# #' @examples g_binary(0, .4); g_binary(1, .8)
# #' @keywords internal
g_binary <- function(y, param_vals){

  stats::dbinom(x = y, size = 1, prob = param_vals)


}

# #' Compute Posterior Hyperparameters for Normal-Inverse Gamma Distribution
# #'
# #' \code{posterior_params_ninvgamma}
# #'
# #' @param y a numeric vector containing the outcomes conditioned on in the posterior distribution.
# #' @param hyper_parameters a numeric vector containing the prior hyperparameters.
# #' @examples posterior_params_ninvgamma(rnorm(10), rep(1, 4))
# #' @keywords internal
# #' @return a numeric vector containing the updated hyperparameters.
# #'
posterior_params_ninvgamma <- function(y, hyper_parameters){

  n <- length(y)

  sum_y <- sum(y)
  sum_y2 <- sum(y^2)

  tau_bar <- 1/(1/(hyper_parameters[2]) + n)

  theta_bar <- (hyper_parameters[1]/hyper_parameters[2] + sum_y)*tau_bar

  a_bar <- hyper_parameters[3] + n/2

  b_bar <- hyper_parameters[4] + .5*(hyper_parameters[1]^2/hyper_parameters[2] + sum_y2 - theta_bar^2/tau_bar)

  c(theta_bar, tau_bar, a_bar, b_bar)

}


# #' Impute Compliance Types with Binary Outcome Model
# #'
# #' \code{impute_c_binary} imputes the missing compliance type variable for subjects under
# #' the binary outcome model
# #' @param dat dataframe containing outcome, treatment assignment, and treatment indicator.
# #' @param theta current parameter value.
# #' @param exclusion_restriction a logical value indicating whether the Exclusion Restriction assumption holds.
# #' @param strong_access a logical value indicating whether the Strong Access Monotonicity assumption holds.
# #' @examples
# #' dat <- data.frame(y = rep(c(0, 1), 5), z = rep(c(0, 1), each = 5), d = c(rep(0, 6), rep(1, 4)) )
# #' theta_k <- c(.9, .1, .5, .5, .5, .5)
# #' impute_c_binary(dat, theta = theta_k, exclusion_restriction = FALSE, strong_access = TRUE)
# #' @keywords internal
# #' @return a character vector of imputed compliance types. "a" = Always Taker, "n" = Never Taker, "c" = Complier
# #'
impute_c_binary <- function(dat, theta, exclusion_restriction = T, strong_access = T){

  output <- rep(NA, dim(dat)[1])

  if(!strong_access){

  omega_c <- theta[1]

  omega_n <- theta[2]

  omega_a <- theta[3]



  if(!exclusion_restriction){



  c_or_n <- dat[dat$z == 0 & dat$d == 0, ]
  output[dat$z == 0 & dat$d == 1] <- "a"
  output[dat$z == 1 & dat$d == 0] <- "n"
  c_or_a <- dat[dat$z == 1 & dat$d == 1, ]



  g_c0 <- g_binary(c_or_n$y, theta[4])
  g_c1 <- g_binary(c_or_a$y, theta[5])
  g_n0 <- g_binary(c_or_n$y, theta[6])
  g_n1 <- g_binary(c_or_n$y, theta[7])
  g_a0 <- g_binary(c_or_a$y, theta[8])
  g_a1 <- g_binary(c_or_a$y, theta[9])


  output[dat$z == 0 & dat$d == 0] <- ifelse(stats::rbinom(dim(c_or_n)[1], 1, omega_c*g_c0/(omega_c*g_c0 + omega_n*g_n0)), "c", "n")

  output[dat$z == 1 & dat$d == 1] <- ifelse(stats::rbinom(dim(c_or_a)[1], 1, omega_c*g_c1/(omega_c*g_c1 + omega_a*g_a1)), "c", "a")



  }

  if(exclusion_restriction){

    c_or_n <- dat[dat$z == 0 & dat$d == 0, ]
    c_or_a <- dat[dat$z == 1 & dat$d == 1, ]
    output[dat$z == 1 & dat$d == 0] <- "n"
    output[dat$z == 0 & dat$d == 1] <- "a"


    g_c0 <- g_binary(c_or_n$y, theta[4])
    g_n0 <- g_binary(c_or_n$y, theta[6])
    g_a1 <- g_binary(c_or_a$y, theta[7])
    g_c1 <- g_binary(c_or_a$y, theta[5])

    output[dat$z == 0 & dat$d == 0] <- ifelse(stats::rbinom(dim(c_or_n)[1], 1, omega_c*g_c0/(omega_c*g_c0 + omega_n*g_n0)), "c", "n")
    output[dat$z == 1 & dat$d == 1] <- ifelse(stats::rbinom(dim(c_or_a)[1], 1, omega_c*g_c1/(omega_c*g_c1 + omega_a*g_a1)), "c", "a")

  }


  return(output)

  }

  if(strong_access){

    omega_c <- theta[1]
    omega_n <- theta[2]

    if(!exclusion_restriction){



      c_or_n <- dat[dat$z == 0, ]
      output[dat$z == 1 & dat$d == 0] <- "n"
      output[dat$z == 1 & dat$d == 1] <- "c"

      g_c0 <- g_binary(c_or_n$y, theta[3])
      g_n0 <- g_binary(c_or_n$y, theta[5])


      output[dat$z == 0] <- ifelse(stats::rbinom(dim(c_or_n)[1], 1, omega_c*g_c0/(omega_c*g_c0 + omega_n*g_n0)), "c", "n")

    }

    if(exclusion_restriction){

      c_or_n <- dat[dat$z == 0, ]

      output[dat$z == 1 & dat$d == 0] <- "n"
      output[dat$z == 1 & dat$d == 1] <- "c"

      g_c0 <- g_binary(c_or_n$y, theta[3])
      g_n0 <- g_binary(c_or_n$y, theta[5])


      output[dat$z == 0] <- ifelse(stats::rbinom(dim(c_or_n)[1], 1, omega_c*g_c0/(omega_c*g_c0 + omega_n*g_n0)), "c", "n")

    }



    return(output)
  }
}

# #' Impute Compliance Types with Normal Outcome Model
# #'
# #' \code{impute_c_normal} imputes the missing compliance type variable for subjects under
# #' the Normal outcome model
# #' @param dat dataframe containing outcome, treatment assignment, and treatment indicator
# #' @param theta current parameter value
# #' @param exclusion_restriction a logical value indicating whether the Exclusion Restriction assumption
# #' should be made
# #' @param strong_access a logical value indicating whether the Strong Access Monotonicity assumption
# #' should be made
# #' @keywords internal
# #'
# #' @examples
# #' dat <- data.frame(y = rnorm(10), z = rep(c(0, 1), each = 5), d = c(rep(0, 7), 1, 1, 1))
# #' theta_k <- c(.5, .5, 0, 1, 0, 1, 0, 1)
# #' impute_c_normal(dat, theta = theta_k, exclusion_restriction = TRUE, strong_access = TRUE)
# #'
# #' @return a character vector of imputed compliance types. "a" = Always Taker, "n" = Never Taker, "c" = Complier
# #'
impute_c_normal <- function(dat, theta, exclusion_restriction = T, strong_access = T){

  output <- rep(NA, dim(dat)[1])

  if(!strong_access){

  omega_c <- theta[1]

  omega_n <- theta[2]

  omega_a <- theta[3]

  if(!exclusion_restriction){



    c_or_n <- dat[dat$z == 0 & dat$d == 0, ]
    output[dat$z == 0 & dat$d == 1] <- "a"
    output[dat$z == 1 & dat$d == 0] <- "n"
    c_or_a <- dat[dat$z == 1 & dat$d == 1, ]



    g_c0 <- g_normal(c_or_n$y, theta[4:5])
    g_c1 <- g_normal(c_or_a$y, theta[6:7])
    g_n0 <- g_normal(c_or_n$y, theta[8:9])
    g_n1 <- g_normal(c_or_n$y, theta[10:11])
    g_a0 <- g_normal(c_or_a$y, theta[12:13])
    g_a1 <- g_normal(c_or_a$y, theta[14:15])


    output[dat$z == 0 & dat$d == 0] <- ifelse(stats::rbinom(dim(c_or_n)[1], 1, omega_c*g_c0/(omega_c*g_c0 + omega_n*g_n0)), "c", "n")

    output[dat$z == 1 & dat$d == 1] <- ifelse(stats::rbinom(dim(c_or_a)[1], 1, omega_c*g_c1/(omega_c*g_c1 + omega_a*g_a1)), "c", "a")



  }

  if(exclusion_restriction){

    c_or_n <- dat[dat$z == 0 & dat$d == 0, ]
    c_or_a <- dat[dat$z == 1 & dat$d == 1, ]
    output[dat$z == 1 & dat$d == 0] <- "n"
    output[dat$z == 0 & dat$d == 1] <- "a"


    g_c0 <- g_normal(c_or_n$y, theta[4:5])
    g_n0 <- g_normal(c_or_n$y, theta[8:9])
    g_a1 <- g_normal(c_or_a$y, theta[10:11])
    g_c1 <- g_normal(c_or_a$y, theta[6:7])

    output[dat$z == 0 & dat$d == 0] <- ifelse(stats::rbinom(dim(c_or_n)[1], 1, omega_c*g_c0/(omega_c*g_c0 + omega_n*g_n0)), "c", "n")
    output[dat$z == 1 & dat$d == 1] <- ifelse(stats::rbinom(dim(c_or_a)[1], 1, omega_c*g_c1/(omega_c*g_c1 + omega_a*g_a1)), "c", "a")

  }

  }

  if(strong_access){

    omega_c <- theta[1]
    omega_n <- theta[2]

    c_or_n <- dat[dat$z == 0,]

    output[dat$z == 1 & dat$d == 1] <- "c"
    output[dat$z == 1 & dat$d == 0] <- "n"

    if(!exclusion_restriction){

      g_c0 <- g_normal(c_or_n$y, theta[3:4])
      g_n0 <- g_normal(c_or_n$y, theta[7:8])

      output[dat$z == 0] <- ifelse(stats::rbinom(dim(c_or_n)[1], 1, omega_c*g_c0/(omega_c*g_c0 + omega_n*g_n0)), "c", "n")

      }

    if(exclusion_restriction){

      g_c0 <- g_normal(c_or_n$y, theta[3:4])
      g_n0 <- g_normal(c_or_n$y, theta[7:8])

      output[dat$z == 0] <- ifelse(stats::rbinom(dim(c_or_n)[1], 1, omega_c*g_c0/(omega_c*g_c0 + omega_n*g_n0)), "c", "n")


      }


  }


  output
}

# #' Compliance Parameter Posterior Draw
# #'
# #' \code{draw_compliance} makes a posterior draw from the compliance type parameters given the imputed compliance types
# #'
# #' @param augmented_data a dataframe with columns y (outcome), z (treatment assignment indicator),
# #' d (observed treatment), and c (imputed compliance type).
# #' @param hyper_parameters a numeric vector containing the parameter values for the prior distribution
# #' on the compliance type parameters.
# #' @param strong_access a logical value indicating whether the Strong Access Monotonicity assumption should be made.
# #'
# #' @keywords internal
# #'
# #' @examples
# #' aug_dat <- data.frame(y = rnorm(10), z = rep(c(0, 1), each = 5), d = c(rep(0, 6), rep(1, 4)),
# #' c = c("c", "n", "c", "c", "n", "n", rep("c", 4)))
# #' hyper_params <- rep(1, 5)
# #' draw_compliance(aug_dat, hyper_parameters = hyper_params, strong_access = TRUE)
# #'
# #' @return a numeric vector containing the draw from the posterior distribution of the compliance type parameters
draw_compliance <- function(augmented_data, hyper_parameters, strong_access = T){

  if(!strong_access){

    N_c <- nrow(augmented_data[augmented_data$c == "c", ])
    N_n <- nrow(augmented_data[augmented_data$c == "n", ])
    N_a <- nrow(augmented_data[augmented_data$c == "a", ])

    omega_draw <- MCMCpack::rdirichlet(n = 1, c(N_c + hyper_parameters[1], N_n + hyper_parameters[2], N_a + hyper_parameters[3]))

    }

  if(strong_access){

    N_c <- nrow(augmented_data[augmented_data$c == "c", ])
    N_n <- nrow(augmented_data[augmented_data$c == "n", ])

    omega_draw <- MCMCpack::rdirichlet(n = 1, c(N_c + hyper_parameters[1], N_n + hyper_parameters[2]))


    }
  omega_draw



}

# #' Binary Model Posterior Draw
# #'
# #' \code{draw_binary} makes a posterior draw from the outcome parameters given the imputed compliance types
# #' when the outcome is binary.
# #'
# #' @param augmented_data a dataframe with columns y (outcome), z (treatment assignment indicator),
# #' d (observed treatment), and c (imputed compliance type).
# #' @param hyper_parameters a vector containing the parameter values that determine the prior distribution.
# #' of the model parameters. See Details for more information on how the values of the vector should be ordered.
# #' @param exclusion_restriction a logical value indicating whether the exclusion restriction should be assumed.
# #' @param strong_access a logical value indicating whether the Strong Access Monotonicity restriction should be assumed.
# #' @examples
# #' aug_dat <- data.frame(y = rbinom(10, 1, .5), z = rep(c(0, 1), each = 5), d = c(rep(0, 6), rep(1, 4)),
# #' c = c("c", "n", "c", "c", "n", "n", rep("c", 4)))
# #' hyper_params <- rep(1, 8)
# #' draw_binary(aug_dat, hyper_parameters = hyper_params, exclusion_restriction = TRUE, strong_access = TRUE)
# #' @keywords internal
# #'
# #' @return a vector containing the draw from the posterior.
# #'
draw_binary <- function(augmented_data, hyper_parameters, exclusion_restriction = T, strong_access = T){

  if(!strong_access){
  never_takers <- augmented_data[augmented_data$c == "n", ]
  always_takers <- augmented_data[augmented_data$c == "a", ]
  compliers <- augmented_data[augmented_data$c == "c", ]

  N <- nrow(augmented_data)
  N_c <- nrow(compliers)
  N_n <- nrow(never_takers)
  N_c1 <- nrow(compliers[compliers$z == 1, ])
  N_c0 <- nrow(compliers[compliers$z == 0, ])
  N_n1 <- nrow(never_takers[never_takers$z == 1, ])
  N_n0 <- nrow(never_takers[never_takers$z == 0, ])
  N_y1_c1 <- nrow(compliers[compliers$y == 1 & compliers$z == 1, ])
  N_y1_c0 <- nrow(compliers[compliers$y == 1 & compliers$z == 0, ])
  N_y1_n1 <- nrow(never_takers[never_takers$y == 1 & never_takers$z == 1, ])
  N_y1_n0 <- nrow(never_takers[never_takers$y == 1 & never_takers$z == 0, ])
  N_a <- nrow(always_takers)
  N_a1 <- nrow(always_takers[always_takers$z == 1, ])
  N_a0 <- nrow(always_takers[always_takers$z == 0, ])
  N_y1_a1 <- nrow(always_takers[always_takers$y == 1 & always_takers$z == 1, ])
  N_y1_a0 <- nrow(always_takers[always_takers$y == 1 & always_takers$z == 0 , ])

  p_c0_draw <- stats::rbeta(1, hyper_parameters[4] + N_y1_c0, hyper_parameters[5] + N_c0 - N_y1_c0)
  p_c1_draw <- stats::rbeta(1, hyper_parameters[6] + N_y1_c1, hyper_parameters[7] + N_c1 - N_y1_c1)


  if(!exclusion_restriction){



  p_n0_draw <- stats::rbeta(1, hyper_parameters[8] + N_y1_n0, hyper_parameters[9] + N_n0 - N_y1_n0)
  p_n1_draw <- stats::rbeta(1, hyper_parameters[10] + N_y1_n1, hyper_parameters[11] + N_n1 - N_y1_n1)

  p_a0_draw <- stats::rbeta(1, hyper_parameters[12] + N_y1_a0, hyper_parameters[13] + N_a0 - N_y1_a0)
  p_a1_draw <- stats::rbeta(1, hyper_parameters[14] + N_y1_a1, hyper_parameters[15] + N_a1 - N_y1_a1)

  return(c(p_c0_draw, p_c1_draw, p_n0_draw, p_n1_draw, p_a0_draw, p_a1_draw))


  }

  if(exclusion_restriction){

    p_n_draw <- stats::rbeta(1, hyper_parameters[8] + N_y1_n1 + N_y1_n0, hyper_parameters[9] + N_n - (N_y1_n1 + N_y1_n0))

    p_a_draw <- stats::rbeta(1, hyper_parameters[10] + N_y1_a1 + N_y1_a0, hyper_parameters[11] + N_a - (N_y1_a1 + N_y1_a0))

    return(c(p_c0_draw, p_c1_draw, p_n_draw, p_a_draw))

  }

  }

  if(strong_access){

    never_takers <- augmented_data[augmented_data$c == "n", ]
    compliers <- augmented_data[augmented_data$c == "c", ]

    N <- nrow(augmented_data)
    N_c <- nrow(compliers)
    N_n <- nrow(never_takers)
    N_c1 <- nrow(compliers[compliers$z == 1, ])
    N_c0 <- nrow(compliers[compliers$z == 0, ])
    N_n1 <- nrow(never_takers[never_takers$z == 1, ])
    N_n0 <- nrow(never_takers[never_takers$z == 0, ])
    N_y1_c1 <- nrow(compliers[compliers$y == 1 & compliers$z == 1, ])
    N_y1_c0 <- nrow(compliers[compliers$y == 1 & compliers$z == 0, ])
    N_y1_n1 <- nrow(never_takers[never_takers$y == 1 & never_takers$z == 1, ])
    N_y1_n0 <- nrow(never_takers[never_takers$y == 1 & never_takers$z == 0, ])

    p_c0_draw <- stats::rbeta(1, hyper_parameters[3] + N_y1_c0, hyper_parameters[4] + N_c0 - N_y1_c0)
    p_c1_draw <- stats::rbeta(1, hyper_parameters[5] + N_y1_c1, hyper_parameters[6] + N_c1 - N_y1_c1)


    if(!exclusion_restriction){



      p_n0_draw <- stats::rbeta(1, hyper_parameters[7] + N_y1_n0, hyper_parameters[8] + N_n0 - N_y1_n0)
      p_n1_draw <- stats::rbeta(1, hyper_parameters[9] + N_y1_n1, hyper_parameters[10] + N_n1 - N_y1_n1)

      return(c(p_c0_draw, p_c1_draw, p_n0_draw, p_n1_draw))


    }

    if(exclusion_restriction){

      p_n_draw <- stats::rbeta(1, hyper_parameters[7] + N_y1_n1 + N_y1_n0, hyper_parameters[8] + N_n - (N_y1_n1 + N_y1_n0))

      return(c(p_c0_draw, p_c1_draw, p_n_draw))

    }






  }

}

# #' Normal Model Posterior Draw
# #'
# #' \code{draw_normal} makes a posterior draw given the imputed compliance types
# #' when the outcome is modeled with the Normal distribution
# #'
# #' @param augmented_data a dataframe with columns y (outcome), z (treatment assignment indicator),
# #' d (observed treatment), and c (imputed compliance type)
# #' @param hyper_parameters a vector containing the parameter values that determine the prior distribution
# #' of the model parameters. See Details for more information on how the values of the vector should be ordered.
# #' @param exclusion_restriction a logical value indicating whether the Exclusion Restriction assumption should be made.
# #' @param strong_access a logical value indicating whether the Strong Access Monotonicity assumption should be made.
# #' @keywords internal
# #'
# #' @examples
# #' aug_dat <- data.frame(y = rnorm(10), z = rep(c(0, 1), each = 5), d = c(rep(0, 6), rep(1, 4)),
# #' c = c("c", "n", "c", "c", "n", "n", rep("c", 4)))
# #' hyper_params <- c(1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1)
# #' draw_normal(aug_dat, hyper_parameters = hyper_params, exclusion_restriction = TRUE, strong_access = TRUE)
# #' @return a vector containing the draw from the posterior
# #'
draw_normal <- function(augmented_data, hyper_parameters, exclusion_restriction = T, strong_access = T){

  if(!strong_access){

  never_takers <- augmented_data[augmented_data$c == "n", ]
  always_takers <- augmented_data[augmented_data$c == "a", ]
  compliers <- augmented_data[augmented_data$c == "c", ]

  if(!exclusion_restriction){
  c0_posterior_params <- posterior_params_ninvgamma(compliers$y[compliers$z ==0 ], hyper_parameters[4:7])
  c1_posterior_params <- posterior_params_ninvgamma(compliers$y[compliers$z == 1], hyper_parameters[8:11])
  n0_posterior_params <- posterior_params_ninvgamma(never_takers$y[never_takers$z == 0], hyper_parameters[12:15])
  n1_posterior_params <- posterior_params_ninvgamma(never_takers$y[never_takers$z == 1], hyper_parameters[16:19])
  a0_posterior_params <- posterior_params_ninvgamma(always_takers$y[always_takers$z == 0], hyper_parameters[20:23])
  a1_posterior_params <- posterior_params_ninvgamma(always_takers$y[always_takers$z == 1], hyper_parameters[24:27])

  s2_c0 <- MCMCpack::rinvgamma(1, c0_posterior_params[3], c0_posterior_params[4])
  s2_c1 <- MCMCpack::rinvgamma(1, c1_posterior_params[3], c1_posterior_params[4])
  s2_n0 <- MCMCpack::rinvgamma(1, n0_posterior_params[3], n0_posterior_params[4])
  s2_n1 <- MCMCpack::rinvgamma(1, n1_posterior_params[3], n1_posterior_params[4])
  s2_a0 <- MCMCpack::rinvgamma(1, a0_posterior_params[3], a0_posterior_params[4])
  s2_a1 <- MCMCpack::rinvgamma(1, a1_posterior_params[3], a1_posterior_params[4])

  mu_c0 <- stats::rnorm(1, mean = c0_posterior_params[1], sd = sqrt(c0_posterior_params[2]*s2_c0))
  mu_c1 <- stats::rnorm(1, mean = c1_posterior_params[1], sd = sqrt(c1_posterior_params[2]*s2_c1))
  mu_n0 <- stats::rnorm(1, mean = n0_posterior_params[1], sd = sqrt(n0_posterior_params[2]*s2_n0))
  mu_n1 <- stats::rnorm(1, mean = n1_posterior_params[1], sd = sqrt(n1_posterior_params[2]*s2_n1))
  mu_a0 <- stats::rnorm(1, mean = a0_posterior_params[1], sd = sqrt(a0_posterior_params[2]*s2_a0))
  mu_a1 <- stats::rnorm(1, mean = a1_posterior_params[1], sd = sqrt(a1_posterior_params[2]*s2_a1))

  return(c(mu_c0, s2_c0, mu_c1, s2_c1, mu_n0, s2_n0, mu_n1, s2_n1, mu_a0, s2_a0, mu_a1, s2_a1))

  }

  if(exclusion_restriction){


    c0_posterior_params <- posterior_params_ninvgamma(compliers$y[compliers$z == 0], hyper_parameters[4:7])
    c1_posterior_params <- posterior_params_ninvgamma(compliers$y[compliers$z == 1], hyper_parameters[8:11])
    n_posterior_params <- posterior_params_ninvgamma(never_takers$y, hyper_parameters[12:15])
    a_posterior_params <- posterior_params_ninvgamma(always_takers$y, hyper_parameters[16:19])

    s2_c0 <- MCMCpack::rinvgamma(1, c0_posterior_params[3], c0_posterior_params[4])
    s2_c1 <- MCMCpack::rinvgamma(1, c1_posterior_params[3], c1_posterior_params[4])
    s2_n <- MCMCpack::rinvgamma(1, n_posterior_params[3], n_posterior_params[4])
    s2_a <- MCMCpack::rinvgamma(1, a_posterior_params[3], a_posterior_params[4])

    mu_c0 <- stats::rnorm(1, mean = c0_posterior_params[1], sd = sqrt(c0_posterior_params[2]*s2_c0))
    mu_c1 <- stats::rnorm(1, mean = c1_posterior_params[1], sd = sqrt(c1_posterior_params[2]*s2_c1))
    mu_n <- stats::rnorm(1, mean = n_posterior_params[1], sd = sqrt(n_posterior_params[2]*s2_n))
    mu_a <- stats::rnorm(1, mean = a_posterior_params[1], sd = sqrt(a_posterior_params[2]*s2_a))

    return(c(mu_c0, s2_c0, mu_c1, s2_c1, mu_n, s2_n, mu_a, s2_a))

  }

  }

  if(strong_access){

    never_takers <- augmented_data[augmented_data$c == "n", ]
    compliers <- augmented_data[augmented_data$c == "c", ]

    if(!exclusion_restriction){
      c0_posterior_params <- posterior_params_ninvgamma(compliers$y[compliers$z == 0], hyper_parameters[3:6])
      c1_posterior_params <- posterior_params_ninvgamma(compliers$y[compliers$z == 1], hyper_parameters[7:10])
      n0_posterior_params <- posterior_params_ninvgamma(never_takers$y[never_takers$z == 0], hyper_parameters[11:14])
      n1_posterior_params <- posterior_params_ninvgamma(never_takers$y[never_takers$z == 1], hyper_parameters[15:18])

      s2_c0 <- MCMCpack::rinvgamma(1, c0_posterior_params[3], c0_posterior_params[4])
      s2_c1 <- MCMCpack::rinvgamma(1, c1_posterior_params[3], c1_posterior_params[4])
      s2_n0 <- MCMCpack::rinvgamma(1, n0_posterior_params[3], n0_posterior_params[4])
      s2_n1 <- MCMCpack::rinvgamma(1, n1_posterior_params[3], n1_posterior_params[4])

      mu_c0 <- stats::rnorm(1, mean = c0_posterior_params[1], sd = sqrt(c0_posterior_params[2]*s2_c0))
      mu_c1 <- stats::rnorm(1, mean = c1_posterior_params[1], sd = sqrt(c1_posterior_params[2]*s2_c1))
      mu_n0 <- stats::rnorm(1, mean = n0_posterior_params[1], sd = sqrt(n0_posterior_params[2]*s2_n0))
      mu_n1 <- stats::rnorm(1, mean = n1_posterior_params[1], sd = sqrt(n1_posterior_params[2]*s2_n1))

      return(c(mu_c0, s2_c0, mu_c1, s2_c1, mu_n0, s2_n0, mu_n1, s2_n1))

    }

    if(exclusion_restriction){


      c0_posterior_params <- posterior_params_ninvgamma(compliers$y[compliers$z == 0], hyper_parameters[3:6])
      c1_posterior_params <- posterior_params_ninvgamma(compliers$y[compliers$z == 1], hyper_parameters[7:10])
      n_posterior_params <- posterior_params_ninvgamma(never_takers$y, hyper_parameters[11:14])

      s2_c0 <- MCMCpack::rinvgamma(1, c0_posterior_params[3], c0_posterior_params[4])
      s2_c1 <- MCMCpack::rinvgamma(1, c1_posterior_params[3], c1_posterior_params[4])
      s2_n <- MCMCpack::rinvgamma(1, n_posterior_params[3], n_posterior_params[4])

      mu_c0 <- stats::rnorm(1, mean = c0_posterior_params[1], sd = sqrt(c0_posterior_params[2]*s2_c0))
      mu_c1 <- stats::rnorm(1, mean = c1_posterior_params[1], sd = sqrt(c1_posterior_params[2]*s2_c1))
      mu_n <- stats::rnorm(1, mean = n_posterior_params[1], sd = sqrt(n_posterior_params[2]*s2_n))

      return(c(mu_c0, s2_c0, mu_c1, s2_c1, mu_n, s2_n))

    }
  }
}

# #' Data Augmentation for Binary Outcome with Non-compliance
# #'
# #' \code{compliance_chain_binary} performs a Bayesian analysis of data from an RCT subject to non-compliance
# #'  with a binary outcome by running a single chain of the data augmentation algorithm.
# #'
# #' @param dat a dataframe with columns y (the outcome), z (the treatment assignment),
# #' and d (the treatment actually received).
# #' @param exclusion_restriction a logical value indicating whether the Exclusion Restriction assumption should be made.
# #' @param strong_access a logical value indicating whether the Strong Access Monotonicity assumption should be made.
# #' @param starting_values the initial parameter values. If NULL, then the initial parameter values are determined
# #' by a random draw from the prior distribution.
# #' @param hyper_parameters a numerical vector containing the values that determine the prior distribution
# #' for the model parameters. If NULL, then the hyper parameters are chosen to give non-informative priors.
# #' @param n_iter number of iterations of the data augmentation algorithm to perform
# #' @param n_burn number of initial iterations to discard
# #'
# #' @examples
# #' dat <- data.frame(y = rep(c(0, 1), 5), z = rep(c(0, 1), each = 5), d = c(rep(0, 6), rep(1, 4)) )
# #' compliance_chain_binary(dat, exclusion_restriction = TRUE, strong_access = TRUE, starting_values = NULL, hyper_parameters = NULL,
# #' n_iter = 10, n_burn = 2)
# #'
# #'
# #' @keywords internal
# #'
# #' @return a matrix containing the posterior draws based on the data augmentation
compliance_chain_binary <- function(dat, exclusion_restriction, strong_access, starting_values, hyper_parameters, n_iter, n_burn){



    if(!strong_access){


      if(is.null(hyper_parameters)){

        if(!exclusion_restriction) hyper_parameters <- rep(1, 15)
        if(exclusion_restriction) hyper_parameters <- rep(1, 11)


      }

      if(is.null(starting_values)){

        if(!exclusion_restriction) starting_values <- c(MCMCpack::rdirichlet(1, hyper_parameters[1:3]), stats::rbeta(6, hyper_parameters[seq(from = 4, to = 14, by = 2)], hyper_parameters[seq(from = 5, to = 15, by = 2)]))
        if(exclusion_restriction) starting_values <- c(MCMCpack::rdirichlet(1, hyper_parameters[1:3]), stats::rbeta(4, hyper_parameters[seq(from = 4, to = 10, by = 2)], hyper_parameters[seq(from = 5, to = 11, by = 2)]))


      }

    }

    if(strong_access){
      if(is.null(hyper_parameters)){

        if(!exclusion_restriction) hyper_parameters <- rep(1, 10)
        if(exclusion_restriction) hyper_parameters <- rep(1, 8)


      }

      if(is.null(starting_values)){

        if(!exclusion_restriction) starting_values <- c(MCMCpack::rdirichlet(1, hyper_parameters[1:2]), stats::rbeta(4, hyper_parameters[seq(from = 3, to = 9, by = 2)], hyper_parameters[seq(from = 4, to = 10, by = 2)]))
        if(exclusion_restriction) starting_values <- c(MCMCpack::rdirichlet(1, hyper_parameters[1:2]), stats::rbeta(3, hyper_parameters[seq(from = 3, to = 7, by = 2)], hyper_parameters[seq(from = 4, to = 8, by = 2)]))


      }



    }


    theta_vals <- matrix(NA, nrow = n_iter + 1, ncol = length(starting_values))
    theta_vals[1, ] <- starting_values

    aug_dat <- dat
    aug_dat$c <- NA





    if(!strong_access){

      for(i in 1:n_iter){


        aug_dat$c <- impute_c_binary(dat = dat, theta = theta_vals[i, ], exclusion_restriction = exclusion_restriction, strong_access = strong_access)
        theta_vals[i + 1, 1:3] <- draw_compliance(augmented_data = aug_dat, hyper_parameters = hyper_parameters[1:3], strong_access = strong_access)

        theta_vals[i+1, 4:ncol(theta_vals) ] <- draw_binary(augmented_data = aug_dat, hyper_parameters = hyper_parameters,
                                                            exclusion_restriction = exclusion_restriction, strong_access = strong_access)






      }

    }

    if(strong_access){

      for(i in 1:n_iter){


        aug_dat$c <- impute_c_binary(dat = dat, theta = theta_vals[i, ], exclusion_restriction = exclusion_restriction)

        theta_vals[i + 1, 1:2] <- draw_compliance(augmented_data = aug_dat, hyper_parameters = hyper_parameters)

        theta_vals[i+1, 3:ncol(theta_vals) ] <- draw_binary(augmented_data = aug_dat, hyper_parameters = hyper_parameters,
                                                            exclusion_restriction = exclusion_restriction)






      }




    }

    theta_vals





}

# #' Data Augmentation for Normal Outcome with Non-compliance
# #'
# #' \code{compliance_chain_normal} performs a Bayesian analysis of data from an RCT subject to non-compliance
# #'  with a Normal outcome by running a single chain of the data augmentation algorithm.
# #'
# #' @param dat a dataframe with columns y (the outcome), z (the treatment assignment),
# #' and d (the treatment actually received).
# #' @param exclusion_restriction a logical value indicating whether the Exclusion Restriction assumption should be made.
# #' @param strong_access a logical value indicating whether the Strong Access Monotonicity assumption should be made.
# #' @param starting_values the initial parameter values. If NULL, then the initial parameter values are determined
# #' by a random draw from the prior distribution.
# #' @param hyper_parameters a numerical vector containing the values that determine the prior distribution
# #' for the model parameters. If NULL, then the hyper parameters are chosen to give non-informative priors.
# #' @param n_iter number of iterations of the data augmentation algorithm to perform
# #' @param n_burn number of initial iterations to discard
# #'
# #' @examples
# #' dat <- data.frame(y = rnorm(10), z = rep(c(0, 1), each = 5), d = c(rep(0, 7), 1, 1, 1))
# #' compliance_chain_normal(dat, exclusion_restriction = TRUE, strong_access = TRUE, starting_values = NULL,
# #' hyper_parameters = NULL, n_iter = 10, n_burn = 2)
# #'
# #' @keywords internal
# #'
# #' @return a matrix containing the posterior draws based on the data augmentation
compliance_chain_normal <- function(dat, exclusion_restriction, strong_access, starting_values, hyper_parameters, n_iter, n_burn){



    reference_prior <- F
    y0_suffstat <- c(mean(dat$y[dat$z == 0]), stats::var(dat$y[dat$z == 0]))
    y1_suffstat <- c(mean(dat$y[dat$z == 1]), stats::var(dat$y[dat$z == 1]))
    if(!strong_access){

      if(is.null(hyper_parameters)){
        reference_prior <- T


        if(!exclusion_restriction) hyper_parameters <- c(1, 1, 1, rep(c(0, 1000, 1, 1), 6) )
        if(exclusion_restriction) hyper_parameters <- c(rep(1, 3), rep(c(0, 1000, 1, 1), 4))

      }

      if(is.null(starting_values)){


        if(!exclusion_restriction){

          starting_values <- c(MCMCpack::rdirichlet(1, hyper_parameters[1:3]), rep(c(y0_suffstat, y1_suffstat), 3))

        }

        if(exclusion_restriction){


          starting_values <- c(MCMCpack::rdirichlet(1, hyper_parameters[1:3]), c(y0_suffstat, rep(y1_suffstat, 3)) )


        }


      }




    }

    if(strong_access){

      if(is.null(hyper_parameters)){

        reference_prior <- T

        if(!exclusion_restriction) hyper_parameters <- c(1, 1, rep(c(0, 10, 1, 1), 4) )
        if(exclusion_restriction) hyper_parameters <- c(rep(1, 2), rep(c(0, 10, 1, 1), 3))

      }

      if(is.null(starting_values)){

        n_hyperparameters <- length(hyper_parameters)
        if(!exclusion_restriction){


          starting_values <- c(MCMCpack::rdirichlet(1, hyper_parameters[1:2]), rep(c(y0_suffstat, y1_suffstat), 2))


        }

        if(exclusion_restriction){


          starting_values <- c(MCMCpack::rdirichlet(1, hyper_parameters[1:2]), c(y0_suffstat, rep(y1_suffstat, 2)))


        }


      }




    }



  theta_vals <- matrix(NA, nrow = n_iter + 1, ncol = length(starting_values))
  theta_vals[1, ] <- starting_values

  aug_dat <- dat
  aug_dat$c <- NA

  if(reference_prior){
    # reset hyperparameters to special case
    if(strong_access){

      if(exclusion_restriction) hyper_parameters <- c(1, 1, rep(c(0, Inf, -1/2, 0), 3))
      if(!exclusion_restriction) hyper_parameters <- c(1, 1, rep(c(0, Inf, -1/2, 0), 4))

      }

    if(!strong_access){

      if(exclusion_restriction) hyper_parameters <- c(1, 1, 1,  rep(c(0, Inf, -1/2, 0), 4))
      if(!exclusion_restriction) hyper_parameters <- c(1, 1, 1,  rep(c(0, Inf, -1/2, 0), 6))


    }

  }


  if(strong_access){



  for(i in 1:n_iter){



    aug_dat$c <- impute_c_normal(dat = dat, theta = theta_vals[i, ],
                                 exclusion_restriction = exclusion_restriction, strong_access = T)

    theta_vals[i + 1, 1:2] <- draw_compliance(augmented_data = aug_dat, hyper_parameters = hyper_parameters, strong_access = strong_access)
    theta_vals[i + 1, 3:ncol(theta_vals)] <- draw_normal(augmented_data = aug_dat, hyper_parameters = hyper_parameters,
                                                         exclusion_restriction = exclusion_restriction, strong_access = strong_access)



          }
  }

  if(!strong_access){


    for(i in 1:n_iter){



      aug_dat$c <- impute_c_normal(dat = dat, theta = theta_vals[i, ],
                                   exclusion_restriction = exclusion_restriction, strong_access = F)

      theta_vals[i + 1, 1:3] <- draw_compliance(augmented_data = aug_dat, hyper_parameters = hyper_parameters, strong_access = strong_access)
      theta_vals[i + 1, 4:ncol(theta_vals)] <- draw_normal(augmented_data = aug_dat, hyper_parameters = hyper_parameters,
                                                           exclusion_restriction = exclusion_restriction, strong_access = strong_access)



    }

  }

  theta_vals

}

#' Data Augmentation for Non-compliance analysis
#'
#' \code{compliance_chain} fits a Bayesian non-compliance model by running a single chain of the data augmentation algorithm
#'
#' @param dat a data frame. The first column of the data frame should be the outcome variable,
#' the second column should be the treatment assignment variable, and the third column should be the
#' indicator for the treatment actually received.
#' @param outcome_model a character string indicating how the outcome should be modeled. Either "normal" for the
#' normal model or "binary" for the binary model.
#' @param exclusion_restriction a logical value indicating whether the exclusion restriction assumption should be made.
#' @param strong_access a logical value indicating whether the strong access monotonicity assumption should be made.
#' @param starting_values the initial parameter values. If NULL, then the initial parameter values are based on either a random draw
#' from the prior distribution (for the binary model) or sample statistics (for the normal model).
#' @param hyper_parameters a numerical vector containing the values that determine the prior distribution
#' for the model parameters. If NULL, then the hyper parameters are chosen to give non-informative or reference priors.
#' @param n_iter number of iterations of the data augmentation algorithm to perform.
#' @param n_burn number of initial iterations to discard.
#'
#' @examples
#' # runs 10 iterations of the data augmentation algorithm on a subset of the vitaminA data
#' set.seed(4923)
#' compliance_chain(vitaminA[sample(1:nrow(vitaminA), 1000),], outcome_model = "binary",
#' exclusion_restriction = TRUE, strong_access = TRUE, n_iter = 10, n_burn = 1)
#'
#' @export
#' @return a matrix containing the draws from the posterior distribution.
compliance_chain <- function(dat, outcome_model = NULL, exclusion_restriction = T, strong_access = T, starting_values = NULL, hyper_parameters = NULL, n_iter = 10000,
                             n_burn = 1000){

  names(dat)[1:3] <- c("y", "z", "d")


  if(is.null(outcome_model)) stop('You must specify an outcome model.')


  if(outcome_model == "binary"){

    theta_vals <- compliance_chain_binary(dat = dat, exclusion_restriction = exclusion_restriction,
                                          strong_access = strong_access, starting_values = starting_values,
                                          hyper_parameters = hyper_parameters, n_iter = n_iter)

  }

  if(outcome_model == "normal"){

    theta_vals <- compliance_chain_normal(dat = dat, exclusion_restriction = exclusion_restriction,
                                          strong_access = strong_access, starting_values = starting_values,
                                          hyper_parameters = hyper_parameters, n_iter = n_iter)

  }




  output <- theta_vals[(n_burn + 1):(n_iter+1), ]

  colnames(output) <- name_columns(outcome_model = outcome_model, exclusion_restriction = exclusion_restriction,
                                   strong_access = strong_access)

  output

}

#' Compute the Posterior Distribution of the CACE
#'
#' \code{cace} takes a sample from the posterior distribution of the model parameters
#' and computes the corresponding posterior distribution of the Complier Average Causal Effect.
#'
#' @param chain a matrix containing the draws from the posterior distribution of the model parameters.
#' The matrix should either be the result of a call to \code{compliance_chain} or have the same
#' structure as one.
#' @param outcome_model a character string indicating which outcome model was used in fitting the model,
#' either "binary" for a dichotomous outcome or "normal" for the Normal model.
#' @param strong_access a logical indicating whether the strong access monotonicity assumption was
#' made when fitting the model
#'
#' @examples
#' # CACE based on a subset of the vitaminA dataset
#' set.seed(4923)
#' chain <- compliance_chain(vitaminA[sample(1:nrow(vitaminA), 1000),], outcome_model = "binary",
#'  exclusion_restriction = TRUE, strong_access = TRUE, n_iter = 10, n_burn = 1)
#'
#' cace(chain, outcome_model = "binary", strong_access = TRUE)
#'
#' # matrix representing the samples from the posterior distribution of the model parameters
#' posterior_mat <- matrix(rnorm(10*8, mean = 10), nrow = 10, ncol = 8)
#' cace(posterior_mat, "normal", strong_access = TRUE)
#'
#' @return a vector containing the draws from the posterior distribution of the CACE
#'
#' @export
cace <- function(chain, outcome_model, strong_access){

  if(outcome_model == "binary"){

    if(strong_access){

      output <- chain[,4] - chain[,3]

    }

    if(!strong_access){

      output <- chain[,5] - chain[,4]

    }


  }

  if(outcome_model == "normal"){

    if(strong_access){

      output <- chain[,5] - chain[,3]


    }

    if(!strong_access){

      output <- chain[,6] - chain[,4]

    }




  }

output

}


#' Posterior Inference based on a Sample from the Posterior
#'
#' \code{summarize_chain} provides posterior summaries based off a sample from the
#' posterior distribution.
#'
#' @param chain a numeric vector containing the samples from the posterior distribution.
#' @param digits the number of decimal places
#'
#' @examples
#' # Suppose the posterior distribution was Normal(15, 5)
#' posterior_chain <- rnorm(100, 15, 5); summarize_chain(posterior_chain)
#'
#' @return a list containing the posterior mean, median, and quantile-based credible intervals
#' calculated from the values in the chain.
#'
#' @export
summarize_chain <- function(chain, digits = 3){

  if(class(chain) == "matrix") stop("You supplied a matrix. Please supply a vector instead!")

  chain_mean <- round(mean(chain), digits = 3)
  chain_median <- round(stats::median(chain), digits = 3)
  chain_50ci <- round(c(stats::quantile(x = chain, .25), stats::quantile(x = chain, .75)), digits = 3)
  chain_90ci <- round(c(stats::quantile(x = chain, .05), stats::quantile(x = chain, .95)), digits = 3)
  chain_95ci <- round(c(stats::quantile(x = chain, .025), stats::quantile(x = chain, .975)), digits = 3)

  cat("Posterior Mean:", chain_mean, "\n")
  cat("Posterior Median:", chain_median, "\n")
  cat("Posterior 50% Credible Interval: (",chain_50ci[1], ", ",chain_50ci[2], ") \n", sep = "")
  cat("Posterior 90% Credible Interval: (", chain_90ci[1], ", ", chain_90ci[2], ") \n", sep = "")
  cat("Posterior 95% Credible Interval: (", chain_95ci[1], ", ", chain_95ci[2], ") \n", sep = "")

  invisible(list(mean = chain_mean, median = chain_median, ci_50 = chain_50ci, ci_90 = chain_90ci, ci_95 = chain_95ci))









}



