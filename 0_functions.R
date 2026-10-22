# ============================================================
# MANAGED SELF-EXCITING ECOLOGICAL PROCESS
#
# Ecological events N:
#
# lambda_N(t) =
#   [mu + alpha * H_N(t) + gamma * H_R(t)]_+
#
# H_N(t) = sum exp[-beta_N (t - t_i)]
# H_R(t) = sum exp[-beta_R (t - r_j)]
#
#
# Management protocol:
#
# "early":
#   Once K events occur within a window of length W,
#   management is scheduled after a short delay D.
#
# "quiet":
#   Once K events occur within W, management becomes "armed".
#   It occurs only after Q event-free time units.
#   Every new ecological event resets the Q-clock.
#
# ============================================================

simulate_managed_hawkes <- function(
    T_end = 500,
    
    # Ecological process
    mu = 0.20,
    alpha = 0.60,
    beta_N = 1.00,
    
    # True management effect
    gamma = 0,
    
    # Duration of management effect
    beta_R = 1.00,
    
    # Management triggering
    protocol = c("early", "quiet"),
    K = 3,
    window = 4,
    
    # Early protocol
    early_delay = 0.10,
    
    # Quiet protocol
    quiet_period = 3,
    
    # Prevent immediate repeated management
    cooldown = 3
) {
  
  protocol <- match.arg(protocol)
  
  if (alpha / beta_N >= 1) {
    stop("alpha / beta_N must be < 1.")
  }
  
  # Current time
  t <- 0
  
  # Hawkes state variables
  H_N <- 0
  H_R <- 0
  
  # Event records
  event_times <- numeric(0)
  management_times <- numeric(0)
  
  # Management state
  management_pending <- FALSE
  scheduled_management <- Inf
  
  last_management <- -Inf
  cooldown_until <- -Inf
  
  while (t < T_end) {
    
    # --------------------------------------------------------
    # Upper bound for ecological intensity
    #
    # Negative gamma can only REDUCE intensity,
    # so ignoring it provides a valid upper bound.
    # --------------------------------------------------------
    
    lambda_upper <-
      mu +
      alpha * H_N +
      max(gamma, 0) * H_R
    
    lambda_upper <- max(lambda_upper, 1e-12)
    
    # Candidate ecological event
    waiting_time <- rexp(1, rate = lambda_upper)
    candidate_time <- t + waiting_time
    
    # --------------------------------------------------------
    # Does scheduled management occur first?
    # --------------------------------------------------------
    
    if (
      management_pending &&
      scheduled_management < candidate_time &&
      scheduled_management <= T_end
    ) {
      
      dt <- scheduled_management - t
      
      # Decay both histories
      H_N <- H_N * exp(-beta_N * dt)
      H_R <- H_R * exp(-beta_R * dt)
      
      t <- scheduled_management
      
      # Management occurs
      management_times <- c(
        management_times,
        t
      )
      
      # H_R jumps by 1
      H_R <- H_R + 1
      
      # Reset management state
      management_pending <- FALSE
      scheduled_management <- Inf
      
      last_management <- t
      cooldown_until <- t + cooldown
      
      next
    }
    
    # No more ecological events before T
    if (candidate_time > T_end) {
      break
    }
    
    # --------------------------------------------------------
    # Move to ecological candidate time
    # --------------------------------------------------------
    
    dt <- candidate_time - t
    
    H_N_candidate <-
      H_N * exp(-beta_N * dt)
    
    H_R_candidate <-
      H_R * exp(-beta_R * dt)
    
    # Actual ecological intensity
    lambda_candidate <- max(
      0,
      mu +
        alpha * H_N_candidate +
        gamma * H_R_candidate
    )
    
    # Ogata thinning
    accepted <-
      runif(1) <
      lambda_candidate / lambda_upper
    
    # Advance time regardless
    t <- candidate_time
    H_N <- H_N_candidate
    H_R <- H_R_candidate
    
    if (!accepted) {
      next
    }
    
    # ========================================================
    # Ecological event occurs
    # ========================================================
    
    event_times <- c(
      event_times,
      t
    )
    
    # Self-excitation jump
    H_N <- H_N + 1
    
    # --------------------------------------------------------
    # Management logic
    # --------------------------------------------------------
    
    if (t < cooldown_until) {
      next
    }
    
    # Count events since both:
    # 1. start of recent window
    # 2. previous management
    lower_limit <- max(
      t - window,
      last_management
    )
    
    recent_events <-
      event_times[
        event_times > lower_limit &
          event_times <= t
      ]
    
    threshold_reached <-
      length(recent_events) >= K
    
    # ========================================================
    # EARLY MANAGEMENT
    # ========================================================
    
    if (protocol == "early") {
      
      if (
        threshold_reached &&
        !management_pending
      ) {
        
        management_pending <- TRUE
        
        scheduled_management <-
          t + early_delay
      }
    }
    
    # ========================================================
    # QUIET-PHASE MANAGEMENT
    #
    # Once threshold reached, intervention is planned.
    # Every new event resets the quiet-period clock.
    # ========================================================
    
    # if (protocol == "quiet") {
    #   
    #   if (threshold_reached) {
    #     
    #     management_pending <- TRUE
    #     
    #     # Reset after EVERY ecological event
    #     scheduled_management <-
    #       t + quiet_period
    #   }
    # }
    
    # changed: once armed, management occurs only after \(Q\) consecutive event-free time units.
    if (protocol == "quiet") {
      
      if (management_pending) {
        scheduled_management <- t + quiet_period
        
      } else if (threshold_reached) {
        management_pending <- TRUE
        scheduled_management <- t + quiet_period
      }
    }
  }
  
  list(
    event_times = event_times,
    management_times = management_times,
    T_end = T_end,
    
    parameters = list(
      mu = mu,
      alpha = alpha,
      beta_N = beta_N,
      gamma = gamma,
      beta_R = beta_R,
      protocol = protocol,
      K = K,
      window = window,
      early_delay = early_delay,
      quiet_period = quiet_period
    )
  )
}



exponential_history <- function(
    target_times,
    source_times,
    beta
) {
  
  target_times <- sort(target_times)
  source_times <- sort(source_times)
  
  if (length(target_times) == 0) {
    return(numeric(0))
  }
  
  H <- numeric(length(target_times))
  
  source_index <- 1L
  
  for (i in seq_along(target_times)) {
    
    t <- target_times[i]
    
    if (source_index == 1L) {
      
      relevant <-
        source_times[source_times < t]
      
    } else {
      
      relevant <-
        source_times[source_times < t]
    }
    
    if (length(relevant) > 0) {
      
      H[i] <-
        sum(
          exp(
            -beta *
              (t - relevant)
          )
        )
    }
  }
  
  H
}


prepare_fit_data <- function(
    event_times,
    management_times,
    T_end,
    beta_N = 1,
    beta_R = 1,
    grid_size = 4000
) {
  
  # Histories at event times
  HN_events <- exponential_history(
    event_times,
    event_times,
    beta_N
  )
  
  HR_events <- exponential_history(
    event_times,
    management_times,
    beta_R
  )
  
  # Midpoint integration grid
  dt <- T_end / grid_size
  
  grid <- (
    seq_len(grid_size) - 0.5
  ) * dt
  
  HN_grid <- exponential_history(
    grid,
    event_times,
    beta_N
  )
  
  HR_grid <- exponential_history(
    grid,
    management_times,
    beta_R
  )
  
  list(
    event_times = event_times,
    management_times = management_times,
    T_end = T_end,
    
    HN_events = HN_events,
    HR_events = HR_events,
    
    HN_grid = HN_grid,
    HR_grid = HR_grid,
    
    grid = grid,
    dt = dt,
    
    beta_N = beta_N,
    beta_R = beta_R
  )
}

## likelihoods

hawkes_signed_nll <- function(
    theta,
    fit_data,
    include_self = TRUE
) {
  
  # Positive background rate
  mu <- exp(theta[1])
  
  if (include_self) {
    
    # Transform alpha so branching ratio < 1
    branching_ratio <-
      0.995 * plogis(theta[2])
    
    alpha <-
      branching_ratio *
      fit_data$beta_N
    
    gamma <- theta[3]
    
  } else {
    
    alpha <- 0
    gamma <- theta[2]
  }
  
  # ----------------------------------
  # Intensity at observed event times
  # ----------------------------------
  
  eta_events <-
    mu +
    alpha * fit_data$HN_events +
    gamma * fit_data$HR_events
  
  # Every observed event must have
  # strictly positive fitted intensity
  if (
    any(eta_events <= 0) ||
    any(!is.finite(eta_events))
  ) {
    return(1e100)
  }
  
  # ----------------------------------
  # Numerical compensator
  # ----------------------------------
  
  eta_grid <-
    mu +
    alpha * fit_data$HN_grid +
    gamma * fit_data$HR_grid
  
  lambda_grid <-
    pmax(
      eta_grid,
      0
    )
  
  compensator <-
    sum(lambda_grid) *
    fit_data$dt
  
  loglik <-
    sum(log(eta_events)) -
    compensator
  
  -loglik
}


## fit correct of misspecified model

fit_signed_hawkes <- function(
    fit_data,
    include_self = TRUE
) {
  
  empirical_rate <-
    length(fit_data$event_times) /
    fit_data$T_end
  
  if (include_self) {
    
    start <- c(
      log(max(empirical_rate / 2, 0.01)),
      
      # branching ratio initial value ~ 0.4
      qlogis(0.4 / 0.995),
      
      # gamma
      0
    )
    
  } else {
    
    start <- c(
      log(max(empirical_rate, 0.01)),
      0
    )
  }
  
  fit <- optim(
    par = start,
    fn = hawkes_signed_nll,
    fit_data = fit_data,
    include_self = include_self,
    method = "Nelder-Mead",
    control = list(
      maxit = 5000,
      reltol = 1e-9
    )
  )
  
  if (include_self) {
    
    mu_hat <- exp(fit$par[1])
    
    branching_hat <-
      0.995 * plogis(fit$par[2])
    
    alpha_hat <-
      branching_hat *
      fit_data$beta_N
    
    gamma_hat <- fit$par[3]
    
    coefficients <- c(
      mu = mu_hat,
      alpha = alpha_hat,
      branching_ratio =
        branching_hat,
      gamma = gamma_hat
    )
    
  } else {
    
    coefficients <- c(
      mu = exp(fit$par[1]),
      gamma = fit$par[2]
    )
  }
  
  list(
    coefficients = coefficients,
    logLik = -fit$value,
    convergence = fit$convergence,
    optim = fit
  )
}

## one replicate test

run_one_replicate <- function(
    gamma_true,
    protocol,
    
    T_end = 500,
    
    mu = 0.20,
    alpha = 0.60,
    beta_N = 1,
    
    beta_R = 1,
    
    K = 3,
    window = 4,
    
    early_delay = 0.1,
    quiet_period = 3,
    
    cooldown = 3,
    
    grid_size = 4000
) {
  
  sim <- simulate_managed_hawkes(
    T_end = T_end,
    
    mu = mu,
    alpha = alpha,
    beta_N = beta_N,
    
    gamma = gamma_true,
    beta_R = beta_R,
    
    protocol = protocol,
    
    K = K,
    window = window,
    
    early_delay = early_delay,
    quiet_period = quiet_period,
    
    cooldown = cooldown
  )
  
  # Need at least a few interventions
  if (
    length(sim$event_times) < 20 ||
    length(sim$management_times) < 3
  ) {
    return(NULL)
  }
  
  fit_data <- prepare_fit_data(
    event_times =
      sim$event_times,
    
    management_times =
      sim$management_times,
    
    T_end = T_end,
    
    beta_N = beta_N,
    beta_R = beta_R,
    
    grid_size = grid_size
  )
  
  # Misspecified: omit self-excitation
  fit_mis <- fit_signed_hawkes(
    fit_data,
    include_self = FALSE
  )
  
  # Correct Hawkes model
  fit_full <- fit_signed_hawkes(
    fit_data,
    include_self = TRUE
  )
  
  # ------------------------------------------
  # Diagnostic quantity:
  #
  # Cov{H_R(t), H_N(t)}
  #
  # This should help explain sign of bias.
  # ------------------------------------------
  
  history_covariance <-
    cov(
      fit_data$HR_grid,
      fit_data$HN_grid
    )
  
  # add this now: variance of H_R
  history_variance <-
    var(
      fit_data$HR_grid
    )
  
  history_correlation <-
    cor(
      fit_data$HR_grid,
      fit_data$HN_grid
    )
  
  data.frame(
    gamma_true = gamma_true,
    protocol = protocol,
    
    n_events =
      length(sim$event_times),
    
    n_management =
      length(sim$management_times),
    
    gamma_hat_mis =
      fit_mis$coefficients["gamma"],
    
    gamma_hat_full =
      fit_full$coefficients["gamma"],
    
    alpha_hat_full =
      fit_full$coefficients["alpha"],
    
    history_covariance =
      history_covariance,
    
    history_variance =
      history_variance,
    
    history_correlation =
      history_correlation,
    
    convergence_mis =
      fit_mis$convergence,
    
    convergence_full =
      fit_full$convergence
  )
}





# ============================================================
# TIMING EXPERIMENT: EARLY/REACTIVE MANAGEMENT
#
# Vary intervention delay relative to Hawkes timescale
#
# d = early_delay * beta_N
#
# ============================================================

run_early_timing_experiment <- function(
    delay_relative = c(
      0.05,
      0.10,
      0.25,
      0.50,
      1,
      2,
      3,
      5
    ),
    gamma_values = c(-0.2, 0, 0.2),
    n_sim = 300,
    seed = 12345,
    
    T_end = 500,
    
    mu = 0.20,
    alpha = 0.60,
    beta_N = 1,
    
    beta_R = 1,
    
    K = 3,
    window = 4,
    cooldown = 3,
    
    grid_size = 4000
) {
  
  set.seed(seed)
  
  output <- list()
  index <- 1L
  
  for (gamma_true in gamma_values) {
    
    for (d_relative in delay_relative) {
      
      # Convert dimensionless delay D * beta_N
      # back to actual time units
      early_delay <-
        d_relative / beta_N
      
      message(
        "gamma = ", gamma_true,
        "; D * beta_N = ", d_relative
      )
      
      for (b in seq_len(n_sim)) {
        
        result <- tryCatch(
          
          run_one_replicate(
            gamma_true = gamma_true,
            protocol = "early",
            
            T_end = T_end,
            
            mu = mu,
            alpha = alpha,
            beta_N = beta_N,
            
            beta_R = beta_R,
            
            K = K,
            window = window,
            
            early_delay = early_delay,
            
            # Irrelevant under early protocol
            quiet_period = 3,
            
            cooldown = cooldown,
            
            grid_size = grid_size
          ),
          
          error = function(e) NULL
        )
        
        if (!is.null(result)) {
          
          result$delay_relative <-
            d_relative
          
          result$early_delay <-
            early_delay
          
          result$replicate <- b
          
          output[[index]] <- result
          index <- index + 1L
        }
      }
    }
  }
  
  if (length(output) == 0) {
    stop("No successful simulations.")
  }
  
  results <- do.call(
    rbind,
    output
  )
  
  rownames(results) <- NULL
  
  results
}