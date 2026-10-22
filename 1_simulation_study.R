
## This script was refined with assistance from
## OpenAI ChatGPT (GPT-5.6 Sol, September 2026).
## The final implementation and outputs were verified by the author.

## 1 Source functions
source('0_functions.R')

## 2 Management protocols

set.seed(1)

sim_early <- simulate_managed_hawkes(
  T_end = 500,
  gamma = 0,
  protocol = "early",
  K = 3,
  window = 4,
  early_delay = 0.1
)

sim_quiet <- simulate_managed_hawkes(
  T_end = 500,
  gamma = 0,
  protocol = "quiet",
  K = 3,
  window = 4,
  quiet_period = 3
)

c(
  events_early = length(sim_early$event_times),
  management_early = length(sim_early$management_times),
  
  events_quiet = length(sim_quiet$event_times),
  management_quiet = length(sim_quiet$management_times)
)

## Plot

plot_managed_process <- function(
    sim,
    xlim = c(0, 100),
    ylab1 ='',
    ylab2='',
    title=''
) {
  
  events <- sim$event_times[
    sim$event_times >= xlim[1] &
      sim$event_times <= xlim[2]
  ]
  
  management <- sim$management_times[
    sim$management_times >= xlim[1] &
      sim$management_times <= xlim[2]
  ]
  
  plot(
    NA,
    xlim = xlim,
    ylim = c(0.5, 2.5),
    yaxt = "n",
    xlab = "Time",
    ylab = "",
    main = title,
  )
  
  segments(
    events,
    0.8,
    events,
    1.2
  )
  
  segments(
    management,
    1.8,
    management,
    2.2,
    lwd = 2
  )
  
  axis(
    2,
    at = c(1, 2),
    labels = c(
      ylab1,
      ylab2
    ),
    las = 1
  )
}

## not used in main manuscript anymore, but in the supplement
pdf(file = 'Plots/plot_sim.pdf', width = 12, height = 4)
par(mfrow = c(1, 2))
# Plot 1: Left = 8, Right = 1 (Total = 9)
par(mar = c(5, 8, 1.5, 1))
plot_managed_process(
  sim_early, ylab1 = "Ecological event", ylab2 = "Management",
  title ='Early',
  xlim = c(0, 50)
)

# Plot 2: Left = 1, Right = 8 (Total = 9)
par(mar = c(5, 1, 1.5, 8))
plot_managed_process(
  sim_quiet, ylab1 = "", ylab2 = "",
  xlim = c(0, 50),
  title = 'Quiet'
)
dev.off()


# define the 6 scenarios

scenario_grid <- expand.grid(
  gamma_true = c(
    -0.20,
    0,
    0.20
  ),
  protocol = c(
    "early",
    "quiet"
  ),
  stringsAsFactors = FALSE
)

scenario_grid$scenario <- with(
  scenario_grid,
  paste0(
    ifelse(
      gamma_true < 0,
      "Negative true effect",
      ifelse(
        gamma_true > 0,
        "Positive true effect",
        "No true effect"
      )
    ),
    " / ",
    protocol
  )
)

scenario_grid

run_simulation_scenarios <- function(
    n_sim = 300,
    seed = 12345,
    ...
) {
  
  set.seed(seed)
  
  scenarios <- expand.grid(
    gamma_true = c(
      -0.20,
      0,
      0.20
    ),
    
    protocol = c(
      "early",
      "quiet"
    ),
    
    stringsAsFactors = FALSE
  )
  
  output <- list()
  counter <- 1L
  
  for (s in seq_len(nrow(scenarios))) {
    
    gamma_true <-
      scenarios$gamma_true[s]
    
    protocol <-
      scenarios$protocol[s]
    
    message(
      "Running: gamma = ",
      gamma_true,
      ", protocol = ",
      protocol
    )
    
    for (b in seq_len(n_sim)) {
      
      result <- tryCatch(
        
        run_one_replicate(
          gamma_true = gamma_true,
          protocol = protocol,
          ...
        ),
        
        error = function(e) {
          NULL
        }
      )
      
      if (!is.null(result)) {
        
        result$replicate <- b
        
        output[[counter]] <- result
        counter <- counter + 1L
      }
    }
  }
  
  do.call(
    rbind,
    output
  )
}


# run

results <- run_simulation_scenarios(
  n_sim = 300,
  
  T_end = 500,
  
  mu = 0.20,
  alpha = 0.60,
  beta_N = 1,
  
  beta_R = 1,
  
  K = 3,
  window = 4,
  
  early_delay = 0.10,
  
  # Important parameter for downward bias
  quiet_period = 3,
  
  cooldown = 3,
  
  grid_size = 4000
)


results_ok <- subset(
  results,
  convergence_mis == 0 &
    convergence_full == 0
)


## bias calculation

results_ok$bias_mis <-
  results_ok$gamma_hat_mis -
  results_ok$gamma_true

results_ok$bias_full <-
  results_ok$gamma_hat_full -
  results_ok$gamma_true


summary_results <- do.call(
  rbind,
  
  lapply(
    split(
      results_ok,
      list(
        results_ok$gamma_true,
        results_ok$protocol
      )
    ),
    
    function(x) {
      
      data.frame(
        gamma_true =
          unique(x$gamma_true),
        
        protocol =
          unique(x$protocol),
        
        n =
          nrow(x),
        
        mean_gamma_mis =
          mean(x$gamma_hat_mis),
        
        mean_gamma_full =
          mean(x$gamma_hat_full),
        
        bias_mis =
          mean(x$bias_mis),
        
        bias_full =
          mean(x$bias_full),
        
        mean_cov_HR_HN =
          mean(
            x$history_covariance
          ),
        
        mean_cor_HR_HN =
          mean(
            x$history_correlation
          )
      )
    }
  )
)

rownames(summary_results) <- NULL

summary_results


# boxplot

results_ok$group <- paste(
  paste0(
    "gamma0 = ",
    results_ok$gamma_true
  ),
  results_ok$protocol,
  sep = "\n"
)

# Create expression vector with formatted math and line breaks
x_labs <- c(
  expression(atop(gamma[0] == -0.2, "early")),
  expression(atop(gamma[0] == -0.2, "quiet")),
  expression(atop(gamma[0] == 0, "early")),
  expression(atop(gamma[0] == 0, "quiet")),
  expression(atop(gamma[0] == 0.2, "early")),
  expression(atop(gamma[0] == 0.2, "quiet"))
)

# Pass directly to the boxplot

pdf(file='Plots/Figure_2_a.pdf', width=6, height=4)
par(mgp = c(3, 2, 0))
boxplot(  bias_mis ~ group,
  data = results_ok,
  pch=19, cex=0.2,
  ylab =
    expression(
      hat(gamma)[mis] - gamma[0]
    ),
  names= x_labs,
  xlab='',
  main =
    "Bias from ignoring self-excitation"
)

abline(
  h = 0,
  lty = 2
)
dev.off()


# compare with correctly specified model

pdf(file='Plots/Figure_2_b.pdf', width=6, height=4)
par(mgp = c(3, 2, 0))
boxplot(
  bias_full ~ group,
  data = results_ok,
  pch=19, cex=0.2,
  ylab =
    expression(
      hat(gamma)[full] - gamma[0]
    ),
  names= x_labs,
  xlab='',
  main =
    "Bias after accounting for self-excitation"
)

abline(
  h = 0,
  lty = 2
)
dev.off()


## useful figure

means <- aggregate(
  cbind(
    gamma_hat_mis,
    gamma_hat_full
  ) ~ gamma_true + protocol,
  
  data = results_ok,
  
  FUN = mean
)

means


# history covariance plot - not used in the end in manuscript

# pdf(file='Plots/covariance.pdf', width=13, height=6)
# plot(
#   results_ok$history_covariance,
#   results_ok$bias_mis,
#   pch=19, cex=0.5,
#   xlab =
#     expression(
#       Cov(H[R](t), H[N](t))
#     ),
#   
#   ylab =
#     expression(
#       hat(gamma)[mis] - gamma[0]
#     )
# )
# 
# abline(
#   h = 0,
#   lty = 2
# )
# 
# abline(
#   v = 0,
#   lty = 2
# )
# dev.off()



## Now the other figures. From timing experiment

# ============================================================
# TIMING EXPERIMENT: EARLY/REACTIVE MANAGEMENT
#
# Vary intervention delay relative to Hawkes timescale
#
# d = early_delay * beta_N
#
# ============================================================



timing_results <- run_early_timing_experiment(
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
  
  gamma_values = c(
    -0.2,
    0,
    0.2
  ),
  
  n_sim = 300,
  
  T_end = 500,
  
  mu = 0.20,
  alpha = 0.60,
  beta_N = 1,
  
  beta_R = 1,
  
  K = 3,
  window = 4,
  
  cooldown = 3,
  
  grid_size = 4000
)


timing_results <- subset(
  timing_results,
  convergence_mis == 0 &
    convergence_full == 0
)

timing_results$bias_mis <-
  timing_results$gamma_hat_mis -
  timing_results$gamma_true

timing_results$bias_full <-
  timing_results$gamma_hat_full -
  timing_results$gamma_true


# Make plotting factors with controlled ordering
timing_results$delay_factor <- factor(
  timing_results$delay_relative,
  levels = c(0.05, 0.10, 0.25, 0.50, 1, 2, 3, 5)
)

timing_results$gamma_factor <- factor(
  timing_results$gamma_true,
  levels = c(-0.2, 0, 0.2),
  labels = c(
    expression(gamma[0] == -0.2),
    expression(gamma[0] == 0),
    expression(gamma[0] == 0.2)
  )
)

###############
#### FIGURE 3
###############

delay_levels <- c(
  0.05,
  0.10,
  0.25,
  0.50,
  1,
  2,
  3,
  5
)

gamma_levels <- c(
  -0.2,
  0,
  0.2
)


gamma_fills <- c(
  "white",
  "grey75",
  "grey40"
)

box_cols <- rep(
  gamma_fills,
  times = length(delay_levels)
)


# List of vectors in desired order:
# for each delay, gamma = -0.2, 0, 0.2
box_data_mis <- list()

for (d in delay_levels) {
  
  for (g in gamma_levels) {
    
    box_data_mis[[length(box_data_mis) + 1]] <-
      timing_results$bias_mis[
        timing_results$delay_relative == d &
          timing_results$gamma_true == g
      ]
  }
}

n_gamma <- length(gamma_levels)

group_centres <- seq(
  from = 2,
  by = 5,
  length.out = length(delay_levels)
)

box_positions <- unlist(
  lapply(
    group_centres,
    function(x) x + c(-1, 0, 1)
  )
)


pdf(file='Plots/Figure_3.pdf', width=10, height=5)
boxplot(
  box_data_mis,
  at = box_positions,
  col = box_cols,
  outline = FALSE,
  xaxt = "n",
  xlab = expression(D * beta[N]),
  ylab = expression(hat(gamma)[mis] - gamma[0]),
  main = "Bias and intervention timing when ignoring self-excitation"
)

axis(
  1,
  at = group_centres,
  labels = delay_levels
)

abline(
  h = 0,
  lty = 2
)

legend(
  "topright",
  legend = c(
    expression(gamma[0] == -0.2),
    expression(gamma[0] == 0),
    expression(gamma[0] == 0.2)
  ),
  fill = gamma_fills,
  bty = "n"
)
dev.off()

###############
#### FIGURE 4
###############

plot3_data <- subset(
  timing_results,
  gamma_true == 0 &
    is.finite(history_covariance) &
    is.finite(bias_mis)
)

delay_levels <- sort(
  unique(plot3_data$delay_relative)
)

plot3_summary <- do.call(
  rbind,
  lapply(
    delay_levels,
    function(d) {
      
      z <- subset(
        plot3_data,
        delay_relative == d
      )
      
      data.frame(
        delay_relative = d,
        
        cov_median =
          median(z$history_covariance/z$history_variance),
        
        cov_q25 =
          quantile(z$history_covariance/z$history_variance, 0.25),
        
        cov_q75 =
          quantile(z$history_covariance/z$history_variance, 0.75),
        
        bias_median =
          median(z$bias_mis),
        
        bias_q25 =
          quantile(z$bias_mis, 0.25),
        
        bias_q75 =
          quantile(z$bias_mis, 0.75)
      )
    }
  )
)



pdf(file='Plots/Figure_4.pdf', width=10, height=5)
plot(
  plot3_data$history_covariance/plot3_data$history_variance,
  plot3_data$bias_mis,
  
  pch = 16,
  cex = 0.45,
  col = adjustcolor("grey40", alpha.f = 0.15),
  
  xlab =
    expression(
      Cov[T](H[R](t), H[N](t))/ Var[T](H[R](t))
    ),
  
  ylab =
    expression(
      hat(gamma)[mis] - gamma[0]
    )
  # ,
  # 
  # main =
  #   expression(gamma[0] == 0)
)

# abline(
#   h = 0,
#   lty = 2
# )
# 
# abline(
#   v = 0,
#   lty = 2
# )

abline(a=0, b=0.6, lty=2)

segments(
  x0 = plot3_summary$cov_q25,
  y0 = plot3_summary$bias_median,
  
  x1 = plot3_summary$cov_q75,
  y1 = plot3_summary$bias_median,
  
  lwd = 2
)

segments(
  x0 = plot3_summary$cov_median,
  y0 = plot3_summary$bias_q25,
  
  x1 = plot3_summary$cov_median,
  y1 = plot3_summary$bias_q75,
  
  lwd = 2
)
points(
  plot3_summary$cov_median,
  plot3_summary$bias_median,
  
  pch = 21,
  bg = "white",
  cex = 1,
  lwd = 1.5
)
# lines(
#   plot3_summary$cov_median,
#   plot3_summary$bias_median,
#   lwd = 1.5
# )

# text(
#   plot3_summary$cov_median,
#   plot3_summary$bias_median,
#   
#   labels =
#     plot3_summary$delay_relative,
#   
#   pos = 4,
#   cex = 0.5
# )
dev.off()


###############
#### FIGURE 5
###############


### BRANCHING RATIO?

run_delay_excitation_grid <- function(
    delay_relative = c(
      0.05,
      0.1,
      0.25,
      0.5,
      1,
      2,
      3,
      5
    ),
    
    branching_ratios = c(
      0,
      0.2,
      0.4,
      0.6,
      0.8
    ),
    
    gamma_true = 0,
    
    n_sim = 300,
    
    beta_N = 1,
    
    ...
) {
  
  output <- list()
  index <- 1L
  
  for (eta in branching_ratios) {
    
    alpha <-
      eta * beta_N
    
    for (d_relative in delay_relative) {
      
      early_delay <-
        d_relative / beta_N
      
      message(
        "eta = ", eta,
        "; D beta = ", d_relative
      )
      
      for (b in seq_len(n_sim)) {
        
        result <- tryCatch(
          
          run_one_replicate(
            gamma_true =
              gamma_true,
            
            protocol =
              "early",
            
            alpha =
              alpha,
            
            beta_N =
              beta_N,
            
            early_delay =
              early_delay,
            
            ...
          ),
          
          error =
            function(e) NULL
        )
        
        if (!is.null(result)) {
          
          result$branching_ratio <-
            eta
          
          result$delay_relative <-
            d_relative
          
          result$replicate <-
            b
          
          output[[index]] <-
            result
          
          index <- index + 1L
        }
      }
    }
  }
  
  result <- do.call(
    rbind,
    output
  )
  
  result$bias_mis <-
    result$gamma_hat_mis -
    result$gamma_true
  
  result$bias_full <-
    result$gamma_hat_full -
    result$gamma_true
  
  result
}

grid_results <- run_delay_excitation_grid(
  delay_relative = c(
    0.05,
    0.1,
    0.25,
    0.5,
    1,
    2,
    3,
    5
  ),
  
  branching_ratios = c(
    0,
    0.2,
    0.4,
    0.6,
    0.8
  ),
  
  gamma_true = 0,
  
  n_sim = 300,
  
  beta_N = 1,
  
  T_end = 500,
  
  mu = 0.20,
  
  beta_R = 1,
  
  K = 3,
  window = 4,
  
  quiet_period = 3,
  cooldown = 3,
  
  grid_size = 4000
)


grid_summary <- do.call(
  rbind,
  
  lapply(
    split(
      grid_results,
      list(
        grid_results$branching_ratio,
        grid_results$delay_relative
      )
    ),
    
    function(x) {
      
      data.frame(
        branching_ratio =
          unique(x$branching_ratio),
        
        delay_relative =
          unique(x$delay_relative),
        
        mean_bias =
          mean(
            x$bias_mis,
            na.rm = TRUE
          ),
        
        median_bias =
          median(
            x$bias_mis,
            na.rm = TRUE
          ),
        
        mean_covariance =
          mean(
            x$history_covariance,
            na.rm = TRUE
          ),
        
        mean_n_management =
          mean(
            x$n_management,
            na.rm = TRUE
          )
      )
    }
  )
)

rownames(grid_summary) <- NULL


etas <-
  sort(
    unique(
      grid_summary$branching_ratio
    )
  )



pdf(file='Plots/Figure_5.pdf', width=10, height=5)
plot(
  NA,
  
  xlim =
    range(
      grid_summary$delay_relative
    ),
  
  ylim =
    range(
      grid_summary$mean_bias
    ),
  
  log = "x",
  
  xlab =
    expression(D * beta[N]),
  
  ylab =
    expression(
      E(
        hat(gamma)[mis] -
          gamma[0]
      )
    )
  # ,
  # 
  # main =
  #   expression(
  #     "True " * gamma[0] == 0
  #   )
)

for (i in seq_along(etas)) {
  
  eta <- etas[i]
  
  temp <- subset(
    grid_summary,
    branching_ratio == eta
  )
  
  lines(
    temp$delay_relative,
    temp$mean_bias,
    
    type = "b",
    
    lty = i,
    
    pch = i
  )
}

abline(
  h = 0,
  lty = 3
)

legend(
  "topright",
  
  legend =
    paste0(
      "eta = ",
      etas
    ),
  
  lty =
    seq_along(etas),
  
  pch =
    seq_along(etas),
  
  bty = "n"
)
dev.off()
