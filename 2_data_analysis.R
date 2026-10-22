

## This script was refined with assistance from
## OpenAI ChatGPT (GPT-5.6 Sol, September 2026).
## The final implementation and outputs were verified by the author.



library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

load(file='Data/wolves_processed.RData')

str(attacks)
str(removals)

range(attacks$date)
range(removals$date)

nrow(attacks)
nrow(removals)

## check


summary(attacks$x)
summary(attacks$y)

summary(removals$x)
summary(removals$y)





## ============================================================
## Matched pseudo-removal analysis
## ============================================================

set.seed(123)

radius <- 10
lookback <- 60
B <- 500

## Function returning attack counts in 5-day bins
get_history <- function(x, y, date, attacks,
                        radius = 10,
                        lookback = 60) {
  
  dist <- sqrt(
    (attacks$x - x)^2 +
      (attacks$y - y)^2
  )
  
  lag <- as.numeric(date - attacks$date)
  
  sapply(seq(0, lookback - 5, by = 5), function(a) {
    
    sum(
      !is.na(dist) &
        !is.na(lag) &
        dist <= radius &
        lag > a &
        lag <= a + 5
    )
  })
}

obs_histories <- t(
  sapply(seq_len(nrow(removals)), function(j) {
    
    get_history(
      removals$x[j],
      removals$y[j],
      removals$date[j],
      attacks,
      radius = radius,
      lookback = lookback
    )
  })
)

dim(obs_histories)


# new control?

management_period <- function(y) {
  ifelse(
    y <= 2014, "2011-2014",
    ifelse(y <= 2018, "2015-2018", "2019-2020")
  )
}

make_same_calendar_date <- function(y, mm, dd) {
  
  candidate <- lubridate::ymd(
    sprintf("%04d-%02d-%02d", y, mm, dd),
    quiet = TRUE
  )
  
  ## Handle 29 February
  if (is.na(candidate) && mm == 2 && dd == 29) {
    candidate <- lubridate::ymd(
      sprintf("%04d-02-28", y)
    )
  }
  
  candidate
}


sample_control_date_same_period <- function(
    d,
    attack_dates,
    removal_dates,
    lookback = 60) {
  
  yy <- year(d)
  mm <- month(d)
  dd <- day(d)
  
  ## Use actual study/removal years, rather than attack years:
  ## avoids accidentally treating a boundary date in 2021 as a full study year.
  study_years <- sort(unique(year(removal_dates)))
  
  candidate_years <- study_years[
    study_years != yy &
      management_period(study_years) == management_period(yy)
  ]
  
  candidate_dates <- as.Date(
    sapply(
      candidate_years,
      function(y)
        make_same_calendar_date(y, mm, dd)
    ),
    origin = "1970-01-01"
  )
  
  ## Require a complete preceding lookback window
  min_allowed <- min(attack_dates, na.rm = TRUE) + days(lookback)
  max_allowed <- max(attack_dates, na.rm = TRUE)
  
  candidate_dates <- candidate_dates[
    !is.na(candidate_dates) &
      candidate_dates >= min_allowed &
      candidate_dates <= max_allowed
  ]
  
  if (length(candidate_dates) == 0) {
    return(as.Date(NA))
  }
  
  candidate_dates[
    sample.int(length(candidate_dates), 1)
  ]
}


## Quick check that every observed removal has at least one eligible control year

n_candidates_same_period <- sapply(
  removals$date,
  function(d) {
    
    yy <- year(d)
    mm <- month(d)
    dd <- day(d)
    
    study_years <- sort(unique(year(removals$date)))
    
    candidate_years <- study_years[
      study_years != yy &
        management_period(study_years) ==
        management_period(yy)
    ]
    
    candidate_dates <- as.Date(
      sapply(
        candidate_years,
        function(y)
          make_same_calendar_date(y, mm, dd)
      ),
      origin = "1970-01-01"
    )
    
    min_allowed <- min(attacks$date) + days(lookback)
    max_allowed <- max(attacks$date)
    
    sum(
      !is.na(candidate_dates) &
        candidate_dates >= min_allowed &
        candidate_dates <= max_allowed
    )
  }
)

table(n_candidates_same_period)
summary(n_candidates_same_period)

## Every observed removal has at least four eligible comparison years, and 239/255 have five. 
## So the corrected matching scheme is well supported

null_curves <- matrix(
  NA,
  nrow = B,
  ncol = 12
)

for (b in seq_len(B)) {
  
  pseudo_histories <- t(
    sapply(seq_len(nrow(removals)), function(j) {
      
      
      pseudo_date <-
        sample_control_date_same_period(
          removals$date[j],
          attacks$date,
          removals$date,
          lookback = lookback
        )
      
      get_history(
        removals$x[j],
        removals$y[j],
        pseudo_date,
        attacks,
        radius = radius,
        lookback = lookback
      )
    })
  )
  
  null_curves[b, ] <-
    colMeans(pseudo_histories)
  print(b)
}


observed_curve <-
  colMeans(obs_histories)

null_mean <-
  colMeans(null_curves)

null_lower <-
  apply(null_curves, 2, quantile, 0.025)

null_upper <-
  apply(null_curves, 2, quantile, 0.975)

plot_data <- data.frame(
  days_before = seq(2.5, 57.5, by = 5),
  observed = observed_curve,
  null = null_mean,
  lower = null_lower,
  upper = null_upper
)

pdf(file='Plots/Figure_6_b.pdf', width=8, height=4)
a <- ggplot(plot_data,
            aes(x = -days_before)) +
  
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.2
  ) +
  
  geom_line(
    aes(y = null),
    linetype = 2
  ) +
  
  geom_line(
    aes(y = observed),
    linewidth = 1
  ) +
  
  geom_point(
    aes(y = observed)
  ) +
  
  labs(
    x = "Days relative to removal",
    y = "Mean attacks per 5-day interval",
    title = "Local attack activity preceding removals",
    subtitle =
      "Observed versus matched pseudo-removals"
    #location- and season-matched pseudo-removals from other years.
  ) +
  
  theme_bw()
a
dev.off()




get_last_attack_delay <- function(x, y, date,
                                  attacks,
                                  radius = 10,
                                  maxlag = 60) {
  
  dist <- sqrt(
    (attacks$x - x)^2 +
      (attacks$y - y)^2
  )
  
  lag <- as.numeric(
    date - attacks$date
  )
  
  valid <-
    !is.na(dist) &
    !is.na(lag) &
    dist <= radius &
    lag >= 1 &
    lag <= maxlag
  
  if (!any(valid))
    return(NA_real_)
  
  min(lag[valid])
}


Q_obs <- sapply(
  seq_len(nrow(removals)),
  function(j)
    get_last_attack_delay(
      removals$x[j],
      removals$y[j],
      removals$date[j],
      attacks
    )
)

# check
c(
  n_removals = length(Q_obs),
  no_attack_60d = sum(is.na(Q_obs)),
  prop_with_attack_60d = mean(!is.na(Q_obs))
)


Q_null <- matrix(
  NA,
  nrow = B,
  ncol = nrow(removals)
)

for (b in seq_len(B)) {
  
  for (j in seq_len(nrow(removals))) {
    
    pseudo_date <-
      sample_control_date_same_period(
        removals$date[j],
        attacks$date,
        removals$date,
        lookback = 60
      )
    
    Q_null[b, j] <-
      get_last_attack_delay(
        removals$x[j],
        removals$y[j],
        pseudo_date,
        attacks
      )
  }
  print(b)
}

## changed here now. 
### For \(q\le60\), a removal with no attack in the previous 60 days should simply contribute FALSE to 
### “occurred within \(q\) days”, rather than disappear from the denominator

T3_obs <- mean(!is.na(Q_obs) & Q_obs <= 3)
T7_obs <- mean(!is.na(Q_obs) & Q_obs <= 7)

T3_null <- apply(Q_null, 1, function(x)
  mean(!is.na(x) & x <= 3))

T7_null <- apply(Q_null, 1, function(x)
  mean(!is.na(x) & x <= 7))

## within 7 days
# T7_obs <- mean(Q_obs <= 7, na.rm = TRUE)
# 
# T7_null <- apply(
#   Q_null,
#   1,
#   function(x)
#     mean(x <= 7, na.rm = TRUE)
# )

mean(T7_null >= T7_obs)


## Fraction of removals within 3 days of an attack
# T3_obs <- mean(Q_obs <= 3, na.rm = TRUE)
# 
# T3_null <- apply(
#   Q_null,
#   1,
#   function(x)
#     mean(x <= 3, na.rm = TRUE)
# )

mean(T3_null >= T3_obs)

data.frame(
  observed = c(T3_obs, T7_obs),
  null = c(mean(T3_null), mean(T7_null)),
  window = c("≤ 3 days", "≤ 7 days")
)

p3 <- (1 + sum(T3_null >= T3_obs)) /
  (length(T3_null) + 1)

p7 <- (1 + sum(T7_null >= T7_obs)) /
  (length(T7_null) + 1)

c(p3 = p3, p7 = p7)

rbind(
  three_days = c(
    observed = T3_obs,
    null_mean = mean(T3_null),
    quantile(T3_null, c(.025, .975))
  ),
  seven_days = c(
    observed = T7_obs,
    null_mean = mean(T7_null),
    quantile(T7_null, c(.025, .975))
  )
)




## escalation

get_escalation <- function(x, y, date,
                           attacks,
                           radius = 10) {
  
  dist <- sqrt(
    (attacks$x - x)^2 +
      (attacks$y - y)^2
  )
  
  lag <- as.numeric(date - attacks$date)
  
  recent <- sum(
    !is.na(dist) &
      !is.na(lag) &
      dist <= radius &
      lag >= 1 &
      lag <= 15
  ) / 15
  
  baseline <- sum(
    !is.na(dist) &
      !is.na(lag) &
      dist <= radius &
      lag >= 30 &
      lag <= 60
  ) / 31
  
  recent - baseline
}

D_obs <- sapply(
  seq_len(nrow(removals)),
  function(j)
    get_escalation(
      removals$x[j],
      removals$y[j],
      removals$date[j],
      attacks
    )
)

mean(D_obs)





library(dplyr)
library(ggplot2)
library(tidyr)

## ---------------------------------------------------------
## Distribution of time since most recent attack:
## observed vs matched pseudo-removals
## ---------------------------------------------------------

q_grid <- 1:30   # show first 30 days; this is where the action is

## Observed empirical CDF
# ecdf_obs <- sapply(q_grid, function(q) {
#   mean(Q_obs <= q, na.rm = TRUE)
# })
ecdf_obs <- sapply(q_grid, function(q)
  mean(!is.na(Q_obs) & Q_obs <= q))

## ECDF for every Monte Carlo pseudo-removal replicate
# ecdf_null <- t(
#   apply(Q_null, 1, function(x) {
#     sapply(q_grid, function(q) {
#       mean(x <= q, na.rm = TRUE)
#     })
#   })
# )

ecdf_null <- t(apply(Q_null, 1, function(x)
  sapply(q_grid, function(q)
    mean(!is.na(x) & x <= q))))

## Null mean + pointwise 95% simulation envelope
plot_Q <- data.frame(
  q = q_grid,
  observed = ecdf_obs,
  null_mean = colMeans(ecdf_null, na.rm = TRUE),
  null_lower = apply(
    ecdf_null, 2, quantile,
    probs = 0.025, na.rm = TRUE
  ),
  null_upper = apply(
    ecdf_null, 2, quantile,
    probs = 0.975, na.rm = TRUE
  )
)


b <- ggplot(plot_Q, aes(x = q)) +
  
  ## Monte Carlo envelope
  geom_ribbon(
    aes(
      ymin = null_lower,
      ymax = null_upper
    ),
    fill = "grey80"
  ) +
  
  ## Mean pseudo-removal distribution
  geom_step(
    aes(y = null_mean),
    linewidth = 0.8,
    # linetype = "dashed",
    colour = "grey35"
  ) +
  
  ## Actual removals
  geom_step(
    aes(y = observed),
    linewidth = 1.1,
    colour = "black"
  ) +
  
  scale_x_continuous(
    breaks = seq(0, 30, by = 5),
    limits = c(0, 30)
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1),
    labels = scales::label_percent(accuracy = 1)
  ) +
  
  geom_vline(
    xintercept = c(3, 7),
    linetype = "dotted",
    linewidth = 0.4,
    colour = "grey60"
  ) +
  labs(
    x = "Days since most recent attack within 10 km",
    y = "Cumulative proportion",
    title = "Removals occur unusually soon after local attacks",
    subtitle = 
      "Actual removals versus matched pseudo-removals ",
    #"(dashed; 95% simulation envelope)"
    
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(), #(face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(size = 11)
  )

pdf(file='Plots/Figure_6_c.pdf', width=7, height=4)
b
dev.off()

Q_null_pool <- as.vector(Q_null)

Q_plot <- bind_rows(
  data.frame(
    Q = Q_obs,
    type = "Actual removals"
  ),
  data.frame(
    Q = Q_null_pool,
    type = "Pseudo-removals"
  )
) %>%
  filter(!is.na(Q), Q <= 30)

ggplot(Q_plot, aes(x = Q, after_stat(density), fill = type)) +
  geom_histogram(
    binwidth = 1,
    boundary = 0,
    position = "identity",
    alpha = 0.4
  ) +
  labs(
    x = "Days since most recent attack within 10 km",
    y = "Proportion",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "top"
  )







## Describe local temporal structure following attacks

radius <- 10
max_lag <- 60

n <- nrow(attacks)

lag_counts <- matrix(
  0,
  nrow = n,
  ncol = max_lag
)

for (i in seq_len(n)) {
  
  dist <- sqrt(
    (attacks$x - attacks$x[i])^2 +
      (attacks$y - attacks$y[i])^2
  )
  
  lag <- as.numeric(
    attacks$date - attacks$date[i]
  )
  
  keep <-
    !is.na(dist) &
    !is.na(lag) &
    dist <= radius &
    lag >= 1 &
    lag <= max_lag
  
  if (any(keep)) {
    lag_counts[i, ] <-
      tabulate(
        lag[keep],
        nbins = max_lag
      )
  }
}

mean_after_attack <-
  colMeans(lag_counts)

plot_cluster <- data.frame(
  lag = 1:max_lag,
  mean_attacks = mean_after_attack
)

c <- ggplot(
  plot_cluster,
  aes(lag, mean_attacks)
) +
  geom_line() +
  labs(
    x = "Days since attack",
    y = "Mean number of subsequent attacks within 10 km",
    title = "Local attack activity following attacks"
  ) +
  theme_classic()

pdf(file='Plots/Figure_6_a.pdf', width=8, height=4)
c
dev.off()

