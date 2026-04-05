min_max_std <- function(a, b) {
  x <- range(a, b)
  std <- diff(x) / 2
  mean <- mean(x)
  c(mean, std)
}

mc_stats <- function(x, ...){
  median_s  <- median(x)
  ci_95    <- quantile(x, c(0.025, 0.975))
  ci_68    <- quantile(x, c(0.16, 0.84))

  log_s <- log10(x)
  mean_log_s <- mean(log_s)
  s_mean <- 10^mean_log_s

  sd_log_s   <- sd(log_s)
  se_log_s   <- sd_log_s / sqrt(length(x))
  sde_range <- 10^c(mean_log_s - se_log_s, mean_log_s + se_log_s)

  list(
    median = median_s |> units::set_units(...),
    mean = s_mean |> units::set_units(...),
    sde = sde_range |> units::set_units(...),
    ir_95 = ci_95 |> units::set_units(...),
    ir_68 = ci_68 |> units::set_units(...),
    samples = x |> units::set_units(...)
  )
}

