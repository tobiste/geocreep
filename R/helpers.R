min_max_std <- function(a, b) {
  x <- range(a, b)
  std <- diff(x) / 2
  mean <- mean(x)
  c(mean, std)
}

mc_stats <- function(x, ...){
  median_s  <- stats::median(x)
  #ci_95    <- stats::quantile(x, c(0.025, 0.975))
  #ci_68    <- stats::quantile(x, c(0.16, 0.84))

  log_s <- log10(x)

  log_s_tt <- t.test(log_s)

  mean_log_s <- log_s_tt$estimate
  CI95_log_s <- log_s_tt$conf.int
  stderr_log <- log_s_tt$stderr

  s_mean <- 10^mean_log_s

  #sd_log_s   <- sd(log_s)
  #se_log_s   <- sd_log_s / sqrt(length(x))
  #sde_range <- 10^c(mean_log_s - se_log_s, mean_log_s + se_log_s)
  conf.int <- 10^CI95_log_s

  list(
    median = median_s |> units::set_units(...),
    mean = s_mean |> units::set_units(...),
    # sde = sde_range |> units::set_units(...),
    # ir_95 = ci_95 |> units::set_units(...),
    # ir_68 = ci_68 |> units::set_units(...),
    conf.int = conf.int |> units::set_units(...),
    stderr_log = stderr_log,
    samples = x |> units::set_units(...)
  )
}

