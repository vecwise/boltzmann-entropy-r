# Entropy experiment on a discrete state space (1..13)
# The script uses reduced units with k_B = 1, so beta = 1 / T.

states <- 1:13
tolerance_beta_zero <- 1e-10

stable_weights <- function(beta) {
  log_w <- -beta * states
  shift <- max(log_w)
  weights <- exp(log_w - shift)

  list(
    weights = weights,
    log_partition = shift + log(sum(weights))
  )
}

probability_from_beta <- function(beta) {
  stable <- stable_weights(beta)
  stable$weights / sum(stable$weights)
}

mean_from_beta <- function(beta) {
  probabilities <- probability_from_beta(beta)
  sum(states * probabilities)
}

find_bracketing_interval <- function(target_mean, initial_bound = 1, max_bound = 1024) {
  lower <- -initial_bound
  upper <- initial_bound
  lower_value <- mean_from_beta(lower) - target_mean
  upper_value <- mean_from_beta(upper) - target_mean

  while (lower_value * upper_value > 0 && upper < max_bound) {
    lower <- lower * 2
    upper <- upper * 2
    lower_value <- mean_from_beta(lower) - target_mean
    upper_value <- mean_from_beta(upper) - target_mean
  }

  if (lower_value * upper_value > 0) {
    stop("Failed to bracket a beta value for the requested mean.")
  }

  c(lower, upper)
}

beta_from_mean <- function(target_mean) {
  min_state <- min(states)
  max_state <- max(states)
  uniform_mean <- mean(states)

  if (target_mean <= min_state || target_mean >= max_state) {
    stop(sprintf("target_mean must be strictly between %g and %g.", min_state, max_state))
  }

  if (isTRUE(all.equal(target_mean, uniform_mean))) {
    return(0)
  }

  mean_residual <- function(beta) {
    mean_from_beta(beta) - target_mean
  }

  interval <- find_bracketing_interval(target_mean)
  uniroot(mean_residual, interval = interval, tol = 1e-12)$root
}

entropy_from_beta <- function(beta) {
  probabilities <- probability_from_beta(beta)
  positive_probabilities <- probabilities[probabilities > 0]
  -sum(positive_probabilities * log(positive_probabilities))
}

entropy_from_mean <- function(target_mean) {
  entropy_from_beta(beta_from_mean(target_mean))
}

temperature_from_beta <- function(beta, tolerance = tolerance_beta_zero) {
  ifelse(abs(beta) < tolerance, NA_real_, 1 / beta)
}

validate_model <- function() {
  sample_means <- c(1.01, 1.1, 7, 12.9, 12.99)
  recovered_betas <- vapply(sample_means, beta_from_mean, numeric(1))
  recovered_means <- vapply(recovered_betas, mean_from_beta, numeric(1))
  entropy_identity <- vapply(
    recovered_betas,
    function(beta) {
      stable <- stable_weights(beta)
      mean_from_beta(beta) * beta + stable$log_partition
    },
    numeric(1)
  )
  direct_entropy <- vapply(recovered_betas, entropy_from_beta, numeric(1))

  stopifnot(max(abs(recovered_means - sample_means)) < 1e-9)
  stopifnot(max(abs(entropy_identity - direct_entropy)) < 1e-10)
  stopifnot(abs(sum(probability_from_beta(0.5)) - 1) < 1e-12)
}

# Example outputs
validate_model()

beta_for_mean_8 <- beta_from_mean(8)
entropy_at_beta_05 <- entropy_from_beta(0.5)
entropy_for_mean_8 <- entropy_from_beta(beta_for_mean_8)

cat("=== Summary Table: Example Results ===\n")
summary_table <- data.frame(
  Metric = c("beta_from_mean(8)", "entropy_from_beta(0.5)", "entropy_from_beta(beta_from_mean(8))"),
  Value = c(beta_for_mean_8, entropy_at_beta_05, entropy_for_mean_8)
)
print(summary_table, row.names = FALSE)

# Build data once, then plot in a compact 2x2 panel.
E_vals <- seq(1.1, 12.9, by = 0.1)
beta_vals <- vapply(E_vals, beta_from_mean, numeric(1))
temp_vals <- temperature_from_beta(beta_vals)

beta_grid <- seq(-4, 4, by = 0.01)
entropy_beta <- vapply(beta_grid, entropy_from_beta, numeric(1))

T_grid <- seq(-20, 20, by = 0.1)
finite_T <- is.finite(T_grid) & abs(T_grid) > 1e-8
T_grid <- T_grid[finite_T]
entropy_T <- vapply(1 / T_grid, entropy_from_beta, numeric(1))

entropy_E <- vapply(E_vals, entropy_from_mean, numeric(1))

# Optional: inspect one probability distribution.
example_prob <- probability_from_beta(0.5)
cat("\n=== Table: Probability Distribution at beta = 0.5 ===\n")
prob_table <- data.frame(State = states, Probability = round(example_prob, 6))
print(prob_table, row.names = FALSE)

# Plot styling
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par), add = TRUE)
par(mfrow = c(2, 2),
    mar = c(4.2, 4.5, 3.2, 1.2),
    oma = c(0, 0, 2, 0),
    col.axis = "#1f2937",
    col.lab = "#111827")

# Panel 1: Beta(E) and Temperature(E)
plot(E_vals, beta_vals,
     type = "l", lwd = 2.2, col = "#1f77b4",
     xlab = "Fixed Mean E", ylab = expression(beta(E)),
     main = "Beta(E) and Temperature(E)")
grid(col = "#d1d5db", lty = 3)
par(new = TRUE)
plot(E_vals, temp_vals,
     type = "l", lwd = 2.2, col = "#d62728",
     axes = FALSE, xlab = "", ylab = "")
axis(side = 4, col.axis = "#d62728", col = "#d62728")
mtext("Temperature T = 1 / beta", side = 4, line = 2.6, col = "#d62728")
legend("topright",
       legend = c(expression(beta(E)), "T(E) = 1/beta"),
       col = c("#1f77b4", "#d62728"),
       lwd = 2.2,
       bty = "n",
       cex = 0.9)

# Panel 2: Entropy vs Beta
plot(beta_grid, entropy_beta,
     type = "l", lwd = 2.2, col = "#2ca02c",
     xlab = expression(beta), ylab = "Entropy S",
     main = "Entropy vs Beta")
grid(col = "#d1d5db", lty = 3)

# Panel 3: Entropy vs Temperature (finite only)
plot(T_grid, entropy_T,
     type = "l", lwd = 2.2, col = "#9467bd",
     xlab = "Temperature T", ylab = "Entropy S",
     main = "Entropy vs Temperature (Finite T)")
grid(col = "#d1d5db", lty = 3)

# Panel 4: Entropy vs Mean E
plot(E_vals, entropy_E,
     type = "l", lwd = 2.2, col = "#ff7f0e",
     xlab = "Fixed Mean E", ylab = "Entropy S",
     main = "Entropy vs Mean")
grid(col = "#d1d5db", lty = 3)

mtext("Maximum-Entropy / Canonical-Form Relationships", outer = TRUE,
      cex = 1.2, col = "#111827")
