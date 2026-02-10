# Entropy experiment on a discrete state space (1..13)
# install.packages("rootSolve")
library(rootSolve)

states <- 1:13

beta_from_mean <- function(target_mean, interval = c(-4, 4)) {
  if (target_mean <= min(states) || target_mean >= max(states)) {
    stop("target_mean must be strictly between 1 and 13.")
  }

  mean_residual <- function(beta) {
    w <- exp(-beta * states)
    z <- sum(w)
    expected_value <- sum(states * w) / z
    expected_value - target_mean
  }

  uniroot(mean_residual, interval)$root
}

entropy_from_beta <- function(beta) {
  w <- exp(-beta * states)
  p <- w / sum(w)
  # Numerical safety for p*log(p) at extreme beta values.
  p <- p[p > 0]
  -sum(p * log(p))
}

probability_from_beta <- function(beta) {
  w <- exp(-beta * states)
  w / sum(w)
}

# Example outputs
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
beta_vals <- sapply(E_vals, beta_from_mean)
temp_vals <- 1 / beta_vals

beta_grid <- seq(-4, 4, by = 0.01)
entropy_beta <- sapply(beta_grid, entropy_from_beta)

T_grid <- seq(-20, 20, by = 0.1)
finite_T <- is.finite(T_grid) & abs(T_grid) > 1e-8
T_grid <- T_grid[finite_T]
entropy_T <- sapply(1 / T_grid, entropy_from_beta)

entropy_E <- sapply(E_vals, function(e) entropy_from_beta(beta_from_mean(e)))

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
     xlab = "Fixed Mean E", ylab = "Value",
     main = "Beta(E) and Temperature(E)")
lines(E_vals, temp_vals, lwd = 2.2, col = "#d62728")
legend("topright",
       legend = c(expression(beta(E)), "T(E) = 1/beta"),
       col = c("#1f77b4", "#d62728"),
       lwd = 2.2,
       bty = "n",
       cex = 0.9)
grid(col = "#d1d5db", lty = 3)

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
