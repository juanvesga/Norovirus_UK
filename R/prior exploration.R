library(ggplot2)
library(gridExtra)
library(dplyr)

# Set up visualization parameters
n_samples <- 10000
x_points <- 1000

# Create a list to store all plots
plot_list <- list()

# 1. factor_beta_1, 2, 4 (all same: Uniform(0.15, 1))
x <- seq(0, 1.2, length.out = x_points)
for (i in c(1, 2, 4)) {
  density <- dunif(x, 0.15, 1)
  plot_list[[paste0("factor_beta_", i)]] <- 
    ggplot(data.frame(x = x, y = density), aes(x, y)) +
    geom_line(linewidth = 1, color = "steelblue") +
    geom_vline(xintercept = c(0.15, 1), linetype = "dashed", color = "red") +
    labs(title = paste0("factor_beta_", i, " ~ Uniform(0.15, 1)"),
         x = paste0("factor_beta_", i), y = "Density") +
    theme_minimal() +
    xlim(0, 1.2)
}

# 2. beta_3 (log_beta_3 transformed)
# log_beta_3 ~ Normal(log(0.18), 0.5)
# beta_3 = exp(log_beta_3) ~ LogNormal(log(0.18), 0.5)
x <- seq(0, 1, length.out = x_points)
density <- dlnorm(x, meanlog = log(0.18), sdlog = 0.5)
median_val <- exp(log(0.18))
q025 <- qlnorm(0.025, meanlog = log(0.18), sdlog = 0.5)
q975 <- qlnorm(0.975, meanlog = log(0.18), sdlog = 0.5)

plot_list[["beta_3"]] <- 
  ggplot(data.frame(x = x, y = density), aes(x, y)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_vline(xintercept = median_val, linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(q025, q975), linetype = "dotted", color = "gray50") +
  labs(title = "beta_3 ~ LogNormal(log(0.18), 0.5)",
       subtitle = sprintf("Median = %.3f, 95%% CI = [%.3f, %.3f]", 
                          median_val, q025, q975),
       x = "beta_3", y = "Density") +
  theme_minimal()

# 3. maternalAB ~ Uniform(1, 365)
x <- seq(0, 400, length.out = x_points)
density <- dunif(x, 1, 365)
plot_list[["maternalAB"]] <- 
  ggplot(data.frame(x = x, y = density), aes(x, y)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_vline(xintercept = c(1, 365), linetype = "dashed", color = "red") +
  labs(title = "maternalAB ~ Uniform(1, 365)",
       subtitle = "Mean = 183 days",
       x = "maternalAB (days)", y = "Density") +
  theme_minimal() +
  xlim(0, 400)

# 4. aduRR ~ Beta(2, 2)
x <- seq(0, 1, length.out = x_points)
density <- dbeta(x, 2, 2)
mean_val <- 2/(2+2)
q025 <- qbeta(0.025, 2, 2)
q975 <- qbeta(0.975, 2, 2)

plot_list[["aduRR"]] <- 
  ggplot(data.frame(x = x, y = density), aes(x, y)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_vline(xintercept = mean_val, linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(q025, q975), linetype = "dotted", color = "gray50") +
  labs(title = "aduRR ~ Beta(2, 2)",
       subtitle = sprintf("Mean = %.3f, 95%% CI = [%.3f, %.3f]", 
                          mean_val, q025, q975),
       x = "aduRR", y = "Density") +
  theme_minimal()

# 5. imm_yr ~ Gamma(shape = 4, scale = 10/4)
x <- seq(0, 15, length.out = x_points)
density <- dgamma(x, shape = 4, scale = 10/4)
mean_val <- 4 * (10/4)
q025 <- qgamma(0.025, shape = 4, scale = 10/4)
q975 <- qgamma(0.975, shape = 4, scale = 10/4)

plot_list[["imm_yr"]] <- 
  ggplot(data.frame(x = x, y = density), aes(x, y)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_vline(xintercept = mean_val, linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(q025, q975), linetype = "dotted", color = "gray50") +
  labs(title = "imm_yr ~ Gamma(4, 2.5)",
       subtitle = sprintf("Mean = %.2f, 95%% CI = [%.2f, %.2f]", 
                          mean_val, q025, q975),
       x = "imm_yr (years)", y = "Density") +
  theme_minimal()

# 6. imm_fac ~ Uniform(1, 5)
x <- seq(0, 6, length.out = x_points)
density <- dunif(x, 1, 5)
plot_list[["imm_fac"]] <- 
  ggplot(data.frame(x = x, y = density), aes(x, y)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_vline(xintercept = c(1, 5), linetype = "dashed", color = "red") +
  labs(title = "imm_fac ~ Uniform(1, 5)",
       subtitle = "Mean = 3",
       x = "imm_fac", y = "Density") +
  theme_minimal() +
  xlim(0, 6)

# 7-9. Reporting ratios (log_ratio_* transformed)
ratio_params <- list(
  list(name = "ratio_0", meanlog = log(2.7), sdlog = 0.5, max_x = 15),
  list(name = "ratio_5", meanlog = log(28), sdlog = 0.5, max_x = 150),
  list(name = "ratio_15", meanlog = log(2.9), sdlog = 0.5, max_x = 15)
)

for (param in ratio_params) {
  x <- seq(0, param$max_x, length.out = x_points)
  density <- dlnorm(x, meanlog = param$meanlog, sdlog = param$sdlog)
  median_val <- exp(param$meanlog)
  q025 <- qlnorm(0.025, meanlog = param$meanlog, sdlog = param$sdlog)
  q975 <- qlnorm(0.975, meanlog = param$meanlog, sdlog = param$sdlog)
  
  plot_list[[param$name]] <- 
    ggplot(data.frame(x = x, y = density), aes(x, y)) +
    geom_line(linewidth = 1, color = "steelblue") +
    geom_vline(xintercept = median_val, linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(q025, q975), linetype = "dotted", color = "gray50") +
    labs(title = paste0(param$name, " ~ LogNormal"),
         subtitle = sprintf("Median = %.2f, 95%% CI = [%.2f, %.2f]", 
                            median_val, q025, q975),
         x = param$name, y = "Density") +
    theme_minimal()
}

# 10. repfac_65p (log_repfac_65p transformed)
x <- seq(0, 300, length.out = x_points)
density <- dlnorm(x, meanlog = log(90), sdlog = 0.3)
median_val <- exp(log(90))
q025 <- qlnorm(0.025, meanlog = log(90), sdlog = 0.3)
q975 <- qlnorm(0.975, meanlog = log(90), sdlog = 0.3)

plot_list[["repfac_65p"]] <- 
  ggplot(data.frame(x = x, y = density), aes(x, y)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_vline(xintercept = median_val, linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(q025, q975), linetype = "dotted", color = "gray50") +
  labs(title = "repfac_65p ~ LogNormal(log(90), 0.3)",
       subtitle = sprintf("Median = %.1f, 95%% CI = [%.1f, %.1f]", 
                          median_val, q025, q975),
       x = "repfac_65p", y = "Density") +
  theme_minimal()

# 11. geno_frac ~ Beta(3, 17)
x <- seq(0, 0.5, length.out = x_points)
density <- dbeta(x, 3, 17)
mean_val <- 3/(3+17)
q025 <- qbeta(0.025, 3, 17)
q975 <- qbeta(0.975, 3, 17)

plot_list[["geno_frac"]] <- 
  ggplot(data.frame(x = x, y = density), aes(x, y)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_vline(xintercept = mean_val, linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(q025, q975), linetype = "dotted", color = "gray50") +
  labs(title = "geno_frac ~ Beta(3, 17)",
       subtitle = sprintf("Mean = %.3f, 95%% CI = [%.3f, %.3f]", 
                          mean_val, q025, q975),
       x = "geno_frac", y = "Density") +
  theme_minimal() +
  xlim(0, 0.5)

# Display all plots in a grid
grid.arrange(grobs = plot_list, ncol = 3)

# Or save individually
for (name in names(plot_list)) {
  ggsave(paste0("prior_", name, ".png"), 
         plot_list[[name]], 
         width = 6, height = 4)
}

# Print summary table
cat("\n=== PRIOR SUMMARY ===\n\n")
cat("Parameter          | Distribution                    | Mean/Median | 95% CI\n")
cat("-------------------|--------------------------------|-------------|------------------\n")
cat("factor_beta_1      | Uniform(0.15, 1)               | 0.575       | [0.15, 1]\n")
cat("factor_beta_2      | Uniform(0.15, 1)               | 0.575       | [0.15, 1]\n")
cat(sprintf("beta_3             | LogNormal(log(0.18), 0.5)      | %.3f       | [%.3f, %.3f]\n", 
            exp(log(0.18) + 0.5^2/2), 
            qlnorm(0.025, log(0.18), 0.5), 
            qlnorm(0.975, log(0.18), 0.5)))
cat("factor_beta_4      | Uniform(0.15, 1)               | 0.575       | [0.15, 1]\n")
cat("maternalAB         | Uniform(1, 365)                | 183 days    | [1, 365]\n")
cat(sprintf("aduRR              | Beta(2, 2)                     | %.3f       | [%.3f, %.3f]\n", 
            0.5, qbeta(0.025, 2, 2), qbeta(0.975, 2, 2)))
cat(sprintf("imm_yr             | Gamma(4, 2.5)                  | %.2f        | [%.2f, %.2f]\n", 
            10, qgamma(0.025, 4, scale=2.5), qgamma(0.975, 4, scale=2.5)))
cat("imm_fac            | Uniform(1, 5)                  | 3           | [1, 5]\n")
cat(sprintf("ratio_0            | LogNormal(log(2.7), 0.5)       | %.2f        | [%.2f, %.2f]\n", 
            exp(log(2.7)), qlnorm(0.025, log(2.7), 0.5), qlnorm(0.975, log(2.7), 0.5)))
cat(sprintf("ratio_5            | LogNormal(log(28), 0.5)        | %.2f        | [%.2f, %.2f]\n", 
            exp(log(28)), qlnorm(0.025, log(28), 0.5), qlnorm(0.975, log(28), 0.5)))
cat(sprintf("ratio_15           | LogNormal(log(2.9), 0.5)       | %.2f        | [%.2f, %.2f]\n", 
            exp(log(2.9)), qlnorm(0.025, log(2.9), 0.5), qlnorm(0.975, log(2.9), 0.5)))
cat(sprintf("repfac_65p         | LogNormal(log(90), 0.3)        | %.1f        | [%.1f, %.1f]\n", 
            exp(log(90)), qlnorm(0.025, log(90), 0.3), qlnorm(0.975, log(90), 0.3)))
cat(sprintf("geno_frac          | Beta(3, 17)                    | %.3f       | [%.3f, %.3f]\n", 
            0.15, qbeta(0.025, 3, 17), qbeta(0.975, 3, 17)))