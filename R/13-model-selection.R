# Load required packages
library(coda)
library(posterior)
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)

# ========================================
# SETUP: Define paths and model names
# ========================================
root           <- here::here()
source(file.path(root, "R", "modify_attach.R"))
# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))
# Set the directory where your RDS files are stored
results_dir <- file.path(root, "output","final","/")  # UPDATE THIS PATH

# Define model names
model_names <- c("Model_A_5pct", "Model_B_5pct", "Model_C_5pct",
                 "Model_A_25pct", "Model_B_25pct", "Model_C_25pct",
                 "Model_A_50pct", "Model_B_50pct", "Model_C_50pct")

# ========================================
# 1. DEVIANCE INFORMATION CRITERION (DIC)
# ========================================

# Function to calculate DIC from MCMC results
calculate_dic <- function(results) {
  # Calculate deviance
  # If results$density contains log-likelihood, convert to deviance
  D <- -2 * results$density
  
  # Make sure D is a vector (flatten if it's a matrix)
  D <- as.vector(D)
  
  # Calculate DIC using Gelman et al. formulation
  # pD = var(D)/2 (effective number of parameters)
  pD <- var(D) / 2
  
  # D_bar = mean(D)
  D_bar <- mean(D)
  
  # DIC = pD + D_bar
  DIC <- pD + D_bar
  
  return(list(
    DIC = DIC, 
    D_bar = D_bar, 
    pD = pD
  ))
}

# Recalculate DIC for all models
dic_results <- list()
for (model_name in model_names) {
  # Load results from RDS file
  results_file <- paste0(results_dir, "results_", model_name, ".rds")
  results <- readRDS(results_file)
  
  dic_results[[model_name]] <- calculate_dic(results)
  cat("Processed DIC for:", model_name, 
      "- DIC =", round(dic_results[[model_name]]$DIC, 2),
      "- pD =", round(dic_results[[model_name]]$pD, 2), "\n")
}

# Create DIC comparison data frame
dic_df <- data.frame(
  Model = names(dic_results),
  DIC = sapply(dic_results, function(x) x$DIC),
  D_bar = sapply(dic_results, function(x) x$D_bar),
  pD = sapply(dic_results, function(x) x$pD),
  stringsAsFactors = FALSE
)

# Add cross-protection level and model type
dic_df$CrossProtection <- case_when(
  grepl("_5pct", dic_df$Model) ~ "5%",
  grepl("_25pct", dic_df$Model) ~ "25%",
  grepl("_50pct", dic_df$Model) ~ "50%"
)

dic_df$ModelType <- case_when(
  grepl("Model_A", dic_df$Model) ~ "Model A",
  grepl("Model_B", dic_df$Model) ~ "Model B",
  grepl("Model_C", dic_df$Model) ~ "Model C"
)

# Sort by DIC
dic_df <- dic_df[order(dic_df$DIC), ]

print(dic_df)

# Show delta DIC relative to best model
dic_df$delta_DIC <- dic_df$DIC - min(dic_df$DIC)
print(dic_df[, c("Model", "DIC", "delta_DIC", "pD")])
# ========================================
# 2. DIC VISUALIZATION
# ========================================

# Bar chart of DIC values
p1 <- ggplot(dic_df, aes(x = ModelType, y = DIC, fill = CrossProtection)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = round(DIC, 1)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  theme_minimal() +
  labs(title = "Deviance Information Criterion (DIC) Comparison",
       subtitle = "Lower DIC indicates better model fit",
       x = "Model", y = "DIC", fill = "Cross-Protection") +
  theme(legend.position = "bottom")

print(p1)

# Alternative: grouped bar chart with delta DIC
best_dic <- min(dic_df$DIC)
dic_df$delta_DIC <- dic_df$DIC - best_dic

p2 <- ggplot(dic_df, aes(x = ModelType, y = delta_DIC, fill = CrossProtection)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = c(2, 7), linetype = "dashed", color = "red", alpha = 0.5) +
  geom_text(aes(label = round(delta_DIC, 1)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  theme_minimal() +
  labs(title = "Delta DIC (Relative to Best Model)",
       subtitle = "ΔDIC > 7 indicates substantially worse fit",
       x = "Model", y = "ΔDIC", fill = "Cross-Protection") +
  theme(legend.position = "bottom")

print(p2)

# ========================================
# 3. RHAT DIAGNOSTIC TABLE
# ========================================

# Function to extract Rhat from posterior package
extract_rhat <- function(sample_df) {
  draws <- as_draws(sample_df)
  rhat_vals <- summarise_draws(draws, "rhat")
  return(rhat_vals)
}

# Create Rhat table for all models
rhat_list <- list()
for (model_name in model_names) {
  # Load sample_df from RDS file
  sample_file <- paste0(results_dir, "sample_df_", model_name, ".qs2")
  sample_df <- qs_read(sample_file)
  
  rhat_list[[model_name]] <- extract_rhat(sample_df)
  cat("Processed Rhat for:", model_name, "\n")
}

# Combine into single table
rhat_table <- bind_rows(rhat_list, .id = "Model")

# Reshape for better visualization
rhat_wide <- rhat_table %>%
  select(Model, variable, rhat) %>%
  pivot_wider(names_from = Model, values_from = rhat)

# Flag problematic Rhat values (> 1.01)
rhat_summary <- rhat_table %>%
  group_by(Model) %>%
  summarise(
    max_rhat = max(rhat, na.rm = TRUE),
    n_problematic = sum(rhat > 1.01, na.rm = TRUE),
    problematic_params = paste(variable[rhat > 1.01], collapse = ", ")
  )

print(rhat_summary)

# Create formatted table
kable(rhat_wide, digits = 4, caption = "Rhat Diagnostics by Model") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"))

# ========================================
# 4. POSTERIOR PARAMETER TABLE WITH 95% CrI
# ========================================

# Function to extract posterior summaries and transform log parameters
extract_posterior_summary <- function(sample_df, model_name) {
  draws <- as_draws(sample_df)

  # Get summary statistics
  summary_stats <- summarise_draws(draws, 
                                   median = median,
                                   q2.5 = ~quantile(.x, 0.025),
                                   q97.5 = ~quantile(.x, 0.975))
  
  # Rename columns for easier handling (the quantile columns will be named "2.5%" and "97.5%")
  names(summary_stats)[names(summary_stats) == "2.5%"] <- "q2.5"
  names(summary_stats)[names(summary_stats) == "97.5%"] <- "q97.5"
  
  # Transform log parameters
  summary_transformed <- summary_stats %>%
    mutate(
      median_trans = ifelse(grepl("^log_", variable), exp(median), median),
      q2.5_trans = ifelse(grepl("^log_", variable), exp(q2.5), q2.5),
      q97.5_trans = ifelse(grepl("^log_", variable), exp(q97.5), q97.5),
      # Create formatted string with CrI
      posterior_summary = sprintf("%.3f (%.3f, %.3f)", 
                                  median_trans, q2.5_trans, q97.5_trans),
      # Clean variable name (remove log_ prefix for display)
      param_name = gsub("^log_", "", variable)
    ) %>%
    select(param_name, posterior_summary)
  
  return(summary_transformed)
}

# Extract posteriors for all models
posterior_list <- list()
for (model_name in model_names) {
  # Load sample_df from RDS file
  sample_file <- paste0(results_dir, "sample_df_", model_name, ".qs2")
  sample_df <- qs_read(sample_file)
  
  posterior_list[[model_name]] <- extract_posterior_summary(sample_df, model_name)
  cat("Processed posteriors for:", model_name, "\n")
}

# Combine into wide format table
posterior_combined <- bind_rows(posterior_list, .id = "Model")

# Parse Model column to get ModelType and CrossProtection
posterior_combined <- posterior_combined %>%
  mutate(
    CrossProtection = case_when(
      grepl("_5pct", Model) ~ "5%",
      grepl("_25pct", Model) ~ "25%",
      grepl("_50pct", Model) ~ "50%"
    ),
    ModelType = case_when(
      grepl("Model_A", Model) ~ "Model A",
      grepl("Model_B", Model) ~ "Model B",
      grepl("Model_C", Model) ~ "Model C"
    )
  )

# Create wide format table matching your template
posterior_wide <- posterior_combined %>%
  select(param_name, ModelType, CrossProtection, posterior_summary) %>%
  pivot_wider(
    names_from = c(CrossProtection, ModelType),
    values_from = posterior_summary,
    names_sep = "_"
  ) %>%
  # Reorder columns to match template
  select(param_name, 
         `5%_Model A`, `5%_Model B`, `5%_Model C`,
         `25%_Model A`, `25%_Model B`, `25%_Model C`,
         `50%_Model A`, `50%_Model B`, `50%_Model C`)

# Handle Imm_fac (only in Model B)
# If parameter doesn't exist for a model, it will be NA

# Create formatted table
kable(posterior_wide, 
      col.names = c("Parameter",
                    rep(c("Model A", "Model B", "Model C"), 3)),
      caption = "Posterior Medians with 95% Credible Intervals (transformed for log-scale parameters)") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  add_header_above(c(" " = 1, 
                     "5% Cross-protection" = 3, 
                     "25% Cross-protection" = 3, 
                     "50% Cross-protection" = 3))

# Save table to CSV
write.csv(posterior_wide, "posterior_summary_table.csv", row.names = FALSE)

# ========================================
# EXPORT RESULTS
# ========================================

# Save DIC results
write.csv(dic_df, paste0(results_dir, "dic_comparison.csv"), row.names = FALSE)

# Save plots
ggsave(paste0(results_dir, "dic_barchart.png"), p1, width = 10, height = 6)
ggsave(paste0(results_dir, "delta_dic_barchart.png"), p2, width = 10, height = 6)

# Save Rhat summary
write.csv(rhat_summary, paste0(results_dir, "rhat_summary.csv"), row.names = FALSE)
write.csv(rhat_wide, paste0(results_dir, "rhat_detailed.csv"), row.names = FALSE)

# Save posterior summary table
write.csv(posterior_wide, paste0(results_dir, "posterior_summary_table.csv"), row.names = FALSE)

cat("\nAll results saved to:", results_dir, "\n")
cat("Files created:\n")
cat("  - dic_comparison.csv\n")
cat("  - dic_barchart.png\n")
cat("  - delta_dic_barchart.png\n")
cat("  - rhat_summary.csv\n")
cat("  - rhat_detailed.csv\n")
cat("  - posterior_summary_table.csv\n")