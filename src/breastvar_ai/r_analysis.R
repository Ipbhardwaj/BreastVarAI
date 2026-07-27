# R Feature Profiling & Distribution Analysis for BreastVarAI
# Validates physicochemical metrics (Grantham scores, weight deltas) across breast cancer gene variants

suppressPackageStartupMessages({
  library(stats)
})

# Load engineered dataset aggregated from Python pipeline
data_files <- list.files(path = "data", pattern = "*_engineered_features.csv", full.names = TRUE)

if (length(data_files) > 0) {
  dataset <- do.call(rbind, lapply(data_files, read.csv))
  
  cat("=== BreastVarAI R Analytics Summary ===\n")
  cat("Total Variants Processed:", nrow(dataset), "\n")
  cat("Genes Analyzed:", length(unique(dataset$gene_symbol)), "\n\n")
  
  cat("=== Grantham Score Distributions by Target Class ===\n")
  print(aggregate(grantham_score ~ target_class, data = dataset, FUN = function(x) c(mean = mean(x), sd = sd(x))))
  
  cat("\n=== Charge Delta Frequency ===\n")
  print(table(dataset$charge_delta))
  
  cat("\n[SUCCESS] R analysis script completed successfully.\n")
} else {
  cat("[WARNING] No engineered feature files found in data/ directory.\n")
}