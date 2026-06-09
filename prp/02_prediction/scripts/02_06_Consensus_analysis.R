# Load required libraries
library(tidyverse)

# --- Configuration ---
# Set the directory where the intermediate CSV files are stored
data_dir <- "../data"

# Define expected file paths
pasta_file <- file.path(data_dir, "pasta2_results.csv")
tango_file <- file.path(data_dir, "tango_results.csv")
appnn_file <- file.path(data_dir, "appnn_results.csv")
amylo_file <- file.path(data_dir, "amylogram_results.csv")
aggre_file <- file.path(data_dir, "aggrescan_results.csv")

# Check if all required files exist before proceeding
required_files <- c(pasta_file, tango_file, appnn_file, amylo_file, aggre_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(paste("Error: Missing required CSV files. Please run the respective scripts first. Missing:",
             paste(basename(missing_files), collapse = ", ")))
}

# --- Data Reading and Preparation ---
# Read the data and select only the required columns (mutation and delta)
df_pasta <- read_csv(pasta_file, show_col_types = FALSE) %>% select(mutation, delta_pasta = delta_energy)
df_tango <- read_csv(tango_file, show_col_types = FALSE) %>% select(mutation, delta_tango = delta_agg)
df_appnn <- read_csv(appnn_file, show_col_types = FALSE) %>% select(mutation, delta_appnn = delta_mean)
df_amylo <- read_csv(amylo_file, show_col_types = FALSE) %>% select(mutation, delta_amylo = delta_prob)
df_aggre <- read_csv(aggre_file, show_col_types = FALSE) %>% select(mutation, delta_aggre = delta_aggre)

# Combine all dataframes into a single consensus table
df_consensus <- list(df_pasta, df_tango, df_amylo, df_aggre, df_appnn) %>%
  reduce(full_join, by = "mutation")

# --- Alignment and Z-standardization ---
df_scored <- df_consensus %>%
  mutate(
    # Align directions: for PASTA, invert the sign so that negative = breaker
    aligned_pasta = -delta_pasta,
    aligned_tango = delta_tango,
    aligned_amylo = delta_amylo,
    aligned_aggre = delta_aggre,
    aligned_appnn = delta_appnn
  ) %>%
  mutate(
    # Apply scale() for Z-standardization of each column
    z_pasta = c(scale(aligned_pasta)),
    z_tango = c(scale(aligned_tango)),
    z_amylo = c(scale(aligned_amylo)),
    z_aggre = c(scale(aligned_aggre)),
    z_appnn = c(scale(aligned_appnn))
  ) %>%
  # Calculate Integrative Score (mean of Z-scores)
  # na.rm = TRUE ensures calculation even if a tool missed a mutation
  rowwise() %>%
  mutate(
    Integrative_Score = mean(c(z_pasta, z_tango, z_amylo, z_aggre, z_appnn), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(Integrative_Score) %>%
  # Factorize for correct ordering on plots
  mutate(mutation = factor(mutation, levels = mutation))

# --- Save Results to CSV and Check Directories ---
# Check if directories exist and create them if not
if (!dir.exists("../data")) {
  dir.create("../data", recursive = TRUE)
}
if (!dir.exists("../figures")) {
  dir.create("../figures", recursive = TRUE)
}

# Save the consensus dataframe for the repository and further analysis
write_csv(df_scored, "../data/consensus_results.csv")

# --- Visualization ---

# Plot 1: Waterfall plot of the integrative score
waterfall_consensus <- ggplot(df_scored, aes(x = mutation, y = Integrative_Score,
                                             fill = Integrative_Score < 0)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), guide = "none") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Consensus Analysis of Amyloidogenicity (5 Tools)",
    subtitle = "Integrative Score (Mean Z-score)",
    x = "Mutation",
    y = "Integrative Z-Score"
  )

# Plot 2: Heatmap of the top 30 best mutations
# The heatmap allows checking if a mutation is recognized as a breaker by ALL tools,
# or if it was pulled up by one extreme value.
top_30_muts <- head(df_scored$mutation, 30)

df_heatmap <- df_scored %>%
  filter(mutation %in% top_30_muts) %>%
  select(mutation, PASTA2.0 = z_pasta, TANGO = z_tango, AmyloGram = z_amylo,
         AGGRESCAN = z_aggre, APPNN = z_appnn) %>%
  pivot_longer(cols = -mutation, names_to = "Tool", values_to = "Z_Score")

heatmap_plot <- ggplot(df_heatmap, aes(x = Tool, y = mutation, fill = Z_Score)) +
  geom_tile(color = "white") +
  # Blue = reduced aggregation (good), Red = increased aggregation (bad)
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C62828", midpoint = 0) +
  theme_minimal() +
  labs(
    title = "Contribution of Each Tool for Top-30 Candidates",
    x = "Prediction Tool",
    y = "Mutation",
    fill = "Z-Score"
  )

# Display plots
print(waterfall_consensus)
print(heatmap_plot)

# Save plots in high resolution
ggsave("../figures/06_consensus_waterfall.png", plot = waterfall_consensus, width = 14, height = 7, units = "in", dpi = 300)
ggsave("../figures/06_consensus_heatmap.png", plot = heatmap_plot, width = 8, height = 10, units = "in", dpi = 300)
