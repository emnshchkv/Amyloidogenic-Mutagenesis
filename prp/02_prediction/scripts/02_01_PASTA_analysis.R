# Load required libraries
library(tidyverse)

# --- Configuration ---
# Path to the directory with unzipped PASTA 2.0 results
data_dir <- "../data/prp_pasta20"
wt_name <- "WT" # Specify the exact name of the wild-type sequence
start_pos <- 148 # Start of the region (170 - 22)
end_pos <- 173   # End of the region (195 - 22)

# Find all files containing the free energy profile
file_paths <- list.files(path = data_dir,
                         pattern = "\\.aggr_profile\\.dat\\.free_energy$",
                         full.names = TRUE)

if (length(file_paths) == 0) {
  stop("Error: No .aggr_profile.dat.free_energy files were found in the specified directory.")
}

# --- Data Reading and Processing ---
results_list <- lapply(file_paths, function(file_path) {
  # Extract mutation name from the file name
  file_name <- basename(file_path)
  mut_name <- str_remove(file_name, "\\.fasta.*$")

  # Read the single column containing energy values
  df <- read_table(file_path, col_names = c("energy"), show_col_types = FALSE)

  # Add sequence position based on the row number
  df <- df %>% mutate(position = row_number())

  # Filter data to keep only the target region
  df_region <- df %>% filter(position >= start_pos & position <= end_pos)

  # Calculate required metrics for the region
  min_energy <- min(df_region$energy, na.rm = TRUE)
  mean_energy <- mean(df_region$energy, na.rm = TRUE)

  tibble(mutation = mut_name, min_energy = min_energy, mean_energy = mean_energy)
})

# Combine all results into a single dataframe
df_results <- bind_rows(results_list)

# --- Calculate Delta Energy Relative to Wild-Type ---
wt_data <- df_results %>% filter(mutation == wt_name)

if (nrow(wt_data) == 0) {
  stop(paste("Error: Data for the wild-type (", wt_name, ") was not found. Please check the 'wt_name' variable."))
}

wt_min_energy <- wt_data$min_energy[1]

df_pasta <- df_results %>%
  filter(mutation != wt_name) %>%
  mutate(
    delta_energy = min_energy - wt_min_energy,
    effect = case_when(
      delta_energy > 0 ~ "Reduced Amyloidogenicity",
      delta_energy < 0 ~ "Increased Amyloidogenicity",
      TRUE ~ "Neutral"
    )
  ) %>%
  arrange(delta_energy) %>%
  mutate(mutation = factor(mutation, levels = mutation))

# --- Save Results to CSV ---
if (!dir.exists("../data")) {
  dir.create("../data", recursive = TRUE)
}
# Save the processed dataframe with delta energies to a CSV file
write_csv(df_pasta, "../data/pasta2_results.csv")

# --- Waterfall Plot Generation ---
waterfall_plot <- ggplot(df_pasta, aes(x = mutation, y = delta_energy, fill = effect)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("Reduced Amyloidogenicity" = "blue",
                               "Increased Amyloidogenicity" = "red",
                               "Neutral" = "#9E9E9E")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  labs(
    title = "Change in Prion Amyloidogenicity Upon Point Mutations",
    subtitle = paste0("Metric: Minimum free energy in region ", start_pos, "-", end_pos, " (relative to WT)"),
    x = "Mutation",
    y = expression(paste(Delta, " Minimum Free Energy (", E[mut] - E[wt], ")"))
  )

# Display the plot
print(waterfall_plot)

# Save the plot in high resolution
if (!dir.exists("../figures")) {
  dir.create("../figures", recursive = TRUE)
}
ggsave("../figures/waterfall_plot_pasta.png", plot = waterfall_plot, width = 14, height = 7, units = "in", dpi = 300)
