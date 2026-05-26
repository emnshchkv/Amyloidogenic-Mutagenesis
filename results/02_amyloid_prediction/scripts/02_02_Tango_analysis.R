# Load required libraries
library(tidyverse)

# --- Configuration ---
# Path to the directory with TANGO results (.txt files)
data_dir <- "../data/prp_tango" 
wt_name <- "WT"  # Specify the exact name of the wild-type file (without .txt)data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAkCAYAAAD7PHgWAAABBklEQVR4Xu2XMQrCQBBFBQvR6wgJHsEDpHVjBDvvoBhbI3bWCkZbFUyhFrYiEat0WgmC6AVkdQqbIVmWZAOi82C64b+/bDWZDEEQP4phTLMaa9d003bTGMgu1psF7JVGNzuWPdzs18GDz443rgrIcndXbvW8g1axGfZKo7P2eBXc+WB74a3FGXtiA1kwzfnpqTF7hL3SwDfAaz+BqvjkwYADe6WhglQwJlQwKVQwKakVTGOoYNL5z4JxwBlUMEwqAu9SwTCpCLxLBcOkIvCusoKT9/WFQ6OkIvCukoJwt5rO0sehUVIReBem6ng+OLBXmnKjn4PbGM5PeKnqgXIlo5vHXoL4Nl4ZYqbbEGA7+wAAAABJRU5ErkJggg==
start_pos <- 148 # Start of the target region (170 - 22)
end_pos <- 173   # End of the target region (195 - 22)

# Find all .txt files in the directory
file_paths <- list.files(path = data_dir, 
                         pattern = "\\.txt$", 
                         full.names = TRUE)

if (length(file_paths) == 0) {
  stop("Error: No .txt files were found in the specified directory.")
}

# --- Data Reading and Processing ---
results_list <- lapply(file_paths, function(file_path) {
  # Extract mutation name (remove .txt extension)
  mut_name <- str_remove(basename(file_path), "\\.txt$")
  
  # Read the data. read_table handles any amount of whitespace well
  df <- read_table(file_path, show_col_types = FALSE)
  
  # Force 'res' column to numeric format to drop leading zeros (01 -> 1)
  df <- df %>% mutate(res = as.numeric(res))
  
  # Filter data by the target region
  df_region <- df %>% filter(res >= start_pos & res <= end_pos)
  
  # Calculate the maximum aggregation value in the region
  max_agg <- max(df_region$Aggregation, na.rm = TRUE)
  
  tibble(mutation = mut_name, max_agg = max_agg)
})

# Combine results into a single dataframe
df_results <- bind_rows(results_list)

# --- Calculate Delta Relative to WT ---
wt_data <- df_results %>% filter(mutation == wt_name)

if (nrow(wt_data) == 0) {
  stop(paste("Error: Data for the wild-type (", wt_name, ") was not found. Please check the 'wt_name' variable."))
}

wt_max_agg <- wt_data$max_agg[1]

df_tango <- df_results %>%
  filter(mutation != wt_name) %>%
  mutate(
    # Calculate the difference
    delta_agg = max_agg - wt_max_agg,
    
    # Note: The logic of effects is INVERTED compared to PASTA
    effect = case_when(
      delta_agg < 0 ~ "Reduced Amyloidogenicity", # Decrease in peak
      delta_agg > 0 ~ "Increased Amyloidogenicity", # Increase in peak
      TRUE ~ "Neutral"
    )
  ) %>%
  arrange(delta_agg) %>%
  # Factorize for correct ordering on the plot
  mutate(mutation = factor(mutation, levels = mutation))

# --- Save Results to CSV and Check Directories ---
# Check if directories exist and create them if not
if (!dir.exists("../data")) {
  dir.create("../data", recursive = TRUE)
}

# Save the processed dataframe
write_csv(df_tango, "../data/tango_results.csv")

# --- Waterfall Plot Generation ---
waterfall_plot <- ggplot(df_tango, aes(x = mutation, y = delta_agg, fill = effect)) +
  geom_bar(stat = "identity", width = 0.8) +
  # Maintain the same colors, but now blue bars represent decreased aggregation
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
    title = "Change in Prion Amyloidogenicity (TANGO Algorithm)",
    subtitle = paste0("Metric: Change in peak aggregation in region ", start_pos, "-", end_pos, " relative to WT"),
    x = "Mutation",
    y = expression(paste(Delta, " Maximum Aggregation (", Agg[mut] - Agg[wt], ")"))
  )

# Display the plot
print(waterfall_plot)

# Save the plot in high resolution
if (!dir.exists("../figures")) {
  dir.create("../figures", recursive = TRUE)
}
ggsave("../figures/waterfall_plot_tango.png", plot = waterfall_plot, width = 14, height = 7, units = "in", dpi = 300)