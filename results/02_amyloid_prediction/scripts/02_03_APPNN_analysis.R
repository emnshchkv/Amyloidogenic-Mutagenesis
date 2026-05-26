# Load required libraries
library(seqinr)
library(appnn)
library(tidyverse)

# --- Configuration ---
# Path to the input FASTA file
fasta_file <- "../data/prp_sequences.fasta"
wt_name <- "WT"
start_pos <- 148 # Start of the target region (170 - 22)
end_pos <- 173   # End of the target region (195 - 22)

# Check if the FASTA file exists
if (!file.exists(fasta_file)) {
  stop(paste("Error: FASTA file not found at", fasta_file))
}

# --- Data Reading and Processing ---
# Read multi-FASTA file
fasta_data <- read.fasta(fasta_file, seqtype = "AA", as.string = TRUE)

# Extract sequences and their names
sequences <- as.character(fasta_data)
seq_names <- getName(fasta_data)

# Run APPNN for all sequences
print("Running APPNN for all sequences. This may take a few minutes...")
all_results <- appnn(sequences)

# --- Calculate Wild-Type (WT) Baseline ---
# Find WT sequence index
wt_idx <- which(seq_names == wt_name)
if (length(wt_idx) == 0) {
  warning("Warning: Wild-type sequence not found by name. Using the first sequence as WT.")
  wt_idx <- 1 
}

# Extract amyloid propensity and calculate the mean for the target window
wt_aminoacids <- all_results[[wt_idx]]$aminoacids
wt_local_window <- wt_aminoacids[start_pos:end_pos]
wt_mean <- mean(wt_local_window, na.rm = TRUE)

# --- Process All Mutations ---
# Using sequential processing (for loop) to ensure Windows compatibility
summary_list <- vector("list", length(all_results))

for (i in seq_along(all_results)) {
  current_aminoacids <- all_results[[i]]$aminoacids
  
  # Check length to avoid out-of-bounds errors (APPNN output might be shorter)
  window_end <- min(end_pos, length(current_aminoacids))
  local_window <- current_aminoacids[start_pos:window_end]
  
  if (length(local_window) > 0) {
    m_max <- max(local_window, na.rm = TRUE)
    m_mean <- mean(local_window, na.rm = TRUE)
    
    # Delta: negative value means a decrease in aggregation potential
    delta_mean <- m_mean - wt_mean
    
    # Store results in a tibble for the list
    summary_list[[i]] <- tibble(
      mutation = seq_names[i],
      local_max = m_max,
      local_mean = m_mean,
      delta_mean = delta_mean
    )
  }
}

# Combine all results into a single dataframe
df_appnn <- bind_rows(summary_list)

# --- Assign Effects and Sort ---
df_appnn <- df_appnn %>%
  filter(mutation != wt_name) %>%
  mutate(
    effect = case_when(
      delta_mean < 0 ~ "Reduced Amyloidogenicity",
      delta_mean > 0 ~ "Increased Amyloidogenicity",
      TRUE ~ "Neutral"
    )
  ) %>%
  arrange(delta_mean) %>%
  mutate(mutation = factor(mutation, levels = mutation))

# Print top 10 mutations
print("Top 10 mutations decreasing mean propensity:")
print(head(df_appnn, 10))

# --- Save Results to CSV and Check Directories ---
# Check if directories exist and create them if not
if (!dir.exists("../data")) {
  dir.create("../data", recursive = TRUE)
}

# Save the processed dataframe
write_csv(df_appnn, "../data/appnn_results.csv")

# --- Waterfall Plot Generation ---
waterfall_plot <- ggplot(df_appnn, aes(x = mutation, y = delta_mean, fill = effect)) +
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
    title = "Change in Prion Amyloidogenicity (APPNN Algorithm)",
    subtitle = paste0("Metric: Change in mean propensity in region ", start_pos, "-", end_pos, " relative to WT"),
    x = "Mutation",
    y = expression(paste(Delta, " Mean Propensity (", P[mut] - P[wt], ")"))
  )

# Display the plot
print(waterfall_plot)

# Save the plot in high resolution
if (!dir.exists("../figures")) {
  dir.create("../figures", recursive = TRUE)
}
ggsave("../figures/waterfall_plot_appnn.png", plot = waterfall_plot, width = 14, height = 7, units = "in", dpi = 300)