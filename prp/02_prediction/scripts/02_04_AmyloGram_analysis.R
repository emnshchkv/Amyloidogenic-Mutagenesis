# Load required libraries
library(seqinr)
library(AmyloGram)
library(tidyverse)

# --- Configuration ---
# Path to the input FASTA file
fasta_file <- "../data/prp_sequences.fasta" 
wt_name <- "WT" # Specify the exact name of the wild-type sequence

# Check if the FASTA file exists
if (!file.exists(fasta_file)) {
  stop(paste("Error: FASTA file not found at", fasta_file))
}

# --- Data Reading and Prediction ---
# Read sequences. seqinr reads sequences as a list of character vectors, 
# which is ideal for the predict() function in AmyloGram.
seqs <- read.fasta(fasta_file, seqtype = "AA")

# Run prediction
# The AmyloGram_model is automatically loaded with the package
print("Running AmyloGram prediction for all sequences. This may take a moment...")
predictions <- predict(AmyloGram_model, seqs)

# Convert results into a convenient dataframe
# The predictions object can be coerced to a data.frame,
# containing 'Name' (protein name) and 'Probability' (overall amyloidogenicity probability)
df_results <- as.data.frame(predictions) %>%
  rename(mutation = Name)

# --- Calculate Delta Probability Relative to Wild-Type ---
wt_data <- df_results %>% filter(mutation == wt_name)

if (nrow(wt_data) == 0) {
  stop(paste("Error: Data for the wild-type (", wt_name, ") was not found. Please check the 'wt_name' variable."))
}

wt_prob <- wt_data$Probability[1]

df_amylo <- df_results %>%
  filter(mutation != wt_name) %>%
  mutate(
    # In AmyloGram: higher probability = higher amyloidogenicity.
    # Therefore, a negative delta means successful disruption (reduced amyloidogenicity).
    delta_prob = Probability - wt_prob,
    effect = case_when(
      delta_prob < 0 ~ "Reduced Amyloidogenicity",
      delta_prob > 0 ~ "Increased Amyloidogenicity",
      TRUE ~ "Neutral"
    )
  ) %>%
  arrange(delta_prob) %>%
  mutate(mutation = factor(mutation, levels = mutation))

# --- Save Results to CSV and Check Directories ---
# Check if directories exist and create them if not
if (!dir.exists("../data")) {
  dir.create("../data", recursive = TRUE)
}

# Save the processed dataframe
write_csv(df_amylo, "../data/amylogram_results.csv")

# --- Waterfall Plot Generation ---
waterfall_plot_amylo <- ggplot(df_amylo, aes(x = mutation, y = delta_prob, fill = effect)) +
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
    title = "Change in Prion Amyloidogenicity (AmyloGram)",
    subtitle = "Metric: Change in probability score relative to WT",
    x = "Mutation",
    y = expression(paste(Delta, " Probability (", P[mut] - P[wt], ")"))
  )

# Display the plot
print(waterfall_plot_amylo)

# Save the plot in high resolution
if (!dir.exists("../figures")) {
  dir.create("../figures", recursive = TRUE)
}
ggsave("../figures/waterfall_plot_amylogram.png", plot = waterfall_plot_amylo, width = 14, height = 7, units = "in", dpi = 300)