# Load required libraries
library(tidyverse)
library(readxl)

# --- Configuration ---
# Path to the actual Excel file (.xls or .xlsx)
file_path <- "../data/prp_aggrescan.xls" 

# Exact name of the wild-type column in your file
wt_name <- "WT" 

# Check if the Excel file exists
if (!file.exists(file_path)) {
  stop(paste("Error: Excel file not found at", file_path))
}

# Read the Excel file
df_raw <- read_excel(file_path)

# --- Data Transformation (Transposing) ---
# Rename the very first column (containing metric names) to "Metric",
# since it might be named "...1" or something else during import.
colnames(df_raw)[1] <- "Metric"

df_aggre <- df_raw %>%
  # 1. Delete the column containing the averaged metric
  select(-`all sequences average`) %>%
  
  # 2. Delete ghost columns
  select(-starts_with("...")) %>%
  
  # 3. Find the row containing the "Na4vSS" metric
  filter(grepl("Na4vSS", Metric)) %>%
  
  # 4. Pivot the table: all columns except "Metric" become rows
  pivot_longer(
    cols = -Metric,          
    names_to = "mutation",   
    values_to = "Na4vSS"     
  ) %>%
  
  # 5. Ensure the value column is numeric
  mutate(Na4vSS = as.numeric(Na4vSS)) %>%
  
  # 6. Remove the Metric column as it is no longer needed
  select(-Metric)

# --- Calculate Delta Relative to WT ---
wt_data <- df_aggre %>% filter(mutation == wt_name)

if (nrow(wt_data) == 0) {
  stop("Error: Wild-type data not found! Check the exact spelling of the WT column in the Excel file.")
}

wt_score <- wt_data$Na4vSS[1]

df_aggre <- df_aggre %>%
  filter(mutation != wt_name) %>%
  mutate(
    delta_aggre = Na4vSS - wt_score,
    # In AGGRESCAN: lower Na4vSS = reduced aggregation (successful breakers)
    effect = case_when(
      delta_aggre < 0 ~ "Reduced Amyloidogenicity",
      delta_aggre > 0 ~ "Increased Amyloidogenicity",
      TRUE ~ "Neutral"
    )
  ) %>%
  arrange(delta_aggre) %>%
  mutate(mutation = factor(mutation, levels = mutation))

# --- Save Results to CSV and Check Directories ---
# Check if directories exist and create them if not
if (!dir.exists("../data")) {
  dir.create("../data", recursive = TRUE)
}

# Save the processed dataframe
write_csv(df_aggre, "../data/aggrescan_results.csv")

# --- Waterfall Plot Generation ---
waterfall_plot_aggre <- ggplot(df_aggre, aes(x = mutation, y = delta_aggre, fill = effect)) +
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
    title = "Change in Prion Amyloidogenicity (AGGRESCAN)",
    subtitle = "Metric: Change in Na4vSS score relative to WT",
    x = "Mutation",
    y = expression(paste(Delta, " Na4vSS Score (", S[mut] - S[wt], ")"))
  )

# Display the plot
print(waterfall_plot_aggre)

# Save the plot in high resolution
if (!dir.exists("../figures")) {
  dir.create("../figures", recursive = TRUE)
}
ggsave("../figures/waterfall_plot_aggrescan.png", plot = waterfall_plot_aggre, width = 14, height = 7, units = "in", dpi = 300)