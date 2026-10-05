# --- Step 1: Set your working directory (optional but recommended) ---
# This tells R where to look for your 'fst_results.csv' file
# Replace 'path/to/your/data' with the actual folder path where you saved the CSV.
setwd("/Users/Charlotte/Documents/Stellenbosch_University/Masters/Ch1 writing")

# --- Step 2: Load your data into R ---
# 'read.csv()' is used to import CSV files.
# 'header=TRUE' means the first row of your CSV contains column names.
fst_data <- read.csv("Full_FST_for_csv.csv", header = TRUE)

#THEN NEUTRAL (Re-run all code)
fst_data <- read.csv("Neutral_FST_for_csv.csv", header = TRUE)

#THEN OUTLIER (Re-run all code)
fst_data <- read.csv("Outlier_FST_for_csv.csv", header = TRUE)

# --- Step 3: View the loaded data (optional, for verification) ---
# This shows the first few rows of your data frame
head(fst_data)

# This shows the structure of your data frame, including column types
str(fst_data)

# --- Step 4: Extract the p-values ---
# Assuming your p-value column is named 'P_Value'
# If it's named something else, change 'P_Value' accordingly
p_values <- fst_data$P_value

# --- Step 5: Apply the Benjamini-Hochberg FDR correction ---
# 'p.adjust()' is the core function for multiple testing corrections.
# 'method = "fdr"' specifies the Benjamini-Hochberg False Discovery Rate correction.
# 'alpha' is your chosen significance level (e.g., 0.05)
fdr_corrected_p_values <- p.adjust(p_values, method = "fdr")

# --- Step 6: Add the FDR corrected p-values back to your data frame ---
fst_data$FDR_P_Value <- fdr_corrected_p_values

# --- Step 7: Identify significant results ---
# You can set your desired significance level (e.g., 0.05) for the FDR-corrected p-values
alpha_level <- 0.05
fst_data$Significant_FDR <- fst_data$FDR_P_Value <= alpha_level

# --- Step 8: View the results ---
# Display the data frame with the new columns
print(fst_data)

# You can also filter to see only the significant comparisons
significant_comparisons <- fst_data[fst_data$Significant_FDR, ]
print("Significant comparisons after FDR correction:")
print(significant_comparisons)

# --- Optional: Export the results to a new CSV file ---
write.csv(fst_data, "fst_results_with_fdr.csv", row.names = FALSE)