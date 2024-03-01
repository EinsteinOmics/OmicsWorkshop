import pandas as pd
import numpy as np

# Read data into Python
het = pd.read_table("het_check.het", header=None)  # Set header to None

# Assuming the columns you want are the 5th and 6th (adjust as needed)
het['HET_RATE'] = (het.iloc[:, 4] - het.iloc[:, 5]) / het.iloc[:, 4]

# Filter rows based on condition
mean_het = het['HET_RATE'].mean()
std_het = het['HET_RATE'].std()
het_fail = het[(het['HET_RATE'] < mean_het - 3 * std_het) | (het['HET_RATE'] > mean_het + 3 * std_het)]

# Calculate Heterozygosity Rate Z-score
het_fail['HET_DST'] = (het_fail['HET_RATE'] - mean_het) / std_het

# Write the filtered DataFrame to a new file
het_fail.to_csv("fail-het-qc.txt", index=False, sep='\t')
