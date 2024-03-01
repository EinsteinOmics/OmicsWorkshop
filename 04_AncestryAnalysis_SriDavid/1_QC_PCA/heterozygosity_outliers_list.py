import pandas as pd
import warnings
warnings.filterwarnings('ignore')

# Read data into Python
het = pd.read_table("results/het_check.het", header=None, delim_whitespace=True)

# Update column names based on the actual DataFrame structure
new_column_names = ['FID', 'IID', 'O(HOM)', 'E(HOM)', 'N(NM)', 'F']
het.columns = new_column_names

# Convert columns to numeric data type
het['N(NM)'] = pd.to_numeric(het['N(NM)'], errors='coerce')
het['O(HOM)'] = pd.to_numeric(het['O(HOM)'], errors='coerce')

# Calculate Heterozygosity Rate
het['HET_RATE'] = (het['N(NM)'] - het['O(HOM)']) / het['N(NM)']

# Filter rows based on condition
mean_het = het['HET_RATE'].mean()
std_het = het['HET_RATE'].std()
het_fail = het[(het['HET_RATE'] < mean_het - 3 * std_het) | (het['HET_RATE'] > mean_het + 3 * std_het)].copy()

# Calculate Heterozygosity Rate Z-score
het_fail['HET_DST'] = (het_fail['HET_RATE'] - mean_het) / std_het

# Write the filtered DataFrame to a new file
het_fail.to_csv("txt_files/fail-het-qc.txt", index=False, sep='\t')

print(het_fail)
