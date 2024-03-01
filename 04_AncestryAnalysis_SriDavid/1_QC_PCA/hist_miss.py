import pandas as pd
import matplotlib.pyplot as plt

# Read data into Python
indmiss = pd.read_table("plink.imiss", header=None, skiprows=1, delim_whitespace=True)  # Skip the header
snpmiss = pd.read_table("plink.lmiss", header=None, skiprows=1, delim_whitespace=True)  # Skip the header

# Plot Histogram for F_MISS in plink.imiss and plink.lmiss
plt.figure(figsize=(12, 6))

# Plot Histogram for F_MISS in plink.imiss
plt.subplot(1, 2, 1)
if len(indmiss.columns) >= 6:
    plt.hist(indmiss.iloc[:, 5], bins=20, color='blue', edgecolor='black')
    plt.xlabel("Frequency of individual missingness rates")
    plt.ylabel("Count")
    plt.title("Histogram individual missingness")
else:
    print("DataFrame 'indmiss' has fewer than 6 columns. Check the structure of your data.")

# Plot Histogram for F_MISS in plink.lmiss
plt.subplot(1, 2, 2)
if len(snpmiss.columns) >= 5:
    plt.hist(snpmiss.iloc[:, 4], bins=20, color='green', edgecolor='black')
    plt.xlabel("Frequency of SNP-based missingness rates ")
    plt.ylabel("Count")
    plt.title("Histogram SNP missingness")
else:
    print("DataFrame 'snpmiss' has fewer than 5 columns. Check the structure of your data.")

# Adjust layout and save the plot
plt.tight_layout()
plt.savefig("histimiss_lmiss_combined.pdf")
plt.show()
