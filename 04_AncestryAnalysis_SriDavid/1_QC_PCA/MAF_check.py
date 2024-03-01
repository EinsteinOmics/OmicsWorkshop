import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
from matplotlib import pyplot
pyplot.rcParams['figure.dpi'] = 200

# Read the data into a DataFrame
maf_freq = pd.read_table("MAF_check.frq", delim_whitespace=True, skiprows=1, names=["CHR", "SNP", "A1", "A2", "MAF", "NCHROBS"])

# Extract MAF values from the 5th column
maf_values = pd.to_numeric(maf_freq["MAF"], errors='coerce')

# Create a histogram
plt.figure(figsize=(5,3))
plt.hist(maf_values.dropna(), edgecolor='black')  # Adjust the number of bins as needed
plt.title("MAF Distribution")
plt.xlabel("MAF")
plt.ylabel("# of SNPs")
plt.show()
