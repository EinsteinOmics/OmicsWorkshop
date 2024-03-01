import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
from matplotlib import pyplot
pyplot.rcParams['figure.dpi'] = 200

# Read data into Python
hwe = pd.read_csv("plink.hwe", delim_whitespace=True)
hwe_zoom = pd.read_csv("plinkzoomhwe.hwe", delim_whitespace=True)

# Plot histogram for HWE
fig, axs = plt.subplots(1, 2, figsize=(7, 3))  # 1 rows, 2 column
axs[0].hist(hwe.iloc[:, 8], color='blue', edgecolor='black')
axs[0].set_title("Histogram HWE")
axs[0].title.set_size(8)
# Plot histogram for strongly deviating SNPs only
axs[1].hist(hwe_zoom.iloc[:, 8], color='blue', edgecolor='black')
axs[1].set_title("Histogram HWE: Strongly deviating SNPs only")
axs[1].title.set_size(8)


# Adjust layout to prevent overlapping
plt.tight_layout()

# Save the combined PDF
plt.savefig("combined_histograms.pdf")
plt.show()

