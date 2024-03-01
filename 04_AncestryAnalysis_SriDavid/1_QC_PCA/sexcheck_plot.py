import pandas as pd
import matplotlib.pyplot as plt

# Read data into Python
gender = pd.read_table("plink.sexcheck", delim_whitespace=True)

# Create a single figure with multiple subplots
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(8, 12))

# Plot Sex Check histogram
axes[0].hist(gender["F"], bins=20, color='blue', edgecolor='black')
axes[0].set_title("Sex Check")
axes[0].set_xlabel("F")

# Plot Men histogram
male = gender[gender["PEDSEX"] == 1]
axes[1].hist(male["F"], bins=20, color='blue', edgecolor='black')
axes[1].set_title("Male")
axes[1].set_xlabel("F")

# Plot Women histogram
female = gender[gender["PEDSEX"] == 2]
axes[2].hist(female["F"], bins=20, color='blue', edgecolor='black')
axes[2].set_title("Female")
axes[2].set_xlabel("F")

# Adjust layout and save the plot
plt.tight_layout()
plt.savefig("Gender_Checks.pdf")
plt.show()
