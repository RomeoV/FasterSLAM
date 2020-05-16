


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# data = pd.read_csv('./benchmarks.csv')
sns.set(style="darkgrid")

# Load an example dataset with long-form data
fmri = sns.load_dataset("fmri")

# Plot the responses for different events and regions
g = sns.lineplot(x="timepoint", y="signal",
             hue="region", style="event",
             data=fmri)
g.figure.savefig("testing.png")

plt.show()