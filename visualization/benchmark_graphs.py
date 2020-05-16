
# go to src/visualizations
# pipenv shell
# python benchmark_graphs.py
# benchmark.csv is at .src/build/core/benchmarks/benchmark.csv

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#data = pd.read_csv('../src/build/core/benchmarks/benchmark.csv')
data = pd.read_csv('./benchmark_sample.csv')
print(data)

sns.set(style="darkgrid")

# Load an example dataset with long-form data
fmri = sns.load_dataset("fmri")

# Plot the responses for different events and regions
g = sns.lineplot(x="timepoint", y="signal",
             hue="region", style="event",
             data=fmri)
g.figure.savefig("testing.png")

plt.show()