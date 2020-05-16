
# go to src/visualizations
# pipenv shell
# python benchmark_graphs.py
# benchmark.csv is at .src/build/core/benchmarks/benchmark.csv

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#data = pd.read_csv('../src/build/core/benchmarks/benchmark.csv')
data = pd.read_csv('./benchmark_sample.csv').values
# format: benchmark name, function name, flops, cycles, performance, speedup 
#print(data)

plot_arrays = []
benchmark_names = []

prev_array = []
prev_key = ""
for x in data:
  key = x[0].split(";")[0]
  if(prev_key == ""):
    prev_key = key
    prev_array.append(x[0].split(";")[1:])
  elif(key != prev_key):
    # print prev array
    print(prev_key, prev_array)
    plot_arrays.append(prev_array)
    benchmark_names.append(prev_key)

    prev_array = []
    prev_array.append(x[0].split(";")[1:])
    prev_key = key
  else:
    prev_array.append(x[0].split(";")[1:])
print(prev_key, prev_array)
plot_arrays.append(prev_array)
benchmark_names.append(prev_key)

sns.set(style="darkgrid")

for index, benchmark in enumerate(plot_arrays):
    titanic = pd.DataFrame(data=plot_arrays[index],
                index=range(0, len(plot_arrays[index])),
                columns=['function','flops','cycles','performance','speedup',''])
    print(titanic)
    sns.catplot(x="function", y="cycles", hue="cycles", kind="bar", data=titanic);
    plt.show()

# Load an example dataset with long-form data
fmri = sns.load_dataset("fmri")

# Plot the responses for different events and regions
g = sns.lineplot(x="timepoint", y="signal",
             hue="region", style="event",
             data=fmri)
g.figure.savefig("testing.png")

plt.show()