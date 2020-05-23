
# go to src/visualizations
# pipenv shell
# python benchmark_graphs.py
# benchmark.csv is at .src/build/core/benchmarks/benchmark.csv

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.cbook import flatten

#data = pd.read_csv('../src/build/core/benchmarks/benchmark.csv')
data = pd.read_csv('./polybox/benchmark.csv').values
# format: benchmark name, function name, flops, cycles, performance, speedup 
#print(data)

plot_arrays = []
benchmark_names = []
optim = []

prev_array = []
prev_key = ""
for i, x in enumerate(data):
  key = x[0].split(";")[0]
  if(prev_key == ""):
    prev_key = key
    arr = x[0].split(";")
    arr[3] = float(arr[3])
    arr[4] = float(arr[4])
    arr[5] = float(arr[5])
    prev_array.append(arr)
  elif(key != prev_key):
    # print prev array
    #print(prev_key, prev_array)
    optim1 = list(prev_array[0])
    optim1[0] = len(benchmark_names)
    optim1[1] = "base"
    optim2 = list(prev_array[1])
    optim2[0] = len(benchmark_names)
    optim2[1] = "ours"
    optim = optim + [optim1] + [optim2]
    plot_arrays.append(prev_array)
    benchmark_names.append(prev_key)

    prev_array = []
    arr = x[0].split(";")
    arr[3] = float(arr[3])
    arr[4] = float(arr[4])
    arr[5] = float(arr[5])
    prev_array.append(arr)
    prev_key = key
  else:
    arr = x[0].split(";")
    arr[3] = float(arr[3])
    arr[4] = float(arr[4])
    arr[5] = float(arr[5])
    prev_array.append(arr)

print(prev_key, prev_array)
plot_arrays.append(prev_array)
benchmark_names.append(prev_key)

sns.set(style="darkgrid")

for index, benchmark in enumerate(plot_arrays):
    print(plot_arrays[index])
    plot_data = pd.DataFrame(data=plot_arrays[index],
                columns=['bench','function','','flops','bytes','cycles','performance','speedup'])
    print(plot_data)
    sns.catplot(x="bench", y="cycles", hue="function", kind="bar", data=plot_data)
    plt.xlabel("")
    plt.savefig("./plots/runtime/{0}_benchmark.png".format(benchmark_names[index]))
    plt.show()

print(optim)
plot_data = pd.DataFrame(data=optim,
            columns=['bench','function','','flops','bytes','cycles','performance','speedup'])
print(plot_data)
sns.catplot(x="bench", y="speedup", hue="function", kind="bar", data=plot_data)
plt.xlabel("")
plt.savefig("./plots/speedup/speedup_benchmark.png")
plt.show()
