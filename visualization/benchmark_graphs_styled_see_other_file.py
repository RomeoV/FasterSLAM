
# go to src/visualizations
# pipenv shell
# python benchmark_graphs.py
# benchmark.csv is at .src/build/core/benchmarks/benchmark.csv

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.cbook import flatten
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import matplotlib

## SETUP MATPLOTLIB
plt.style.use('seaborn')
matplotlib.rcParams['axes.labelweight'] = 'bold'
matplotlib.rcParams['axes.titleweight'] = 'bold'
matplotlib.rcParams['font.weight'] = 'semibold'
sns.set_color_codes()
sns.set(style="darkgrid")

#data = pd.read_csv('../src/build/core/benchmarks/benchmark.csv')
data = pd.read_csv('./benchmark_sample.csv').values
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
    prev_array.append(arr)
    prev_key = key
  else:
    arr = x[0].split(";")
    arr[3] = float(arr[3])
    prev_array.append(arr)

print(prev_key, prev_array)
plot_arrays.append(prev_array)
benchmark_names.append(prev_key)

def setup_axis(ax):
    ax.xaxis.set_label_coords(1, -0.075)
    ax.yaxis.set_label_coords(0, 1.005)
    ax.set_ylabel(ax.get_ylabel(), rotation='horizontal')
    ax.labelweight = 'bold'
    ax.set_xlabel("function\n", horizontalalignment='right')
    ax.set_ylabel("cycles\n", rotation='horizontal', horizontalalignment='left')
    ax.get_xaxis().set_major_formatter(StrMethodFormatter('{x:3.2f}'))
    ax.get_yaxis().set_major_formatter(StrMethodFormatter('{x:2.2f}'))
    ax.tick_params('x', direction='in')
    ax.legend(loc='lower left')
    return ax

def plot_in_style(elements):
    def new_plot(f):
      for i, el in enumerate(elements):
        fig, ax = plt.subplots(figsize=(6,4))
        f(el, ax)
        setup_axis(ax)
        fig.tight_layout()
        print(i)
        print(benchmark_names[i])
        fig.savefig("./{0}_performance.png".format(benchmark_names[i]))
    return new_plot

@plot_in_style(range(len(benchmark_names)))
def plot_performance(i, ax):
    print("123")
    plot_data = pd.DataFrame(data=plot_arrays[i],
      columns=['bench','function','flops','cycles','performance','speedup',''])
    print(plot_data)

    plot_data.plot.scatter('bench', 'cycles', s=80,
      title="{0}".format(benchmark_names[i]),
      loglog=False, ax=ax)
    # hide dots for data points
    '''
    bench_df = plot_data
    bench_active = bench_df[bench_df['function'].str.contains('active') | (bench_df['function'].str.contains('fast') & ~bench_df['function'].str.contains('base')) | bench_df['function'].str.contains('dag')]
    bench_base = bench_df[bench_df['function'].str.contains('base')]

    bench_df = bench_df[(~bench_df.index.isin(bench_base.index)) & (~bench_df.index.isin(bench_active.index))]
    bench_df.plot.scatter('bench', 'cycles', s=80,
                     title="{0}".format(benchmark_names[i]),
                     loglog=False, ax=ax)

    sns.lineplot(x="bench", y="cycles", data=bench_base)
    sns.lineplot(x="bench", y="cycles", data=bench_active)'''

    sns.catplot(x="bench", y="cycles", hue="function", kind="bar", data=plot_data)
    plt.show()

plt.tight_layout()
plt.show()

print(optim)
plot_data = pd.DataFrame(data=optim,
            columns=['bench','function','flops','cycles','performance','speedup',''])
print(plot_data)
sns.catplot(x="bench", y="speedup", hue="function", kind="bar", data=plot_data)
plt.show()

# Load an example dataset with long-form data
fmri = sns.load_dataset("fmri")

# Plot the responses for different events and regions
'''g = sns.lineplot(x="timepoint", y="signal",
             hue="region", style="event",
             data=fmri)
g.figure.savefig("testing.png")'''

plt.show()