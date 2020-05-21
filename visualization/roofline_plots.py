import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import matplotlib
from math import floor, sqrt, ceil
from functools import wraps
import seaborn as sns

## SETUP MATPLOTLIB
plt.style.use('seaborn')
matplotlib.rcParams['axes.labelweight'] = 'bold'
matplotlib.rcParams['axes.titleweight'] = 'bold'
matplotlib.rcParams['font.weight'] = 'semibold'
matplotlib.rcParams['lines.markersize'] = 12
sns.set_color_codes()

## SETUP CONSTANTS
peak_perf_seq = 2.0
peak_perf_avx = 8.0
peak_memory = 6.0

## PARSE INPUTS
parser = argparse.ArgumentParser(description='Generate FasterSlam roofline plots from benchmark data')
parser.add_argument('-i', '--input')

args = parser.parse_args()
print(args.input)
infile_path = Path(args.input)
assert(infile_path.exists())
assert(infile_path.is_file())
assert(infile_path.parts[-1].endswith('.csv'))

## SETUP DATAFRAME
df = pd.read_csv(infile_path, sep=';',
                           names=['benchmark', 'version', 'work', 'memory', 'time', 'performance', 'speedup', '_'])
df['memory'] = 1
df['op_intensity'] = df['work'] / df['memory']
benchmarks = df['benchmark'].unique()
num_benches = benchmarks.size
print(df.head())

def get_memory_bounds_points(ax):
    x_min = min(min(l.get_xdata()) for l in ax.lines)
    x_max = max(max(l.get_xdata()) for l in ax.lines)
    x_min = x_min - (0.5*(x_max - x_min) + 0.1*x_min)
    y_max = 2*peak_perf_avx

    ret_min = (x_min, peak_memory * x_min)
    ret_max = (1/peak_memory * y_max, y_max)
    return (ret_min, ret_max)

def setup_axis(ax):
    ax.xaxis.set_label_coords(1, -0.005)
    ax.yaxis.set_label_coords(0, 1.005)
    ax.set_aspect(1/5)
    ax.set_ylabel(ax.get_ylabel(), rotation='horizontal')
    ax.labelweight = 'bold'
    l1 = ax.axhline(peak_perf_seq, linestyle='--', c='m')
    plt.annotate(f"Peak perf. seq\n({peak_perf_seq:2.1f} flops/cycle)",
                    (1, 0.), 
                    xycoords=l1, 
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    fontstyle='italic')
    l2 = ax.axhline(peak_perf_avx, linestyle='--', c='m')
    plt.annotate(f"Peak perf. AVX\n({peak_perf_avx:2.1f} flops/cycle)",
                    (1, 0.), 
                    xycoords=l2, 
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    fontstyle='italic')
    (ret_min, ret_max) = get_memory_bounds_points(ax)
    # ax.plot([ret_min[0], ret_max[0]], [ret_min[1], ret_max[1]], c='r')
    ax.set_xlabel("operational intensity\n[ops/byte]", horizontalalignment='right')
    ax.set_ylabel("performance\n[ops/cycle]", rotation='horizontal', horizontalalignment='left')
    ax.get_xaxis().set_major_formatter(StrMethodFormatter('{x:2.0f}'))
    ax.get_yaxis().set_major_formatter(StrMethodFormatter('{x:2.2f}'))
    ax.legend(loc='lower left')
    return ax

## DECORATORS
def plot_in_style(elements):
    def new_plot(f):
        for i, el in enumerate(elements):
            fig, ax = plt.subplots()
            f(el, ax)
            setup_axis(ax)
    return new_plot

def plot_in_grid(elements):
    def repeat_plot(f):
        num_elements = len(elements)
        M = int(floor(sqrt(num_elements)))
        N = int(ceil(num_elements/M))
        assert(M*N >= num_elements)
        plt.figure()
        for i, el in enumerate(elements):
            ax = plt.subplot(M, N, i+1)
            f(el, ax)
            setup_axis(ax)
    return repeat_plot

## PLOT DATA
@plot_in_grid(range(4))
def plot_roofline(i, ax):
    df[df['benchmark'] == benchmarks[i]].plot.scatter('op_intensity', 'performance', loglog=True, ax=ax)
@plot_in_grid([14])
def plot_roofline(i, ax):
    df[df['benchmark'] == benchmarks[i]]. \
        plot.scatter('op_intensity', 'performance', 
                     title=benchmarks[i].split(' ')[0],
                     loglog=True, ax=ax)

@plot_in_style([i for i in range(num_benches) if 'KF_' in benchmarks[i]])
def plot_roofline(i, ax):
    bench_df = df[df['benchmark'] == benchmarks[i]]
    bench_base = bench_df[bench_df['version'].str.contains('base')]
    bench_active = bench_df[bench_df['version'].str.contains('active')]
    bench_df = bench_df[(~bench_df.index.isin(bench_base.index)) & (~bench_df.index.isin(bench_active.index))]
    bench_df.plot.scatter('op_intensity', 'performance', s=40,
                     title=benchmarks[i].split(' ')[0],
                     loglog=True, ax=ax)
    bench_active.plot.scatter('op_intensity', 'performance', marker='^', c='g', s=40,
                          loglog=True, ax=ax, label='active')
    bench_base.plot.scatter('op_intensity', 'performance', marker='v', c='c', s=40,
                          loglog=True, ax=ax, label='base')

plt.tight_layout()
plt.show()