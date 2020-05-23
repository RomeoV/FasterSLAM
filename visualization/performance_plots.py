# python performance_plots.py -i polybox/fastslam1_feature.csv -t features
#  python performance_plots.py -i polybox/fastslam1_particle.csv -t particles

# import warnings; warnings.simplefilter('error', UserWarning)
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import matplotlib
from math import floor, sqrt, ceil, log, exp, atan2, pi
from functools import wraps
import seaborn as sns
import numpy as np

## SETUP MATPLOTLIB
plt.style.use('seaborn')
matplotlib.rcParams['axes.labelweight'] = 'bold'
matplotlib.rcParams['axes.titleweight'] = 'bold'
matplotlib.rcParams['font.weight'] = 'semibold'
sns.set_color_codes()

## SETUP CONSTANTS
peak_perf_seq = 4.0
peak_perf_avx = 16.0
peak_memory = 6.0

## PARSE INPUTS
parser = argparse.ArgumentParser(description='Generate FasterSlam performance plots from benchmark data')
parser.add_argument('-i', '--input', required=True, help="Path to benchmark.csv file")
parser.add_argument('-t', '--type', required=True, help="Choose features/particles")

args = parser.parse_args()
print(args.input)
infile_path = Path(args.input)
assert(infile_path.exists())
assert(infile_path.is_file())
assert(infile_path.parts[-1].endswith('.csv'))

## SETUP DATAFRAME
df = pd.read_csv(infile_path, sep=';',
                           names=['benchmark', 'version', 'size', 'work', 'memory', 'time', 'performance', 'speedup'])
df.fillna(1, axis='columns', inplace=True)
versions = df['version'].unique()
num_versions = versions.size
print(num_versions)
print(df.head())

def setup_axis(ax):
    ax.xaxis.set_label_coords(1, -0.075)
    ax.yaxis.set_label_coords(0, 1.005)
    ax.set_ylabel(ax.get_ylabel(), rotation='horizontal')
    ax.labelweight = 'bold'
    '''l1 = ax.axhline(peak_perf_seq, linestyle='--', c='m')
    plt.annotate(f"Peak perf. seq\n({peak_perf_seq:2.1f} flops/cycle)",
                    (1, 0.), 
                    xycoords=l1, 
                    horizontalalignment='right',
                    verticalalignment='top',
                    fontstyle='italic')
    l2 = ax.axhline(peak_perf_avx, linestyle='--', c='m')
    plt.annotate(f"Peak perf. AVX\n({peak_perf_avx:2.1f} flops/cycle)",
                    (1, 0.), 
                    xycoords=l2, 
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    fontstyle='italic')'''
    ax.set_xlabel("{0}\n".format(args.type), horizontalalignment='right')
    ax.set_ylabel("performance\n[ops/cycle]", rotation='horizontal', horizontalalignment='left')
    ax.get_xaxis().set_major_formatter(StrMethodFormatter('{x:3.2f}'))
    ax.get_yaxis().set_major_formatter(StrMethodFormatter('{x:2.2f}'))
    ax.tick_params('x', direction='in')
    ax.legend(loc='lower left')
    # ax.set_aspect('equal', adjustable='box')

## DECORATORS
def plot_in_style(elements):
    def new_plot(f):
        for i, el in enumerate(elements):
            fig, ax = plt.subplots(figsize=(5,4))
            f(el, ax)
            setup_axis(ax)
            fig.tight_layout()
            fig.savefig(f"./{0}_performance.png".format(args.type))
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

@plot_in_style(range(1))
def plot_performance(i, ax):
    bench_df = df
    # hide dots for data points
    bench_active = bench_df[bench_df['version'].str.contains('active')]
    bench_base = bench_df[bench_df['version'].str.contains('base')]

    bench_df = bench_df[(~bench_df.index.isin(bench_base.index)) & (~bench_df.index.isin(bench_active.index))]
    bench_df.plot.scatter('size', 'performance', s=80,
                     title="Performance/{0}".format(args.type),
                     loglog=False, ax=ax)

    sns.lineplot(x="size", y="performance", data=bench_base)
    sns.lineplot(x="size", y="performance", data=bench_active)

plt.tight_layout()
plt.show()