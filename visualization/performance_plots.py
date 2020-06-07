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
parser.add_argument('-i', '--input', required=False, help="Path to benchmark.csv file")
parser.add_argument('-t', '--title', required=False, help="Choose features/particles")

args = parser.parse_args()

if(args.input):
    print(args.input)
    infile_path = Path(args.input)
    assert(infile_path.exists())
    assert(infile_path.is_file())
    assert(infile_path.parts[-1].endswith('.csv'))
    file_names = [args.input]
    file_titles = [args.title]
else:
    file_names = "polybox/fastslam1_particle.csv", "polybox/fastslam1_feature.csv", "polybox/resample_particles_scale_particles.csv" 
    #file_names = ["polybox/observe_update_scale_particles.csv", "polybox/predict_update_scale_particles.csv"]
    file_titles = ["Fastslam1 Scaling with Particles", "Fastslam1 Scaling with Features", "resample_particles  Scaling with Particles"]
    fun_names = ["observe_update", "predict_update"]
## SETUP DATAFRAME
dfs = []

for i in range(len(file_names)):
    df = pd.read_csv(file_names[i], sep=';',
                            names=['benchmark', 'version', 'size', 'work', 'memory', 'time', 'performance', 'speedup'], \
                            dtype={'benchmark': str, 'version':str, 'size':int, 'work':int, 'memory':int, 'time':float, 'performance':float, 'speedup':float})
    df.fillna(1, axis='columns', inplace=True)
    dfs.append(df)
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
    ax.set_xlabel("input size\n", horizontalalignment='right')
    ax.set_ylabel("performance\n[ops/cycle]", rotation='horizontal', horizontalalignment='left')
    # ax.get_xaxis().set_major_formatter(StrMethodFormatter('{x:3.2f}'))
    ax.get_yaxis().set_major_formatter(StrMethodFormatter('{x:2.2f}'))
    ax.tick_params('x', direction='in')
    ax.set_ylim([0,6])
    #ax.set_xlim([0,14000])
    ax.legend(loc='lower left')
    # ax.set_aspect('equal', adjustable='box')

## DECORATORS
def plot_in_style(elements):
    def new_plot(f):
        
        for i, el in enumerate(elements):
            fig, ax = plt.subplots(figsize=(10,6))
            # fig, ax = plt.subplots(figsize=(10,6))
            f(el, ax)
            setup_axis(ax)
            fig.tight_layout()
            fig.savefig(f"./{0}_performance.png".format(i))
    return new_plot

@plot_in_style(range(len(file_names)))
def plot_performance(i, ax):
    bench_df = dfs[i]
    # hide dots for data points
    bench_active = bench_df[bench_df['version'].str.contains('active') | (bench_df['version'].str.contains('fast') & ~bench_df['version'].str.contains('base')) | bench_df['version'].str.contains('dag')]
    bench_base = bench_df[bench_df['version'].str.contains('base')]

    #bench_df = bench_df[(~bench_df.index.isin(bench_base.index)) & (~bench_df.index.isin(bench_active.index))]
    # bench_df.plot.scatter('size', 'performance', s=80,
    #                  title="{0}".format(file_titles[i]),
    #                  loglog=False, ax=ax)
    sns.lineplot(x="size", y="performance", marker="o",data=bench_base)
    sns.lineplot(x="size", y="performance", marker="o", data=bench_active)
    ax.set_title(file_titles[i])
    #ax.set_yscale("log")
    #sns.lineplot(x="size", y="performance", color=sns.color_palette("Paired")[2*i], marker="o", label=fun_names[i]+" base", data=bench_base)
    #sns.lineplot(x="size", y="performance", color=sns.color_palette("Paired")[2*i+1], marker="o", label=fun_names[i]+" fast", data=bench_active)
    #sns.lineplot(x="size", y="performance", hue="version",data=bench_df)
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.tight_layout()

plt.show()