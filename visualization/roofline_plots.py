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
parser = argparse.ArgumentParser(description='Generate FasterSlam roofline plots from benchmark data')
parser.add_argument('-i', '--input', required=True, help="Path to benchmark.csv file")
parser.add_argument('-o', '--output', help="Output directory")

args = parser.parse_args()
print(args.input)
infile_path = Path(args.input)
assert(infile_path.exists())
assert(infile_path.is_file())
assert(infile_path.parts[-1].endswith('.csv'))

output_path = Path(args.output)
assert(output_path.exists())
assert(output_path.is_dir())

## SETUP DATAFRAME
df = pd.read_csv(infile_path, sep=';',
                           names=['benchmark', 'version', 'size', 'work', 'memory', 'time', 'performance', 'speedup'])
df.fillna(1, axis='columns', inplace=True)
df['op_intensity'] = df['work'] / df['memory']
benchmarks = df['benchmark'].unique()
num_benches = benchmarks.size
print(df.head())

def label_line(ax, line, label):
    """Add a label to a line, at the proper angle.

    Arguments
    ---------
    line : matplotlib.lines.Line2D object,
    label : str
    x : float
        x-position to place center of text (in data coordinated
    y : float
        y-position to place center of text (in data coordinates)
    color : str
    size : float
    """
    xdata, ydata = line.get_data()
    x1 = xdata[0]
    x2 = xdata[-1]
    y1 = ydata[0]
    y2 = ydata[-1]

    text = ax.annotate(label, xy=(0.1, 0.15), xycoords=line,
                       xytext=(-10, 0),
                       textcoords='offset pixels',
                       horizontalalignment='left',
                       verticalalignment='bottom',
                       fontstyle='italic')

    sp1 = ax.transData.transform_point((x1, y1))
    sp2 = ax.transData.transform_point((x2, y2))

    rise = (sp2[1] - sp1[1])
    run = (sp2[0] - sp1[0])

    slope_degrees = np.degrees(np.arctan2(rise, run))
    trans_angle = ax.transData.transform_angles(np.array((slope_degrees,)),
                                                   np.reshape(text.get_position(), (1,2)))[0]
    text.set_rotation(slope_degrees)
    return text

def get_memory_bounds_points(ax):
    x_min = min(min(l.get_xdata()) for l in ax.get_lines())  # will be from hline
    x_max = max(max(l.get_xdata()) for l in ax.get_lines())  # will be from hline
    # x_min = 1/10 * x_min - (0.1*x_max+0.0*x_min)
    x_min = 1 * ax.get_xlim()[0]
    y_max = 2*peak_perf_avx

    ret_min = (x_min, peak_memory * x_min)
    ret_max = (1/peak_memory * y_max, y_max)
    return (ret_min, ret_max)

def setup_axis(ax):
    ax.xaxis.set_label_coords(1, -0.075)
    ax.yaxis.set_label_coords(0, 1.005)
    ax.set_ylabel(ax.get_ylabel(), rotation='horizontal')
    ax.labelweight = 'bold'
    l1 = ax.axhline(peak_perf_seq, linestyle='--', c='m')
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
                    fontstyle='italic')
    (ret_min, ret_max) = get_memory_bounds_points(ax)
    l3, = ax.plot([ret_min[0], ret_max[0]], [ret_min[1], ret_max[1]], '--', c='r')
    #plt.annotate(f"Read/write limit\n({peak_memory:2.1f} bytes/cycle)",
    #                (0.0, 0.), 
    #                rotation=atan2((ret_max[1]-ret_min[1]),(ret_max[0]-ret_min[0])) * 180/pi,
    #                xycoords=l3, 
    #                horizontalalignment='left',
    #                verticalalignment='bottom',
    #                fontstyle='italic')
    ax.set_xlabel("operational intensity\n[ops/byte]", horizontalalignment='right')
    ax.set_ylabel("performance\n[ops/cycle]", rotation='horizontal', horizontalalignment='left')
    ax.get_xaxis().set_major_formatter(StrMethodFormatter('{x:3.2f}'))
    ax.get_yaxis().set_major_formatter(StrMethodFormatter('{x:2.2f}'))
    ax.tick_params('x', direction='in')
    ax.legend(loc='lower left')
    # ax.set_aspect('equal', adjustable='box')
    return l3

## DECORATORS
def plot_in_style(elements):
    def new_plot(f):
        for i, el in enumerate(elements):
            fig, ax = plt.subplots(figsize=(5,4))
            f(el, ax)
            memory_line = setup_axis(ax)
            fig.tight_layout()
            label_line(ax, memory_line, f"Read/write limit\n({peak_memory:2.1f} bytes/cycle)")
            fig.savefig(f"{output_path}/{benchmarks[i]}.png")
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
#@plot_in_grid(range(4))
#def plot_roofline(i, ax):
    #df[df['benchmark'] == benchmarks[i]].plot.scatter('op_intensity', 'performance', loglog=True, ax=ax)
#@plot_in_grid([14])
#def plot_roofline(i, ax):
    #df[df['benchmark'] == benchmarks[i]]. \
        #plot.scatter('op_intensity', 'performance', 
                     #title=benchmarks[i].split(' ')[0],
                     #loglog=True, ax=ax)

#@plot_in_style([i for i in range(num_benches) if 'KF_' in benchmarks[i]])
#def plot_roofline(i, ax):
#    bench_df = df[df['benchmark'] == benchmarks[i]]
#    bench_base = bench_df[bench_df['version'].str.contains('base')]
#    bench_active = bench_df[bench_df['version'].str.contains('active')]
#    bench_df = bench_df[(~bench_df.index.isin(bench_base.index)) & (~bench_df.index.isin(bench_active.index))]
#    bench_df.plot.scatter('op_intensity', 'performance', s=40,
#                     title=benchmarks[i].split(' ')[0],
#                     loglog=True, ax=ax)
#    bench_active.plot.scatter('op_intensity', 'performance', marker='^', c='g', s=40,
#                          loglog=True, ax=ax, label='active')
#    bench_base.plot.scatter('op_intensity', 'performance', marker='v', c='c', s=40,
#                          loglog=True, ax=ax, label='base')

@plot_in_style(range(num_benches))
def plot_roofline(i, ax):
    bench_df = df[df['benchmark'] == benchmarks[i]]
    bench_base = bench_df[bench_df['version'].str.contains('base')]
    bench_active = bench_df[bench_df['version'].str.contains('active')]
    bench_df = bench_df[(~bench_df.index.isin(bench_base.index)) & (~bench_df.index.isin(bench_active.index))]
    bench_df.plot.scatter('op_intensity', 'performance', s=80,
                     title=benchmarks[i].split(' ')[0],
                     loglog=True, ax=ax)
    bench_active.plot.scatter('op_intensity', 'performance', marker='^', c='g', s=80,
                          loglog=True, ax=ax, label='active')
    bench_base.plot.scatter('op_intensity', 'performance', marker='v', c='c', s=80,
                          loglog=True, ax=ax, label='base')

plt.tight_layout()
#plt.show()