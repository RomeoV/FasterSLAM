import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

total = 1900
SELF = np.array([
    441,
    279,
    108,
    102,
    100,
    99,
    84,
    83,
    75,
    74,
    71,
    65], dtype=np.float64)/total

CALLED = np.array([
    4.1,
    1.7,
    3.3,
    3.6,
    1.7,
    1.7,
    3.6,
    2.7,
    1.4,
    0.2,
    1.75,
    1.76])

NAME = np.array([
    'mul',
    'sincos',
    'copy',
    'random_r',
    'predict_base',
    '???',
    'random',
    'add',
    'transpose',
    'compute_jacobians',
    'multivariate_gauss',
    ''])

df = pd.DataFrame({'time':SELF, 'called':CALLED, 'name':NAME})
ax = df.plot.bar('name', 'time', width = 0.9, legend=False, alpha=0.6, color='royalblue')
ax.set_ylabel('percentage of total runtime')
ax2 = ax.twinx()
df.plot.bar('name', 'called', ax=ax2, legend=False, alpha=1.0, color='darkorange')
ax2.set_ylabel('Nb of calls in millions')
ax.figure.legend()
plt.title('Most costly functions (rel. runtime and # of calls)')
plt.show()
