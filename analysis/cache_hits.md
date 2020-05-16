## on Leonhard
- currently unavailable
- check status at: https://scicomp.ethz.ch/wiki/System_status

## on Ubuntu
- Intel i5-8400 (Coffee Lake microarchitecture - optimization of Skylake)
- (9M Cache, up to 4.00 GHz)
- valgrind --tool=callgrind --callgrind-out-file=results.txt ./build/fastslam1/FasterSlamExe ./input_data/example_webmap.mat
- valgrind --tool=callgrind --callgrind-out-file=results_cache.txt --simulate-cache=yes ./build/fastslam1/FasterSlamExe ./input_data/example_webmap.mat

## results 
- see results.txt
- see results_cache.txt + output (below)
==23077== I   refs:      4,328,205,179
==23077== I1  misses:          175,961
==23077== LLi misses:            3,111
==23077== I1  miss rate:          0.00%
==23077== LLi miss rate:          0.00%
==23077== 
==23077== D   refs:      1,566,284,649  (970,594,366 rd + 595,690,283 wr)
==23077== D1  misses:        8,860,225  (  7,570,145 rd +   1,290,080 wr)
==23077== LLd misses:           14,127  (      9,389 rd +       4,738 wr)
==23077== D1  miss rate:           0.6% (        0.8%   +         0.2%  )
==23077== LLd miss rate:           0.0% (        0.0%   +         0.0%  )
==23077== 
==23077== LL refs:           9,036,186  (  7,746,106 rd +   1,290,080 wr)
==23077== LL misses:            17,238  (     12,500 rd +       4,738 wr)
==23077== LL miss rate:            0.0% (        0.0%   +         0.0%  )
