from compare_polang import *
import matplotlib as mpl

# This is needed if you want to write out lots of plots
mpl.rcParams['agg.path.chunksize'] = 10000

# compare_polang(prefix='mfi',date='201905')
#compare_polang(prefix='mfi302',date='201904',use_variance=False)
#compare_polang(prefix='mfi',date='201907')
# compare_polang(prefix='mfi',date='201911',datestr='nov2019')
# compare_polang(prefix='mfi',date='201911',datestr='nov2019_half2')

# compare_polang(prefix='mfi',date='202003',datestr='mar2020')
compare_polang(prefix='mfi',date='202004',datestr='apr2020')
