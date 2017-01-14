from argparse import ArgumentParser

def parse_channels(s):
    splt = (':'+s).strip().split(':')
    chs = [ch for ch in splt[-1].strip().split(',')]
    title = splt[-2].strip()
    return {'title': title, 'chs': chs}

parser = ArgumentParser()
parser.add_argument('file', help='Resampled power NetCDF file.')
parser.add_argument(
    '-c', '--chs', dest='spec', nargs='+',
    default=None, type=parse_channels,
    help='Title and channels for each plot, title:ch1,ch2,ch3.',
)
parser.add_argument(
    '--minor', action='store_true', default=False,
    help='Include minor axis grid.',
)
args = parser.parse_args()


import os
import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['image.cmap'] = 'viridis'

import vstools

fpath = os.path.normpath(args.file)

pwr = xr.open_dataset(fpath, engine='h5netcdf')

if args.spec is None:
    spec = [dict(title='', chs=pwr.data_vars.keys())]
else:
    spec = args.spec
plotvars = [d['chs'] for d in spec]
titles = [d['title'] for d in spec]

pt = pwr.load().median(dim='f')
pt_db = pt.apply(lambda x: 10*np.log10(x), keep_attrs=True)

fig, axarr = vstools.plot_dataset(
    pt_db, plotvars=plotvars, titles=titles,
    xlabel='Time (UTC)', ylabel='Power (dB)',
)
for ax in axarr:
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    if args.minor:
        ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.grid(True, which='major', linewidth=1.5)
    ax.grid(True, which='minor', linewidth=0.75)
plt.tight_layout()
plt.show()
