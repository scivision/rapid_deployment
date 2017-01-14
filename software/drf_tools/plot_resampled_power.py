from argparse import ArgumentParser

def parse_channels(s):
    splt = (':'+s).strip().split(':')
    chs = [ch for ch in splt[-1].strip().split(',')]
    title = splt[-2].strip()
    return {'title': title, 'chs': chs}

parser = ArgumentParser()
parser.add_argument('file', help='Resampled power NetCDF file.')
parser.add_argument(
    '-c', '--chs', dest='spec',
    default=None, type=parse_channels,
    help='Title and channels to plot, title:ch1,ch2,ch3.',
)
parser.add_argument(
    '--pmin', default=None, type=float,
    help='Lower power limit in dB for time series plot.',
)
parser.add_argument(
    '--pmax', default=None, type=float,
    help='Upper power limit in dB for time series plot.',
)
parser.add_argument(
    '--vmin', default='2%',
    help='Lower power limit in dB for spectrograms. (default: %(default)s)',
)
parser.add_argument(
    '--vmax', default='98%',
    help='Upper power limit in dB for spectrograms. (default: %(default)s)',
)
parser.add_argument(
    '--samescale', action='store_true', default=False,
    help='Set vmin/vmax to match time series power scale.',
)
parser.add_argument(
    '--minor', action='store_true', default=False,
    help='Include minor axis grid.',
)
args = parser.parse_args()


import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['image.cmap'] = 'viridis'

import vstools

fpath = os.path.normpath(args.file)

pwr = xr.open_dataset(fpath, engine='h5netcdf')

if args.spec is None:
    spec = dict(title='', chs=pwr.data_vars.keys())
else:
    spec = args.spec

pwr = pwr.assign_coords(f=lambda x: x.f/1e6)

kwargs = {}
if args.vmin is not None:
    if not args.vmin.endswith('%'):
        args.vmin = float(args.vmin)
    kwargs['vmin'] = args.vmin
if args.vmax is not None:
    if not args.vmax.endswith('%'):
        args.vmax = float(args.vmax)
    kwargs['vmax'] = args.vmax

axlist = vstools.spectrograms_and_timeseries(
    pwr[spec['chs']], ax=None,
    ptransform=lambda p: 10*np.log10(p),
    flabel='Frequency (MHz)', plabel='Power (dB)',
    pmin=args.pmin, pmax=args.pmax, set_vminmax=args.samescale,
    minor=args.minor, **kwargs
)
axlist[0].set_title(spec['title'])
axlist[-1].set_xlabel(axlist[-1].get_xlabel() + ', UTC')
plt.tight_layout()
plt.show()
