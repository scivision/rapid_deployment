from argparse import ArgumentParser

def parse_channels(s):
    splt = (':'+s).strip().split(':')
    chs = [ch for ch in splt[-1].strip().split(',')]
    title = splt[-2].strip()
    return {'title': title, 'chs': chs}

parser = ArgumentParser()
parser.add_argument(
    'file', nargs='+',
    help='File or files containing processed power.'
)
parser.add_argument(
    '-n', '--groupsize', type=int, default=1,
    help='Number of consecutive sorted files to group into the same plot. (default: %(default)s)',
)
parser.add_argument(
    '-c', '--chs', dest='spec',
    default=None, type=parse_channels,
    help='Title and channels to plot, title:ch1,ch2,ch3.',
)
parser.add_argument(
    '--nrows', default=None, type=int,
    help='Number of rows in output image grid. (default: as many as necessary)',
)
parser.add_argument(
    '--ncols', default=None, type=int,
    help='Number of columns in output image grid. (default: 1)',
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
    '--width', default=14, type=float,
    help='Figure width in inches. (default: %(default)s)',
)
parser.add_argument(
    '--height', default=8, type=float,
    help='Figure height in inches. (default: %(default)s)',
)
parser.add_argument(
    '--ext', default='png',
    help='Extension and file type of output image. (default: %(default)s)',
)
args = parser.parse_args()


import glob
import os
import numpy as np
import xarray as xr
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['image.cmap'] = 'viridis'

import vstools

kwargs = dict(ylabel='Frequency (MHz)', clabel='Power (dB)')
if args.nrows is not None:
    kwargs['nrows'] = args.nrows
if args.ncols is not None:
    kwargs['ncols'] = args.ncols
if args.vmin is not None:
    if not args.vmin.endswith('%'):
        args.vmin = float(args.vmin)
    kwargs['vmin'] = args.vmin
if args.vmax is not None:
    if not args.vmax.endswith('%'):
        args.vmax = float(args.vmax)
    kwargs['vmax'] = args.vmax

pwrfiles = sorted(args.file)

for k, idx in enumerate(range(0, len(pwrfiles), args.groupsize)):
    print('Plotting chunk {0}/{1}'.format(k+1, (len(pwrfiles) - 1)//args.groupsize + 1))
    pfiles = pwrfiles[idx:idx+args.groupsize]
    with xr.open_mfdataset(pfiles, engine='h5netcdf') as pwr:
        if args.spec is None:
            spec = dict(title='', chs=pwr.data_vars.keys())
        else:
            spec = args.spec

        pwr = pwr[spec['chs']]
        pwr = pwr.assign_coords(f=lambda x: x.f/1e6)
        pwr_db = pwr.apply(lambda x: 10*np.log10(x), keep_attrs=True)

        fig, axlist = vstools.spectrograms(
            pwr_db, title=spec['title'],
            fig_kw=dict(figsize=(args.width, args.height)),
            **kwargs
        )
        fig.tight_layout()

        dirname, fname = os.path.split(os.path.splitext(pfiles[0])[0])
        figdir = os.path.join(dirname, 'spectrograms')
        if not os.path.exists(figdir):
            os.makedirs(figdir)
        figname = '{c}_{f}.{e}'.format(c=','.join(spec['chs']), f=fname, e=args.ext)
        figpath = os.path.join(figdir, figname)
        fig.savefig(figpath, bbox_inches='tight', pad_inches=0.05)

        plt.close(fig)
