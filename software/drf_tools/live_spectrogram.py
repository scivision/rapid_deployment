import os
import re
import json
import multiprocessing
import signal
import traceback
from argparse import ArgumentParser

def parse_channels(s):
    splt = (':'+s).strip().split(':')
    chs = [ch for ch in splt[-1].strip().split(',')]
    title = splt[-2].strip()
    return {'title': title, 'chs': chs}

parser = ArgumentParser()
parser.add_argument(
    'dirs', nargs='+',
    help='Data directory or directories.',
)
parser.add_argument(
    '-f', '--filesperchunk', type=int, default=1,
    help='Number of DigitalRF files to process in each chunk.',
)
parser.add_argument(
    '-w', '--window', default='blackmanharris',
    help='Window function, in form of argument to scipy.signal.get_window. (default: %(default)s)',
)
parser.add_argument(
    '-s', '--nperseg', type=int, default=256,
    help='Length of segment for estimating spectrum. (default: %(default)s)',
)
parser.add_argument(
    '-o', '--noverlap', type=int, default=None,
    help='Number of overlapping samples when estimating spectrum. (default: half of nperseg)',
)
parser.add_argument(
    '-n', '--nfft', type=int, default=None,
    help='Length of fft for calculating spectrum. (default: equal to nperseg)',
)
parser.add_argument(
    '-i', '--minints', type=int, default=20,
    help='Number of integration segments to use per chunk of data. (default: %(default)s)',
)
parser.add_argument(
    '-p', '--pfp', type=float, default=1e-4,
    help='Probability of false positives to use for spike removal.',
)
parser.add_argument(
    '--noremovespikes', action='store_true',
    help='Do not remove spikes when processing. (default: %(default)s)',
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
    '--ncols', default=1, type=int,
    help='Number of columns in output image grid. (default: %(default)s)',
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

args = parser.parse_args()


import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits import axes_grid1
import seaborn as sns
plt.rcParams['image.cmap'] = 'viridis'

import vstools
import plotting
import time_utils

datadirs = [os.path.normpath(d) for d in args.dirs]
if args.spec is None:
    title = ''
    chs = None
else:
    title = args.spec['title']
    chs = args.spec['chs']
if args.vmin is not None and not args.vmin.endswith('%'):
    args.vmin = float(args.vmin)
if args.vmax is not None and not args.vmax.endswith('%'):
    args.vmax = float(args.vmax)

def calc_spectrograms(drf_read_kwargs):
    x = vstools.drf_read_ds(chs=chs, **drf_read_kwargs)
    if not args.noremovespikes:
        y = x.apply(vstools.remove_spikes_da, keep_attrs=True, pfp=args.pfp)
    else:
        y = x
    pwr = vstools.acs_ds(
        y, window=args.window, nperseg=args.nperseg, noverlap=args.noverlap,
        nfft=args.nfft, minints=args.minints, cross=False,
    )
    pwr = pwr.assign_coords(f=lambda x: x.f/1e6)
    pwr_db = pwr.apply(lambda x: 10*np.log10(x), keep_attrs=True)

    return pwr_db


g = vstools.drf_live_chunk_gen(datadirs, filesperchunk=args.filesperchunk)

print('Beginning processing. Type Ctrl + C to send SIGINT and quit.')

try:
    d = g.next()
    pwr_db = calc_spectrograms(d)

    if chs is None:
        chs = pwr_db.data_vars.keys()
    if args.nrows is None:
        args.nrows = int(np.ceil(len(chs)/float(args.ncols)))

    plt.ion()
    fig = plt.figure(figsize=(args.width, args.height))
    gs = mpl.gridspec.GridSpec(args.nrows, args.ncols)
    h = [axes_grid1.Size.Scaled(1.0)]
    v = [axes_grid1.Size.Scaled(1.0)]
    axlist = []
    ax0 = None
    for ss in gs:
        div = axes_grid1.SubplotDivider(fig, ss, horizontal=h, vertical=v)
        loc = div.new_locator(nx=0, ny=0)
        if ax0 is None:
            ax = fig.add_axes(loc(None, None))
            ax0 = ax
        else:
            ax = fig.add_axes(loc(None, None), sharex=ax0, sharey=ax0)
        ax.set_axes_locator(loc)
        ax.get_subplotspec = loc.get_subplotspec
        axlist.append(ax)

    images = []
    for ax, ch in zip(axlist, chs):
        ds = pwr_db[[ch]]
        da = ds.values()[0]
        img = vstools.plot_spectrogram(
            da.transpose('t', 'f').values, ds.t.values, ds.f.values,
            ax=ax, ylabel='Frequency (MHz)', clabel='Power (dB)',
            vmin=args.vmin, vmax=args.vmax,
        )
        ax.set_title(ch)
        images.append(img)
    fig.tight_layout()
    plt.show()
    plt.pause(0.01)
    sys.stdout.write('.')
    sys.stdout.flush()

    for d in g:
        pwr_db = calc_spectrograms(d)
        for ax, img, ch in zip(axlist, images, chs):
            ds = pwr_db[[ch]]
            da = ds.values()[0]
            x = da.transpose('t', 'f').values
            img.set_data(np.fft.fftshift(x, axes=1).T)
            # extract vmin/vmax string percentage arguments into floats in dict
            # indexed by 'vmin' or 'vmax', respectively
            pct = {}
            if isinstance(args.vmin, str) and re.match(r'[0-9]{1,3}\.?[0-9]*\%', args.vmin):
                pct['vmin'] = float(args.vmin.rstrip('%'))
            if isinstance(args.vmax, str) and re.match(r'[0-9]{1,3}\.?[0-9]*\%', args.vmax):
                pct['vmax'] = float(args.vmax.rstrip('%'))
            # get desired percentiles of data and assign them to vmin/vmax kwargs
            vminmax = dict(zip(pct.keys(), np.percentile(x, pct.values())))
            if vminmax:
                img.set_clim(**vminmax)
            ax.xaxis.set_label_text('')
            ts = ds.t.values[0]
            te = ds.t.values[-1]
            epoch = ts.astype('datetime64[D]').astype(ts.dtype)
            tsf = time_utils.datetime_to_float(ts, epoch)
            tef = time_utils.datetime_to_float(te, epoch)
            step = (tef - tsf)/(ds.t.values.shape[0] - 1)
            old_ext = img.get_extent()
            img.set_extent((tsf - step/2.0, tef + step/2.0, old_ext[2], old_ext[3]))
            plotting.timeticks(ax.xaxis, ts, te, epoch)
        plt.draw()
        plt.pause(0.01)
        sys.stdout.write('.')
        sys.stdout.flush()
except KeyboardInterrupt:
    pass
