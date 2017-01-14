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
    '--pmin', default='0%',
    help='Lower power limit in dB for spectra. (default: %(default)s)',
)
parser.add_argument(
    '--pmax', default='100%',
    help='Upper power limit in dB for spectra. (default: %(default)s)',
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
if args.pmin is not None and not args.pmin.endswith('%'):
    args.pmin = float(args.pmin)
if args.pmax is not None and not args.pmax.endswith('%'):
    args.pmax = float(args.pmax)

def calc_spectra(drf_read_kwargs):
    x = vstools.drf_read_ds(chs=chs, **drf_read_kwargs)
    if not args.noremovespikes:
        y = x.apply(vstools.remove_spikes_da, keep_attrs=True, pfp=args.pfp)
    else:
        y = x
    pwr = vstools.acs_ds(
        y, window=args.window, nperseg=args.nperseg, noverlap=args.noverlap,
        nfft=args.nfft, minints=1, cross=False,
    )
    pwr = pwr.isel(t=0).drop('t').assign_coords(f=lambda x: x.f/1e6)
    pwr_db = pwr.apply(lambda x: 10*np.log10(x), keep_attrs=True)

    return pwr_db


g = vstools.drf_live_chunk_gen(datadirs, filesperchunk=args.filesperchunk)

print('Beginning processing. Type Ctrl + C to send SIGINT and quit.')

try:
    d = g.next()
    pwr_db = calc_spectra(d)

    if chs is None:
        chs = pwr_db.data_vars.keys()

    plt.ion()
    fig, ax = plt.subplots(figsize=(args.width, args.height))

    ds = pwr_db[chs]
    df = vstools.ds_fftshift(ds).to_dataframe()
    df.plot(ax=ax)
    pct = {}
    if isinstance(args.pmin, str) and re.match(r'[0-9]{1,3}\.?[0-9]*\%', args.pmin):
        pct['bottom'] = float(args.pmin.rstrip('%'))
    if isinstance(args.pmax, str) and re.match(r'[0-9]{1,3}\.?[0-9]*\%', args.pmax):
        pct['top'] = float(args.pmax.rstrip('%'))
    # get desired percentiles of data and assign them to pmin/pmax kwargs
    pminmax = dict(zip(pct.keys(), np.percentile(df.values, pct.values())))
    if pminmax:
        ax.set_ylim(**pminmax)
    ax.set_title(title)
    ax.set_xlabel('Frequency (MHz)')
    ax.set_ylabel('Power (dB)')
    fig.tight_layout()
    plt.show()
    plt.pause(0.01)
    sys.stdout.write('.')
    sys.stdout.flush()

    for d in g:
        pwr_db = calc_spectra(d)

        ds = pwr_db[chs]
        df = vstools.ds_fftshift(ds).to_dataframe()
        for line, ch in zip(ax.lines, chs):
            line.set_ydata(df[ch].values)
        bottom, top = ax.get_ylim()
        # get desired percentiles of data and assign them to pmin/pmax kwargs
        pminmax = dict(zip(pct.keys(), np.percentile(df.values, pct.values())))
        if pminmax.has_key('bottom'):
            bottom = min(bottom, pminmax['bottom'])
        if pminmax.has_key('top'):
            top = max(top, pminmax['top'])
        ax.set_ylim(bottom=bottom, top=top)
        plt.draw()
        plt.pause(0.01)
        sys.stdout.write('.')
        sys.stdout.flush()
except KeyboardInterrupt:
    pass
