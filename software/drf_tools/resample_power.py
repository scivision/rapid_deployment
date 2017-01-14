from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('dir', help='Processed power data directory.')
parser.add_argument(
    'file', nargs='?', default=None,
    help="""Resampled power NetCDF save file.
            (default: {dir}/pwr_resampled@{t0}--{t1}.nc)""",
)
parser.add_argument(
    '-f', '--freq', default='30s',
    help='Time resample period. (default: %(default)s)',
)
parser.add_argument(
    '-m', '--method', default='mean',
    help='Time resample method. (default: %(default)s)',
)
parser.add_argument(
    '-t', '--index', default=None, metavar='T0[--T1]',
    help="""Time index to select before resampling.
            Pass "?" to print start and end of data. (default: %(default)s)""",
)
args = parser.parse_args()


import os
import sys
import numpy as np
import xarray as xr

datadir = os.path.normpath(args.dir)

pwr = xr.open_mfdataset(os.path.join(datadir, 'pwr@*.nc'), concat_dim='t', engine='h5netcdf')

if args.index == '?':
    tidx = pwr.indexes['t']
    print('Data goes from {0} to {1}.'.format(tidx[0], tidx[-1]))
    sys.exit()
elif args.index is not None:
    sections = args.index.split('--')
    if len(sections) == 1:
        index = sections[0].strip()
    elif len(sections) == 2:
        start = sections[0].strip()
        stop = sections[1].strip()
        index = slice(start, stop)
    else:
        raise ValueError('Index must be of the form "x" or "x to y"')

    pwr = pwr.sel(t=index)

pwr_rs = pwr.resample(args.freq, dim='t', how=args.method)

if args.file is None:
    tidx = pwr.indexes['t']
    t0, t1 = [t.strftime('%Y%m%d%H%M%S') for t in tidx[0], tidx[-1]]
    fname = 'pwr_resampled@{0}--{1}.nc'.format(t0, t1)
    fname = fname.replace(' ', '_')
    fpath = os.path.join(datadir, fname)
else:
    fpath = os.path.normpath(args.file)

pwr_rs.to_netcdf(fpath, engine='h5netcdf')
