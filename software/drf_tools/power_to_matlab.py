from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument(
    'file', nargs='+',
    help='File or files containing processed power.'
)
parser.add_argument(
    '-n', '--groupsize', type=int, default=1,
    help='Number of consecutive sorted files to group into the same Matlab file. (default: %(default)s)',
)
args = parser.parse_args()


import glob
import os
import xarray as xr
import numpy as np
from scipy.io import savemat

import vstools

pwrfiles = sorted(args.file)
for k, idx in enumerate(range(0, len(pwrfiles), args.groupsize)):
    print('Saving chunk {0}/{1}'.format(k+1, (len(pwrfiles) - 1)//args.groupsize + 1))
    pfiles = pwrfiles[idx:idx+args.groupsize]
    ds = xr.open_mfdataset(pfiles, engine='h5netcdf')
    S = vstools.ds_acs_to_matrix(ds)
    chs = ds.autovars.astype(np.object_)
    f = ds.f.values
    t = ds.t.values
    t0 = t[0].astype(np.int64)/1e9
    dt = (t - t[0])/np.timedelta64(1, 's')

    dirname, fname = os.path.split(os.path.splitext(pfiles[0])[0])
    matpath = os.path.join(dirname, fname + '.mat')
    savemat(matpath, dict(S=S, chs=chs, f=f, t0=t0, dt=dt), oned_as='column')
