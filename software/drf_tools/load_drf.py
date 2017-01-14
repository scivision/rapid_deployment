import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['image.cmap'] = 'viridis'
import os
from argparse import ArgumentParser
from scipy.stats import chi2

import digital_rf_hdf5 as drf
import vstools

parser = ArgumentParser()
parser.add_argument('dirs', nargs='+', help='Data directory or directories.')

args = parser.parse_args()

datadirs = [os.path.normpath(d) for d in args.dirs]

chunkgen = vstools.drf_chunk_gen(datadirs, filesperchunk=10)

def vltgen():
    for d in chunkgen:
        yield vstools.drf_read_ds(**d)

g = vltgen()
vlt = g.next()

def plot_power(vlt, **kwargs):
    pwr = vlt.apply(lambda x: 20*np.log10(np.abs(x)), keep_attrs=True)
    t = pwr.t.data.astype(np.int_)
    tus = (t-t[0])/1e6
    pwr = pwr.assign_coords(t=('t', tus))

    kw = dict(plotvars=[pwr.data_vars.keys()], xlabel='Time (us)', ylabel='Power (dB)')
    kw.update(kwargs)
    fig, axlist = vstools.plot_dataset(pwr, **kw)

    plt.tight_layout()

    return fig, axlist

def plot_spectrograms(vlt, window='blackmanharris', nperseg=1024, noverlap=None, nfft=None, minints=1000, **kwargs):
    acs = vstools.acs_ds(vlt, window=window, nperseg=nperseg, noverlap=noverlap, nfft=nfft, minints=minints, cross=False)
    spec = 10*acs.apply(np.log10, keep_attrs=True)
    spec = spec.assign_coords(f=spec.f/1e6)

    kw = dict(ylabel='Frequency (MHz)', clabel='Power (dB)')
    kw.update(kwargs)
    fig, axlist = vstools.spectrograms(spec, **kw)

    plt.tight_layout()

    return fig, axlist

def plot_spectra(vlt, window='blackmanharris', nperseg=1024, noverlap=None, nfft=None, minints=1000, **kwargs):
    acs = vstools.acs_ds(vlt, window=window, nperseg=nperseg, noverlap=noverlap, nfft=nfft, minints=minints, cross=False)
    pf = acs.mean(dim='t', keep_attrs=True)
    pf_db = 10*pf.apply(np.log10, keep_attrs=True)
    pf_db = pf_db.assign_coords(f=pf_db.f/1e6)

    kw = dict(plotvars=[pf_db.data_vars.keys()], xlabel='Frequency (MHz)', ylabel='Power (dB)')
    kw.update(kwargs)
    fig, axlist = vstools.plot_dataset(vstools.ds_fftshift(pf_db), **kw)

    plt.tight_layout()

    return fig, axlist
