import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['image.cmap'] = 'viridis'
import os
from argparse import ArgumentParser

from vstools import plot_dataset, plot_spectrogram, spectrograms, ds_fftshift

parser = ArgumentParser()
parser.add_argument('dir', help='Processed power data directory.')

args = parser.parse_args()

datadir = os.path.normpath(args.dir)

pwrs = xr.open_mfdataset(os.path.join(datadir, 'pwr@*.nc'), concat_dim='t', engine='h5netcdf')
#pwrs = xr.open_mfdataset(os.path.join(datadir, 'pwr@*.nc'), concat_dim='t')
pwrs = pwrs.assign_coords(f=lambda x: x.f/1e6)
pxx = pwrs[pwrs.autovars]
try:
    pxy = pwrs[pwrs.crossvars]
except:
    pass

vschs = ['dA', 'lA', 'dB', 'lB', 'dC', 'lC']

print('Dataset loaded in pwrs, pxx, and pxy.')
print('Spectrogram:')
print('    pft = pxx.isel(t=slice(0, 2000))')
print('    pft_db = pft.apply(lambda x: 10*np.log10(x), keep_attrs=True)')
print('    spectrograms(pft_db, ncols=3, ylabel="Frequency (MHz)", clabel="Power (dB)", vmin=0, vmax=40)')
print('Averaged spectra:')
print('    pf = pft.mean(dim="t", keep_attrs=True)')
print('    pf_db = pf.apply(lambda x: 10*np.log10(x), keep_attrs=True)')
print('    plot_dataset(ds_fftshift(pf_db), plotvars=[["dA", "dB", "dC"], ["lA", "lB", "lC"]], titles=["Dipoles", "Loops"], xlabel="Frequency (MHz)", ylabel="Power (dB)")')
print('Power timeseries, resampled:')
print('    pxxrs = pxx.resample("30s", dim="t", how="mean")')
print('    pt = pxxrs.load().median(dim="f")')
print('    pt_db = pt.apply(lambda x: 10*np.log10(x), keep_attrs=True)')
print('    plot_dataset(pt_db, plotvars=[["dA", "dB", "dC"], ["lA", "lB", "lC"]], titles=["Dipoles", "Loops"], ylabel="Power (dB)")')
