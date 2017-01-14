import numpy as np
import pandas
import xarray as xr
import time
import os
import glob
import re
from scipy.signal import get_window
from scipy.stats import chi2
from collections import namedtuple, OrderedDict
from itertools import combinations
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits import axes_grid1

import digital_rf_hdf5 as drf
from plotting import implot, rotate_ticklabels
from time_utils import datetime_to_float

noise_pwr_rv = chi2(2)
med_pwr_est_factor = noise_pwr_rv.mean()/noise_pwr_rv.median()

def interval_range(start, stop, step=1):
    iterable = iter(xrange(start, stop, step))
    a = b = iterable.next()
    for b in iterable:
        yield (a, b)
        a = b
    yield (b, stop)

def drf_chunk_gen(datadirs, t0=None, t1=None, filesperchunk=1):
    d = drf.read_hdf5(datadirs)

    chs = d.get_channels()
    [ss, se] = d.get_bounds(chs[0])
    fs = d.get_metadata(chs[0])['sample_rate'].value
    nperfile = int(d.get_rf_file_metadata(chs[0])['samples_per_file'][0])

    # expect user to provide time within 1e-3, so round within that below
    decimals = int(np.ceil(-np.log10((1e-3*fs)/nperfile)))

    if t0 is None:
        s0 = ss
    else:
        # floor start sample to prior file boundary, inclusive of rounded t0
        s0 = int(np.floor(np.round((t0*fs - ss)/nperfile, decimals=decimals)))*nperfile + ss
    if t1 is None:
        s1 = se
    else:
        # add 1 to floor of end sample to find subsequent file boundary,
        # inclusive of rounded t1
        s1 = int(np.floor(np.round((t1*fs - ss)/nperfile, decimals=decimals)) + 1)*nperfile + ss

    chunksize = filesperchunk*nperfile
    nchunks = int(np.ceil(float(s1 - s0)/chunksize))

    for k, (ks, ke) in enumerate(interval_range(s0, s1, chunksize)):
        print('Yielding data chunk {0}/{1}'.format(k+1, nchunks))
        nsamples = ke - ks

        yield dict(reader=d, unix_sample=ks, vector_length=nsamples)

def drf_live_chunk_gen(datadirs, filesperchunk=1):
    d = drf.read_hdf5(datadirs)

    chs = d.get_channels()
    [ss, se] = d.get_bounds(chs[0])
    fs = d.get_metadata(chs[0])['sample_rate'].value
    nperfile = int(d.get_rf_file_metadata(chs[0])['samples_per_file'][0])

    chunksize = filesperchunk*nperfile

    while True:
        d.reload()
        [ss, se] = d.get_bounds(chs[0])

        ks = max(se - chunksize, ss)
        nsamples = se - ks

        yield dict(reader=d, unix_sample=ks, vector_length=nsamples)

def drf_read_ds(reader, unix_sample, vector_length, chs=None):
    if chs is None:
        chs = reader.get_channels()
    fs = reader.get_metadata(chs[0])['sample_rate'].value
    f0 = reader.get_metadata(chs[0])['center_frequencies'][0]
    stream_args = reader.get_metadata(chs[0])['usrp_stream_args'].value
    attrs = dict(fs=fs, f0=f0, stream_args=stream_args)

    try:
        z = [reader.read_vector(unix_sample, vector_length, ch).ravel() for ch in chs]
    except IOError:
        print('Data does not exist at t={0}, reloading...'.format(unix_sample/fs))
        time.sleep(1)
        reader.reload()
        try:
            z = [reader.read_vector(unix_sample, vector_length, ch).ravel() for ch in chs]
        except IOError as e:
            print('Could not find data at t={0}.'.format(unix_sample/fs))
            raise e

    ts_ns = 1e9/fs
    freq = ts_ns*pandas.tseries.offsets.Nano()
    tidx = pandas.DatetimeIndex(
        start=unix_sample*ts_ns, periods=vector_length, freq=freq, unit='ns',
    )

    data_attrs = attrs.copy()
    data_attrs['label'] = 'Voltage'
    data_vars = OrderedDict((ch, ('t', x, data_attrs)) for ch, x in zip(chs, z))
    coords = dict(t=('t', tidx, {'label': 'Time (UTC)'}))
    data = xr.Dataset(data_vars=data_vars, coords=coords, attrs=attrs)

    return data

def _acs_shape_helper(n, nperseg=256, noverlap=None, minints=1):
    if noverlap is None:
        noverlap = nperseg//2

    # samples to step to advance to the next segment
    step = nperseg - noverlap
    # number of segments per integration, given minints
    nsegsperint = (n - noverlap)//(step*minints)
    # samples to step to advance to the next integration window
    intstep = nsegsperint*step
    # number of integrations
    nints = (n - noverlap)//intstep

    return noverlap, step, nsegsperint, intstep, nints

def acs_divisions(n, nperseg=256, noverlap=None, minints=1):
    """Number of samples used out of n when calling `acs`."""
    noverlap, step, nsegsperint, intstep, nints = _acs_shape_helper(
        n, nperseg=nperseg, noverlap=noverlap, minints=minints,
    )
    nused = nints*intstep + noverlap

    Integration = namedtuple('Integration', ['nints', 'intstep', 'noverlap'])
    Segment = namedtuple('Segment', ['nsegsperint', 'step', 'noverlap'])

    return nused, Integration(nints, intstep, noverlap), Segment(nsegsperint, step, noverlap)

def acs(x, fs=1.0, window='hanning', nperseg=256, noverlap=None, nfft=None, minints=1, cross=True):
    """Average covariance spectrum.

    """
    noverlap, step, nsegsperint, intstep, nints = _acs_shape_helper(
        x.shape[-1], nperseg=nperseg, noverlap=noverlap, minints=minints,
    )

    if nfft is None:
        nfft = nperseg

    try:
        win = get_window(window, nperseg)
    except ValueError:
        win = np.asarray(window)
        if len(win.shape) != 1:
            raise ValueError('window must be 1-D')
        if win.shape[0] != nperseg:
            raise ValueError('window must have length of nperseg')

    # normalize window so integrated power keeps same units as sample power
    win = win/np.linalg.norm(win)
    # incorporate division by nsegsperint into window so that we get mean from summation
    win = win/nsegsperint

    # shape and strides for each integration segment
    shape = x.shape[:-1] + (nsegsperint, nperseg)
    strides = x.strides[:-1] + (step*x.strides[-1], x.strides[-1])

    # number of samples needed for each integration window, including the
    # noverlap samples at the end that will also be included at the beginning
    # of the next integration window
    intlen = intstep + noverlap

    # array of end samples of integration windows
    # (+1 needed on x.shape[-1] because end values need to possibly include x.shape[-1])
    ends = np.arange(intlen, x.shape[-1]+1, intstep)

    # make time index using beginning of integration windows
    t = ((ends-intlen)*(1e9/fs)).astype('timedelta64[ns]')
    # make frequency index
    f = np.fft.fftfreq(nfft, 1/fs)

    P = np.empty(x.shape[:-2] + (x.shape[-2], x.shape[-2], nints, nfft), dtype=np.complex128)
    P.fill(np.nan)

    if cross is True:
        einsumstr = '...ilk,...jlk->...ijk'
        outview = P
    else:
        einsumstr = '...ilk,...ilk->...ik'
        outview = np.einsum('...iijk->...ijk', P)

    for kint, kend in enumerate(ends):
        intseg = x[..., (kend-intlen):kend]
        strided = np.lib.stride_tricks.as_strided(intseg, shape=shape, strides=strides)
        windowed = strided*win
        spec = np.fft.fft(windowed, nfft)
        np.einsum(einsumstr, spec, spec.conj(), out=outview[..., kint, :])

    return t, f, P

def acs_ds(ds, cross=True, **kwargs):
    x = ds.to_array(dim='ch')
    attrs = ds.attrs.copy()
    fs = attrs.pop('fs')
    t, f, P = acs(x.values, fs=fs, cross=cross, **kwargs)
    f = f + attrs['f0']

    chs = list(x.ch.values)

    t0 = x.t[0].values
    tidx = pandas.DatetimeIndex(data=t0+t)
    Pxx = np.einsum('...iijk->...ijk', P).real

    data_attrs = attrs.copy()
    data_attrs['label'] = 'Power'
    data_vars = OrderedDict((ch, (('t', 'f'), Pxx[kch], data_attrs)) for kch, ch in enumerate(chs))
    coords = dict(t=('t', tidx, {'label': 'Time (UTC)'}),
                  f=('f', f, {'label': 'Frequency (Hz)'}),
             )

    if cross is True:
        (ia, ib) = np.triu_indices(P.shape[0], 1)
        Pxy = P[..., ia, ib, :, :]
        xchs = [a + '-' + b for a, b in combinations(chs, 2)]

        data_vars.update((ch, (('t', 'f'), Pxy[kch], data_attrs)) for kch, ch in enumerate(xchs))

    ds = xr.Dataset(data_vars=data_vars, coords=coords, attrs=attrs)

    ds.attrs['autovars'] = chs
    if cross is True:
        ds.attrs['crossvars'] = xchs
    return ds

def remove_spikes(x, pfp=1e-4, mean=None):
    pwr = np.abs(x)**2
    if mean is None:
        mean = med_pwr_est_factor*np.median(pwr, axis=-1, keepdims=True)
    outlier_thresh = mean*noise_pwr_rv.ppf(1-pfp)
    x[pwr > outlier_thresh] = 0

    return x

def remove_spikes_da(x, pfp=1e-4, mean=None):
    pwr = np.abs(x.values)**2
    if mean is None:
        axis = x.dims.index('t')
        mean = med_pwr_est_factor*np.median(pwr, axis=axis, keepdims=True)
    outlier_thresh = mean*noise_pwr_rv.ppf(1-pfp)
    x.values[pwr > outlier_thresh] = 0

    return x

def save_ds_chunk(x, savedir, name='x'):
    t0_seconds = (x.t.values[0] - np.datetime64(0, 's'))/np.timedelta64(1, 's')
    fname = '{0}@{1:.3f}.nc'.format(name, t0_seconds)
    fpath = os.path.join(savedir, fname)

    if not os.path.exists(savedir):
        os.makedirs(savedir)

    x.to_netcdf(fpath, format='NETCDF4', engine='h5netcdf')

def ds_acs_to_matrix(ds):
    da = ds.to_array(dim='ch')

    d = da.sel(ch=da.autovars).values
    if hasattr(da, 'crossvars'):
        x = da.sel(ch=da.crossvars).values
        dtype = x.dtype
    else:
        dtype = d.dtype

    S = np.zeros(d.shape[:-3] + (d.shape[-3],) + d.shape[-3:], dtype=dtype)
    (ia, ib) = np.diag_indices(S.shape[-3])
    S[..., ia, ib, :, :] = d
    if hasattr(da, 'crossvars'):
        (ia, ib) = np.triu_indices(S.shape[-3], 1)
        S[..., ia, ib, :, :] = x
        S[..., ib, ia, :, :] = x.conj()

    return S

#def combine_mfdataset(paths):
    #paths = sorted(glob.glob(paths))

def ds_fftshift(ds, dim='f'):
    shift = (len(ds[dim]) + 1)//2
    return ds.roll(**{dim: shift})

def ds_ifftshift(ds, dim='f'):
    n = len(ds[dim])
    shift = n - (n + 1)//2
    return ds.roll(**{dim: shift})

def plot_dataset(ds, nrows=None, ncols=1, plotvars=None, titles=None,
                 fig_kw=dict(figsize=(14, 8)), plotfun=None, **kwargs):
    if plotvars is None:
        titles = ds.data_vars.keys()
        plotvars = [[k] for k in titles]
    if titles is None:
        titles = [None]*len(plotvars)

    if nrows is None:
        nrows = int(np.ceil(len(plotvars)/float(ncols)))

    if plotfun is None:
        def plotfun(ds, ax, xlabel=None, ylabel=None, **kwargs):
            df = ds.to_dataframe()
            df.plot(ax=ax, **kwargs)
            ax.xaxis.get_major_formatter().set_useOffset(False)
            ax.yaxis.get_major_formatter().set_useOffset(False)
            if xlabel is not None:
                ax.set_xlabel(xlabel)
            if ylabel is not None:
                ax.set_ylabel(ylabel)
            return ax

    fig = plt.figure(**fig_kw)
    gs = mpl.gridspec.GridSpec(nrows, ncols)
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

    for ax, vars_subset, title in zip(axlist, plotvars, titles):
        # ensure that vars_subset is a list, so dataset subsetting always occurs
        if not isinstance(vars_subset, list):
            vars_subset = [vars_subset]
        pds_subset = ds[vars_subset]
        plotfun(pds_subset, ax=ax, **kwargs)
        if title is not None:
            ax.set_title(title)

    return fig, axlist

def plot_spectrogram(x, t, f, ax=None, vmin=None, vmax=None, **kwargs):
    if ax is None:
        fig, ax = plt.subplots()

    cmap = plt.get_cmap()
    cmap.set_bad(color=cmap(0))

    implot_kwargs = dict(
        xistime=True, exact_ticks=False, ylabel='Frequency', clabel='Power',
        vmin=vmin, vmax=vmax,
    )
    implot_kwargs.update(kwargs)

    # extract vmin/vmax string percentage arguments into floats in dict
    # indexed by 'vmin' or 'vmax', respectively
    pct = {}
    if isinstance(vmin, str) and re.match(r'[0-9]{1,3}\.?[0-9]*\%', vmin):
        pct['vmin'] = float(vmin.rstrip('%'))
    if isinstance(vmax, str) and re.match(r'[0-9]{1,3}\.?[0-9]*\%', vmax):
        pct['vmax'] = float(vmax.rstrip('%'))
    # get desired percentiles of data and assign them to vmin/vmax kwargs
    vminmax = dict(zip(pct.keys(), np.percentile(x, pct.values())))
    implot_kwargs.update(vminmax)

    img = implot(
        np.fft.fftshift(x, axes=1), t, np.fft.fftshift(f), ax=ax,
        **implot_kwargs
    )
    ax.grid(False)
    ax.xaxis.set_tick_params(size=4.0)
    ax.xaxis.tick_bottom()
    ax.yaxis.set_tick_params(size=4.0)
    ax.yaxis.tick_left()
    plt.setp(ax.get_xticklabels(), visible=True)
    plt.setp(ax.get_yticklabels(), visible=True)
    ax.yaxis.get_major_formatter().set_useOffset(False)
    return img

def spectrograms(pds, nrows=None, ncols=1,
                 fig_kw=dict(figsize=(14, 8)), **kwargs):
    def plotfun(ds, ax, **kwargs):
        da = ds.values()[0]
        plot_spectrogram(
            da.transpose('t', 'f').values, ds.t.values, ds.f.values,
            ax=ax, **kwargs
        )
        return ax

    fig, axlist = plot_dataset(
        pds, nrows=nrows, ncols=ncols, fig_kw=fig_kw, plotfun=plotfun, **kwargs
    )

    return fig, axlist

def spectrograms_and_timeseries(pds, ax=None, fig_kw=dict(figsize=(14, 8)),
                                ptransform=lambda p: 10*np.log10(p),
                                pmin=None, pmax=None, set_vminmax=True,
                                flabel='Frequency', plabel='Power (dB)',
                                pad=0.15, cpad=0.05, cwidth=0.125, tsfrac=0.4,
                                minor=False, **spec_kwargs):
    nvars = len(pds.data_vars)

    ax_width = axes_grid1.Size.Scaled(1)
    cax_pad = axes_grid1.Size.Fixed(cpad)
    cax_width = axes_grid1.Size.Fixed(cwidth)
    h = [ax_width, cax_pad, cax_width]

    tsax_height = axes_grid1.Size.Scaled(tsfrac)
    ax_height = axes_grid1.Size.Scaled(1)
    ax_pad = axes_grid1.Size.Fixed(pad)
    v = [ax_height, ax_pad]*nvars
    v.append(tsax_height)

    if ax is None:
        fig = plt.figure(**fig_kw)
        div = axes_grid1.SubplotDivider(fig, 1, 1, 1, horizontal=h, vertical=v)
        tsloc = div.new_locator(nx=0, ny=2*nvars)
        tsax = fig.add_axes(tsloc(None, None))
    else:
        fig = ax.get_figure()
        allloc = ax.get_axes_locator()
        ss = allloc.get_subplotspec()
        div = axes_grid1.SubplotDivider(fig, ss, horizontal=h, vertical=v)

        tsloc = div.new_locator(nx=0, ny=2*nvars)
        tsax = ax
    tsax.set_axes_locator(tsloc)
    tsax.get_subplotspec = tsloc.get_subplotspec # need for tight_layout

    pt = pds.median(dim='f')
    pt_xf = pt.apply(ptransform, keep_attrs=True)
    t = pt_xf.t.values
    tepoch = t[0].astype('datetime64[D]').astype(t[0].dtype)
    t_float = datetime_to_float(t, epoch=tepoch)
    for var in pt_xf.data_vars:
        tsax.plot(t_float, pt_xf[var].values, label=var)
    tsax.set_ylim(bottom=pmin, top=pmax)
    plt.setp(tsax.get_xticklabels(), visible=False)
    tsax.locator_params(axis='y', nbins=8)
    tsax.tick_params(axis='y', labelright=True)
    tsax.get_yaxis().get_major_formatter().set_useOffset(False)
    if minor:
        tsax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    tsax.grid(True, which='major', linewidth=1.5)
    tsax.grid(True, which='minor', linewidth=0.75)
    tsax.set_ylabel(plabel)
    tsax.legend(loc='best', ncol=(nvars//3 + 1))

    s_kwargs = dict()
    if set_vminmax:
        vmin, vmax = tsax.get_ylim()
        spec_kwargs['vmin'] = vmin
        spec_kwargs['vmax'] = vmax
    s_kwargs.update(spec_kwargs)
    axlist = [tsax]
    for k, var in enumerate(pds.data_vars):
        loc = div.new_locator(nx=0, ny=2*(nvars-(k+1)))
        if len(axlist) == 1:
            ax = fig.add_axes(loc(None, None), sharex=tsax)
        else:
            ax = fig.add_axes(loc(None, None), sharex=tsax, sharey=axlist[1])
        ax.set_axes_locator(loc)
        ax.get_subplotspec = loc.get_subplotspec # need for tight_layout
        axlist.append(ax)

        cbloc = div.new_locator(nx=2, ny=2*(nvars-(k+1)))
        cax = fig.add_axes(cbloc(None, None))
        cax.set_axes_locator(cbloc)
        cax.get_subplotspec = cbloc.get_subplotspec # need for tight_layout

        da = pds[var]

        img = plot_spectrogram(
            ptransform(da.transpose('t', 'f').values),
            da.t.values, da.f.values,
            ax=ax, cbar=False,
            ylabel=flabel, **s_kwargs
        )
        cb = fig.colorbar(img, cax=cax, ax=ax, orientation='vertical')
        cb.set_label(var + ' ' + plabel)
        img.colorbar = cb
        ax.get_yaxis().get_major_formatter().set_useOffset(False)

        if k < nvars-1:
            plt.setp(ax.get_xticklabels(), visible=False)
            ax.set_xlabel('')
            ax.xaxis.set_tick_params(size=0)
    rotate_ticklabels(ax.xaxis, rotation=0)

    return axlist
