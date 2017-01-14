import os
import re
import json
import multiprocessing
import signal
import traceback
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument(
    'dirs', nargs='+',
    help='Data directory or directories.',
)
parser.add_argument(
    'savedir',
    help='Directory for saving results.',
)
parser.add_argument(
    '-0', '--t0', type=float, default=None,
    help='Start time, seconds since epoch. (default: beginning of data)',
)
parser.add_argument(
    '-1', '--t1', type=float, default=None,
    help='End time, seconds since epoch. (default: end of data)',
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
    '--nocross', action='store_true',
    help='Do not calculate cross spectra between channels. (default: %(default)s)',
)
parser.add_argument(
    '--noremovespikes', action='store_true',
    help='Do not remove spikes when processing. (default: %(default)s)',
)
parser.add_argument(
    '--clear', action='store_true',
    help='Remove previous results (*.nc files, argfile) from savedir if they exist. (default: %(default)s)',
)
parser.add_argument(
    '--overwrite', action='store_true',
    help='Overwrite previous results (*.nc files, argfile) in savedir if they conflict. (default: %(default)s)',
)

args = parser.parse_args()


import vstools

datadirs = [os.path.normpath(d) for d in args.dirs]
savedir = os.path.normpath(args.savedir)
argsfile = os.path.join(savedir, 'pwr_args')

if os.path.exists(argsfile):
    if args.clear:
        for f in os.listdir(savedir):
            if re.search('pwr@[0-9]+\.[0-9]*\.nc', f):
                os.remove(os.path.join(savedir, f))
        os.remove(argsfile)
    elif not args.overwrite:
        raise ValueError('Previous results round in save directory. Pass "--clear" to erase or "--overwrite" to add/overwrite.')
elif not os.path.exists(savedir):
    os.makedirs(savedir)

with open(argsfile, 'w') as f:
    json.dump(args.__dict__, f)

def process_data(drf_read_kwargs):
    x = vstools.drf_read_ds(**drf_read_kwargs)
    if not args.noremovespikes:
        y = x.apply(vstools.remove_spikes_da, keep_attrs=True, pfp=args.pfp)
    else:
        y = x
    ds = vstools.acs_ds(
        y, window=args.window, nperseg=args.nperseg, noverlap=args.noverlap,
        nfft=args.nfft, minints=args.minints, cross=not args.nocross,
    )
    vstools.save_ds_chunk(ds, savedir, name='pwr')

# wrapper function to print full traceback if exception is raised
def run_process_data(d):
    try:
        return process_data(d)
    except:
        raise Exception(''.join(traceback.format_exception(*sys.exc_info())))

g = vstools.drf_chunk_gen(datadirs, t0=args.t0, t1=args.t1, filesperchunk=args.filesperchunk)

print('Beginning processing. Type Ctrl + \ to send SIGQUIT and quit.')

# process first chunk serially so we get any errors
d = g.next()
process_data(d)

# now in parallel
# ignore interrupts in each worker process
def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)
pool = multiprocessing.Pool(initializer=init_worker)
resiter = pool.imap_unordered(run_process_data, g, chunksize=1)
pool.close()
# catch interrupts here and terminate
try:
    for res in resiter:
        pass
except KeyboardInterrupt:
    pool.terminate()
