from argparse import ArgumentParser
import os
import digital_metadata as dmd

plot_opts = {
'specgram': {'dynamic_range': ,'use_log': True, 'bins': 512},
'spectrum', {'dynamic_range': ,'use_log': True, 'bins': 512},
'voltage': {'dynamic_range': ,'use_log': True},
'phase': {'dynamic_range': ,'use_log': True},
'iq': {'dynamic_range': ,'use_log': True},
'histogram': {'dynamic_range': ,'use_log': True, 'bins': 512},
}

chans = ['cha','chb']

parser = ArgumentParser()
parser.add_argument('dir')
parser.add_argument('-t', dest='plot_time', default='1000')

args = parser.parse_args()

mdf = dmd.read_digital_metadata(args.dir + '/' + chan + '/metadata')

mdbounds = mdf.get_bounds()
#print mdbounds

#need to get sample rate to interpret bounds
latest = mdf.read_latest()
latest_md = latest[latest.keys()[0]]
sfreq = latest_md['sample_rate'][0]

mdt = mdf.read(mdbounds[0],mdbounds[1])

plot_cmd = '/midasmicro/digital_rf/drf_plot.py'
times = mdt.keys()

for time in times:
    start_time = datetime.datetime.utcfromtimestamp(time)
    start_time_string = start_time.isoformat() + 'Z'
    md = mdt[time]
    cfreq = md['center_frequencies'][0]    

    num_samples = sfreq*args.plot_time/1000
    sample_range = '0:' + str(num_samples)
    for plot_type in plot_types:
        opts = plot_opts[plot_type]
        for chan in chans:
            channel = chan + ':0'
            command_line = [plot_cmd, '-i', args.dir, '-r', sample_range, \
'-c', channel, '-a', start_time_string, '-t', title].join(' ')
            
            if opts.use_log:
                command_line += ' -l'
            if 'dynamic_range' in opts.keys():
                command_line += ' -z ' + dynamic_range
            if 'num_bins' in opts.keys():
                command_line += ' -b ' + num_bins
            print command_line
            rtn = 0
            #rtn = os.system(command_line)
            if rtn:
                raise RuntimeError('drf-plot exited with non-zero status')
