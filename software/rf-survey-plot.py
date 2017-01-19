from argparse import ArgumentParser
import os
import digital_metadata as dmd
import datetime

plot_opts = {
'specgram': {'use_log': True, 'bins': '1024'},
'spectrum': {'use_log': True, 'bins': '1024'},
'voltage': {},
'phase': {'bins': '128', 'use_log': True},
'iq': {},
'histogram': {'use_log': True, 'bins': '128'}
}

chans = ['cha','chb']

parser = ArgumentParser()
parser.add_argument('dir')
parser.add_argument('-t', dest='plot_time', default='1000')
parser.add_argument('-s', dest='save_dir', default='')
args = parser.parse_args()

#plot_cmd = '~/midasmicro/digital_rf/drf_plot.py'
plot_cmd = 'python ~/rapid_deployment/software/drf_plot.py'
i = 0
for chan in chans:
    chan_string = chan + ':0'
    mdf = dmd.read_digital_metadata(args.dir + chan + '/metadata')

    mdbounds = mdf.get_bounds()
    print mdbounds

    #need to get sample rate to interpret bounds
    latest = mdf.read_latest()
    latest_md = latest[latest.keys()[0]]
    sfreq = latest_md['sample_rate']

    mdt = mdf.read(mdbounds[0],mdbounds[1])
    times = mdt.keys()

    for time in times:
        print "time: " + str(time)
        start_time = datetime.datetime.utcfromtimestamp(time/sfreq)
        start_time_string = start_time.isoformat()
        md = mdt[time]
        cfreq = int(md['center_frequencies'].tolist()[0][0])
        print "cfreq: " + str(cfreq)
        print "sfreq: " + str(sfreq)
        num_samples = int(sfreq*int(args.plot_time)/1000)
        sample_range = '0:{:d}'.format(num_samples)
        for plot_type in plot_opts.keys():
            opts = plot_opts[plot_type]
            command_line = ' '.join([plot_cmd, '-i', args.dir, '-p', plot_type, '-r', \
sample_range, '-c', chan_string, '-a', start_time_string])
            
            if 'use_log' in opts.keys() and opts['use_log']:
                command_line += ' -l'
            if 'dynamic_range' in opts.keys():
                command_line += ' -z ' + opts['dynamic_range']
            if 'num_bins' in opts.keys():
                command_line += ' -b ' + opts['num_bins']
            if args.save_dir != "":
                fname = '_'.join([plot_type,chan,str(cfreq)])
                command_line += ' -s ' + args.save_dir + fname + '.png'
            print command_line
            #rtn = 0
            #i += 1
            #if i > 5:
            #    quit()
            rtn = os.system(command_line)
            if rtn:
                raise RuntimeError('drf-plot exited with non-zero status')
