"""
 drf_plot.py

 $Id$

 Simple program to load 16 bit IQ data and make some basic plots. Command
 line options are supported and data frames may be filtered from the output. The
 program can offset into a data file to limit the memory usage when plotting
 a subset of a data file.

"""

import sys
import os
import string
import struct
import getopt
import traceback
import glob
import time
import calendar

import numpy
import matplotlib
matplotlib.use('GTKAgg')
import matplotlib.pylab

import digital_rf_hdf5
import digital_metadata as dmd

import pandas as pd
import xarray as xr

def voltage_process(data, sfreq, toffset, modulus, integration, log_scale, title):
    """ Break voltages by modulus and display each block. Integration here acts
    as a pure average on the voltage level data.
    """
    if modulus:
        block = 0
        block_size = integration*modulus
        block_toffset = toffset
        while block < len(data) / (block_size):
            dblock = data[block*block_size:block*block_size+modulus]
            # complete integration
            for idx in range(1,integration):
                dblock += data[block*block_size+idx*modulus:block*block_size+idx*modulus+modulus]

            dblock /= integration

            voltage_plot(dblock, sfreq, block_toffset, log_scale, title)

            block += 1
            block_toffset += block_size / sfreq

    else:
        voltage_plot(data, sfreq, toffset, log_scale, title)

"""

    Plot the real and imaginary voltage from IQ data.

"""
def voltage_plot(data, sfreq, toffset, log_scale, title):

    print("voltage")

    t_axis = numpy.arange(0,len(data))/sfreq + toffset

    matplotlib.pylab.subplot(211)
    matplotlib.pylab.plot(t_axis,data.real)
    matplotlib.pylab.grid(1)
    maxr = numpy.max(data.real)
    minr = numpy.min(data.real)

    if minr == 0.0 and maxr == 0.0:
        minr = -1.0
        maxr = 1.0

    matplotlib.pylab.axis([t_axis[0],t_axis[len(t_axis)-1],minr, maxr]);
    matplotlib.pylab.ylabel('I sample value (A/D units)')

    matplotlib.pylab.subplot(212)
    matplotlib.pylab.plot(t_axis,data.imag)
    matplotlib.pylab.grid(1)
    maxi = numpy.max(data.imag)
    mini = numpy.min(data.imag)

    if mini == 0.0 and maxi == 0.0:
        mini = -1.0
        maxi = 1.0

    matplotlib.pylab.axis([t_axis[0],t_axis[len(t_axis)-1],mini,maxi]);

    matplotlib.pylab.xlabel('time (seconds)')

    matplotlib.pylab.ylabel('Q sample value (A/D units)')
    matplotlib.pylab.title(title)

    #matplotlib.pylab.show()



def power_process(data, sfreq, toffset, modulus, integration, log_scale, zscale, title):
    """ Break power by modulus and display each block. Integration here acts
    as a pure average on the power level data.
    """
    if modulus:
        block = 0
        block_size = integration*modulus
        block_toffset = toffset
        while block < len(data) / block_size:

            vblock = data[block*block_size:block*block_size+modulus]
            pblock = vblock * numpy.conjugate(vblock)

            # complete integration
            for idx in range(1,integration):

                vblock = data[block*block_size+idx*modulus:block*block_size+idx*modulus+modulus]
                pblock += vblock * numpy.conjugate(vblock)

            pblock /= integration

            power_plot(pblock, sfreq, block_toffset, log_scale, zscale, title)

            block += 1
            block_toffset += block_size / sfreq

    else:
        pdata = data * numpy.conjugate(data)
        power_plot(pdata, sfreq, toffset, log_scale, zscale, title)



"""

    Plot the computed power of the iq data.

"""
def power_plot(data, sfreq, toffset, log_scale, zscale, title):

    print("power")

    t_axis = numpy.arange(0,len(data))/sfreq + toffset

    if log_scale:
        lrxpwr = 10*numpy.log10(data + 1E-12)
    else:
        lrxpwr = data

    zscale_low, zscale_high = zscale

    if zscale_low == 0 and zscale_high == 0:
        if log_scale:
            zscale_low = numpy.min(lrxpwr[numpy.where(lrxpwr.real != -numpy.Inf)])
            zscale_high = numpy.max(lrxpwr) + 3.0
        else:
            zscale_low = numpy.min(lrxpwr)
            zscale_high = numpy.max(lrxpwr)

    matplotlib.pylab.plot(t_axis,lrxpwr.real)
    matplotlib.pylab.grid(1)
    matplotlib.pylab.axis([toffset,t_axis[len(t_axis)-1],zscale_low,zscale_high]);

    matplotlib.pylab.xlabel('time (seconds)')
    matplotlib.pylab.ylabel('power (dB)')
    matplotlib.pylab.title(title)

    #matplotlib.pylab.show()


def iq_process(data, sfreq,  toffset, modulus, integration, log_scale, title):
    """ Break voltages by modulus and display each block. Integration here acts
    as a pure average on the voltage level data prior to iq plotting.
    """
    if modulus:
        block = 0
        block_size = integration*modulus
        block_toffset = toffset
        while block < len(data) / block_size:
            dblock = data[block*block_size:block*block_size+modulus]
            # complete integration
            for idx in range(1,integration):
                dblock += data[block*block_size+idx*modulus:block*block_size+idx*modulus+modulus]

            dblock /= integration

            iq_plot(dblock, block_toffset, log_scale, title)

            block += 1
            block_toffset += block_size / sfreq

    else:
        iq_plot(data, toffset, log_scale, title)


"""

    Plot an IQ circle from the data in linear or log scale.

"""
def iq_plot(data, toffset, log_scale, title):

    print("iq")

    if log_scale:
        rx_raster_r = numpy.sign(data.real) * numpy.log10(numpy.abs(data.real)+1E-30) / numpy.log10(2.)
        rx_raster_i = numpy.sign(data.imag) * numpy.log10(numpy.abs(data.imag)+1E-30) / numpy.log10(2.)
    else:
        data *= 1.0 / 32768.0
        rx_raster_r = data.real
        rx_raster_i = data.imag

    matplotlib.pylab.plot(rx_raster_r, rx_raster_i,'.')

    axmx = numpy.max([numpy.max(rx_raster_r), numpy.max(rx_raster_i)])

    matplotlib.pylab.axis([-axmx, axmx, -axmx, axmx])
    matplotlib.pylab.grid(1)
    matplotlib.pylab.xlabel('I')
    matplotlib.pylab.ylabel('Q')
    matplotlib.pylab.title(title)

    #matplotlib.pylab.show()

def phase_process(data, sfreq,  toffset, modulus, integration, log_scale, title):
    """ Break voltages by modulus and display the phase of each block. Integration here acts
    as a pure average on the voltage level data prior to iq plotting.
    """
    if modulus:
        block = 0
        block_size = integration*modulus
        block_toffset = toffset
        while block < len(data) / block_size:
            dblock = data[block*block_size:block*block_size+modulus]
            # complete integration
            for idx in range(1,integration):
                dblock += data[block*block_size+idx*modulus:block*block_size+idx*modulus+modulus]

            dblock /= integration

            phase_plot(dblock, block_toffset, log_scale, title)

            block += 1
            block_toffset += block_size / sfreq

    else:
        phase_plot(data, toffset, log_scale, title)

"""

    Plot the phase of the data in linear or log scale.

"""
def phase_plot(data, toffset, log_scale, title):

    print("phase")

    phase = numpy.angle(data) / numpy.pi

    matplotlib.pylab.plot(phase)

    axmx = numpy.max(phase)

    #matplotlib.pylab.axis([-axmx, axmx, -axmx, axmx])
    matplotlib.pylab.grid(1)
    matplotlib.pylab.xlabel('time')
    matplotlib.pylab.ylabel('phase')
    matplotlib.pylab.title(title)

    #matplotlib.pylab.show()



def spectrum_process(data, sfreq, cfreq, toffset, modulus, integration, bins, log_scale, zscale, detrend, title, clr):
    """ Break spectrum by modulus and display each block. Integration here acts
    as a pure average on the spectral data.
    """

    if detrend:
        dfn = matplotlib.mlab.detrend_mean
    else:
        dfn = matplotlib.mlab.detrend_none

    win = numpy.blackman(bins)

    if modulus:
        block = 0
        block_size = integration*modulus
        block_toffset = toffset
        while block < len(data) / block_size:

            vblock = data[block*block_size:block*block_size+modulus]
            pblock,freq = matplotlib.mlab.psd(vblock,NFFT=bins,Fs=sfreq,detrend=dfn,window=win,scale_by_freq=False)

            # complete integration
            for idx in range(1,integration):

                vblock = data[block*block_size+idx*modulus:block*block_size+idx*modulus+modulus]
                pblock_n,freq = matplotlib.mlab.psd(vblock,NFFT=bins,Fs=sfreq,detrend=dfn,window=matplotlib.mlab.window_hanning,scale_by_freq=False)
                pblock += pblock_n

            pblock /= integration

            spectrum_plot(pblock, freq, block_toffset, log_scale, zscale, title, clr)

            block += 1
            block_toffset += block_size / sfreq

    else:
        pdata,freq = matplotlib.mlab.psd(data,NFFT=bins,Fs=sfreq,detrend=dfn,window=win,scale_by_freq=False)
        spectrum_plot(pdata, freq, cfreq, toffset, log_scale, zscale, title, clr)


"""
    Plot a spectrum from the data for a given fft bin size.

"""
def spectrum_plot(data, freq, cfreq, toffset, log_scale, zscale, title,clr):

    print("spectrum")
    tail_str = ''
    if log_scale:
#        pss = 10.0*numpy.log10(data / numpy.max(data))
	pss = 10.0*numpy.log10(data + 1E-12)
	tail_str = ' (dB)'
    else:
        pss = data

    print freq
    freq_s = freq/1.0E6 + cfreq /1.0E6
    print freq_s
    zscale_low, zscale_high = zscale

    if zscale_low == 0 and zscale_high == 0:
        if log_scale:
            zscale_low = numpy.median(numpy.min(pss[numpy.where(pss.real != -numpy.Inf)])) - 3.0
            zscale_high = numpy.median(numpy.max(pss)) + 3.0
        else:
            zscale_low = numpy.median(numpy.min(pss))
            zscale_high = numpy.median(numpy.max(pss))

    matplotlib.pylab.plot(freq_s,pss,clr)
    print freq_s[0], freq_s[-1], zscale_low, zscale_high
    matplotlib.pylab.axis([freq_s[0], freq_s[-1], zscale_low, zscale_high])
    matplotlib.pylab.grid(1)
    matplotlib.pylab.xlabel('frequency (MHz)')
    matplotlib.pylab.ylabel('power spectral density' + tail_str, fontsize=12)
    matplotlib.pylab.title(title)

    #matplotlib.pylab.show()

def histogram_process(data, sfreq, toffset, modulus, integration, bins, log_scale, title):
    """ Break voltages by modulus and display each block. Integration here acts
    as a pure average on the voltage level data prior to the histogram.
    """
    if modulus:
        block = 0
        block_size = integration*modulus
        block_toffset = toffset
        while block < len(data) / block_size:
            dblock = data[block*block_size:block*block_size+modulus]
            # complete integration
            for idx in range(1,integration):
                dblock += data[block*block_size+idx*modulus:block*block_size+idx*modulus+modulus]

            dblock /= integration

            histogram_plot(dblock, sfreq, block_toffset, bins, log_scale, title)

            block += 1
            block_toffset += block_size / sfreq

    else:
        histogram_plot(data, sfreq, toffset, bins, log_scale, title)


"""
    Plot a histogram of the data for a given bin size.

"""
def histogram_plot(data, sfreq, toffset, bins, log_scale, title):

    print("histogram")

    matplotlib.pylab.hist(numpy.real(data),bins,log=log_scale,histtype='bar',color=['green'])
    matplotlib.pylab.hold(1)
    matplotlib.pylab.hist(numpy.imag(data),bins,log=log_scale,histtype='bar',color=['blue'])
    matplotlib.pylab.grid(1)
    matplotlib.pylab.xlabel('adc value')
    matplotlib.pylab.ylabel('frequency')
    matplotlib.pylab.title(title)
    matplotlib.pylab.hold(0)

    #matplotlib.pylab.show()

def specgram_process(data, sfreq, cfreq, toffset, modulus, integration, bins, detrend, log_scale, zscale, title):
    """ Break spectrum by modulus and display each block. Integration here acts
    as a pure average on the spectral data.
    """

    if detrend:
        dfn = matplotlib.mlab.detrend_mean
    else:
        dfn = matplotlib.mlab.detrend_none

    noverlap = int(bins * 0.9)

    if modulus:
        block = 0
        block_size = integration*modulus
        block_toffset = toffset
        while block < len(data) / block_size:

            vblock = data[block*block_size:block*block_size+modulus]
            pblock,freq, tm = matplotlib.mlab.specgram(vblock, NFFT=bins, Fs=sfreq, noverlap=noverlap, detrend=dfn,window=matplotlib.mlab.window_hanning, scale_by_freq=False)

            # complete integration
            for idx in range(1,integration):

                vblock = data[block*block_size+idx*modulus:block*block_size+idx*modulus+modulus]
                pblock_n,freq, tm = matplotlib.mlab.specgram(vblock, NFFT=bins, Fs=sfreq, noverlap=noverlap, detrend=dfn,window=matplotlib.mlab.window_hanning, scale_by_freq=False)
                #print idx, numpy.shape(pblock), numpy.shape(pblock_n)
                pblock += pblock_n

            pblock /= integration

            extent = [block_toffset, (block_size/sfreq + block_toffset), -sfreq/2.0E6+cfreq/1.0E6, sfreq/2.0E6 + cfreq/1.0E6]

            specgram_plot(pblock, extent, log_scale, zscale, title)

            block += 1
            block_toffset += block_size / sfreq

    else:
        pdata,freq,tm =  matplotlib.mlab.specgram(data, NFFT=bins, Fs=sfreq, noverlap=noverlap, detrend=dfn,window=matplotlib.mlab.window_hanning, scale_by_freq=False)

        extent = [toffset, (len(data)/sfreq + toffset), -sfreq/2.0E6+cfreq/1.0E6, sfreq/2.0E6+cfreq/1.0E6]

        specgram_plot(pdata, extent, log_scale, zscale, title)


"""

    Plot a specgram from the data for a given fft size.

"""
def specgram_plot(data, extent, log_scale, zscale, title):

    print("specgram")

    # set to log scaling
    if log_scale:
        Pss = 10.0*numpy.log10(data + 1E-12)
    else:
        Pss = data

    # scale for zero centered kilohertz
    # determine image x-y extent

    # determine image color extent in log scale units
    zscale_low, zscale_high = zscale

    if zscale_low == 0 and zscale_high == 0:
        if log_scale:
            zscale_low = numpy.median(numpy.min(Pss)) - 3.0
            zscale_high = numpy.median(numpy.max(Pss)) + 10.0
        else:
            zscale_low = numpy.median(numpy.min(Pss))
            zscale_high = numpy.median(numpy.max(Pss))

    vmin = zscale_low
    vmax = zscale_high

    matplotlib.pylab.imshow(Pss, cmap=matplotlib.cm.jet, origin='lower', extent=extent, interpolation='nearest', vmin=vmin, vmax=vmax, aspect='auto')
    matplotlib.pylab.colorbar()
    matplotlib.pylab.xlabel('time (seconds)')
    matplotlib.pylab.ylabel('frequency (MHz)', fontsize=12)
    # matplotlib.pylab.title(title)
    #matplotlib.pylab.show()

def rti_process(data, sfreq, toffset, modulus, integration, detrend, log_scale, zscale, title):
    """ Break power by modulus and make an RTI stripe for each block. Integration here acts
    as a pure average on the power level data.
    """
    if modulus:
        block = 0
        block_size = integration*modulus
        block_toffset = toffset

        rti_bins = len(data) / block_size

        RTIdata = numpy.zeros([modulus,rti_bins],numpy.complex64)
        RTItimes = numpy.zeros([rti_bins,])

        while block < len(data) / block_size:

            vblock = data[block*block_size:block*block_size+modulus]
            pblock = vblock * numpy.conjugate(vblock)

            # complete integration
            for idx in range(1,integration):

                vblock = data[block*block_size+idx*modulus:block*block_size+idx*modulus+modulus]
                pblock += vblock * numpy.conjugate(vblock)

            pblock /= integration
            if detrend:
                pblock -= numpy.mean(pblock)

            # load RTI stripe
            # preclear the row to be written to eliminate old values from longer rasters
            RTIdata[:,block] = 0.0
            # write the latest row of data
            RTIdata[0:len(pblock),block] = pblock
            RTItimes[block] = float(block_toffset)

            block += 1
            block_toffset += block_size / sfreq
    else:
        print "Must have a modulus for an RTI!"
        return

    # create time axis
    tick_spacing = numpy.arange(0, rti_bins, rti_bins/len(RTItimes))
    tick_labels = []

    for s in tick_spacing:
        tick_time = RTItimes[s]

        if tick_time == 0:
            tick_string = ''
        else:
            tick_string = '%04.3f' % (tick_time)

        tick_labels.append(tick_string)


    # create a range axis
    rx_axis = numpy.arange(0,modulus)*0.15 # km per microsecond

    range_scale = 1.0E6/sfreq # sampling period in microseconds

    rx_axis *= range_scale # km range scale

    # determine image x-y extent
    extent = 0, rti_bins, 0, numpy.max(rx_axis)

    rti_plot(RTIdata.real,extent, tick_spacing, tick_labels, log_scale, zscale, title)


def rti_plot(data, extent, tick_spacing, tick_labels, log_scale, zscale, title):


    # set to log scaling
    if log_scale:
        RTId = 10.0*numpy.log10(data)
    else:
        RTId = data

    zscale_low, zscale_high = zscale
    if zscale_low == 0 and zscale_high == 0:
        if log_scale:
            zscale_low = numpy.median(numpy.min(RTId[numpy.where(RTId.real != -numpy.Inf)])) - 3.0
            zscale_high = numpy.median(numpy.max(RTId)) + 10.0
        else:
            zscale_low = numpy.median(numpy.min(RTId))
            zscale_high = numpy.median(numpy.max(RTId))

    vmin = zscale_low
    vmax = zscale_high


    matplotlib.pylab.imshow(RTId, cmap=matplotlib.cm.jet, origin='lower', extent=extent, interpolation='nearest', vmin=vmin, vmax=vmax,aspect='auto')

    # plot dates

    matplotlib.pylab.xticks(tick_spacing, tick_labels, rotation=-45, fontsize=10)
    matplotlib.pylab.colorbar()
    matplotlib.pylab.xlabel('time (seconds)', fontsize=12)
    matplotlib.pylab.ylabel('range (km)', fontsize=12)
    matplotlib.pylab.title(title)
    #matplotlib.pylab.show()

def sti_process(data, sfreq, toffset, modulus, integration, bins, detrend, log_scale, zscale, title):
    """ Break data by modulus and make an STI stripe for each block. Integration here acts
    as a pure average on the spectrum level data.
    """

    if detrend:
        dfn = matplotlib.mlab.detrend_mean
    else:
        dfn = matplotlib.mlab.detrend_none

    if modulus:
        block = 0
        block_size = integration*modulus
        block_toffset = toffset

        sti_bins = len(data) / block_size

        STIdata = numpy.zeros([bins,sti_bins],numpy.complex64)
        STItimes = numpy.zeros([sti_bins,])

        while block < len(data) / block_size:

            vblock = data[block*block_size:block*block_size+modulus]
            pblock,freq = matplotlib.mlab.psd(vblock,NFFT=bins,Fs=sfreq,detrend=dfn,window=matplotlib.mlab.window_hanning,scale_by_freq=False)

            # complete integration
            for idx in range(1,integration):

                vblock = data[block*block_size+idx*modulus:block*block_size+idx*modulus+modulus]
                pblock_n,freq = matplotlib.mlab.psd(vblock,NFFT=bins,Fs=sfreq,detrend=dfn,window=matplotlib.mlab.window_hanning,scale_by_freq=False)
                pblock += pblock_n

            pblock /= integration

            # load RTI stripe
            # preclear the row to be written to eliminate old values from longer rasters
            STIdata[:,block] = 0.0
            # write the latest row of data
            STIdata[0:len(pblock),block] = pblock
            STItimes[block] = float(block_toffset)

            block += 1
            block_toffset += block_size / sfreq
    else:
        print "Must have a modulus for an STI!"
        return

    # create time axis
    tick_spacing = numpy.arange(0, sti_bins, sti_bins/len(STItimes))
    tick_labels = []

    for s in tick_spacing:
        tick_time = STItimes[s]

        if tick_time == 0:
            tick_string = ''
        else:
            tick_string = '%04.3f' % (tick_time)

        tick_labels.append(tick_string)

    axf = -sfreq/2000.0
    extent = [0, sti_bins, -axf, axf]

    sti_plot(STIdata.real, freq, extent, tick_spacing, tick_labels, log_scale, zscale, title)


def sti_plot(data, freq, extent, tick_spacing, tick_labels, log_scale, zscale, title):

    pss = data
    freq_s = freq/1000.0

    # set to log scaling
    if log_scale:
        STId = 10.0*numpy.log10(pss)
    else:
        STId = pss

    zscale_low, zscale_high = zscale
    if zscale_low == 0 and zscale_high == 0:
        if log_scale:
            zscale_low = numpy.median(numpy.min(STId)) - 3.0
            zscale_high = numpy.median(numpy.max(STId)) + 10.0
        else:
            zscale_low = numpy.median(numpy.min(STId))
            zscale_high = numpy.median(numpy.max(STId))

    axf = numpy.max(numpy.abs(freq_s))


    vmin = zscale_low
    vmax = zscale_high


    matplotlib.pylab.imshow(STId, cmap=matplotlib.cm.jet, origin='lower', extent=extent, interpolation='nearest', vmin=vmin, vmax=vmax,aspect='auto')

    # plot dates

    matplotlib.pylab.colorbar()
    matplotlib.pylab.xticks(tick_spacing, tick_labels, rotation=-45, fontsize=10)
    matplotlib.pylab.xlabel('time (seconds)', fontsize=12)
    matplotlib.pylab.ylabel('frequency (kHz)', fontsize=12)
    matplotlib.pylab.title(title)
    #matplotlib.pylab.show()


def hex2vec(h, ell):
    """hex2vec(h, ell) generates sign vector of length ell from the hex string h.
    ell must be <= 4*len(h) (excluding the optional leading "0x")
    """

    if h[0:2] in [ '0x', '0X' ]:
        h = h[2:]

    nybble = numpy.array([
                          [0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 0, 1, 1],
                          [0, 1, 0, 0], [0, 1, 0, 1], [0, 1, 1, 0], [0, 1, 1, 1],
                          [1, 0, 0, 0], [1, 0, 0, 1], [1, 0, 1, 0], [1, 0, 1, 1],
                          [1, 1, 0, 0], [1, 1, 0, 1], [1, 1, 1, 0], [1, 1, 1, 1]])

    vec = numpy.ravel(numpy.array([ nybble[int(x, 16)] for x in h ]))

    if len(vec) < ell:
        raise ValueError('hex string too short')
    return vec[ len(vec)-ell : ]

def apply_msl_filter(data, msl_code_length, msl_baud_length):

    code_table = {2:[0,1],  # barker codes
                  3:[0,0,1],
                  4:[0,1,0,0],
                  5:[0,0,0,1,0],
                  7:[0,0,0,1,1,0,1],
                  11:[0,0,0,1,1,1,0,1,1,0,1],
                  13:[0,0,0,0,0,1,1,0,0,1,0,1,0],
                  14:[0,1,0,1,0,0,1,0,0,0,0,0,1,1], # MSL codes from here down
                  15:[0,0,1,1,0,0,0,0,0,1,0,1,0,1,1],
                  16:[0,1,1,0,1,0,0,0,0,1,1,1,0,1,1,1],
                  17:[0,0,1,1,1,0,1,1,1,0,1,0,0,1,0,1,1],
                  18:[0,1,1,0,0,1,0,0,0,0,1,1,1,1,0,1,0,1],
                  19:[1,0,1,1,0,1,1,1,0,1,1,1,0,0,0,1,1,1,1],
                  20:[0,1,0,1,0,0,0,1,1,0,0,0,0,0,0,1,1,0,1,1],
                  21:[1,0,1,1,0,1,0,1,1,1,0,1,1,1,0,0,0,0,0,1,1],
                  22:[0,0,1,1,1,0,0,1,1,0,1,1,0,1,0,1,0,1,1,1,1,1],
                  23:[0,1,1,1,0,0,0,1,1,1,1,1,1,0,1,0,1,0,0,1,0,0,1],
                  24:[0,1,1,0,0,1,0,0,1,0,1,0,1,1,1,1,1,1,1,0,0,0,1,1],
                  25:[1,0,0,1,0,0,1,0,1,0,1,0,0,0,0,0,0,1,1,0,1,1,0,0,1],
                  26:[1,0,0,0,1,1,1,0,0,0,0,0,0,0,1,0,1,0,1,1,0,1,1,0,0,1],
                  27:[0,1,0,0,1,0,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,0,0,0,1,1,1],
                  28:hex2vec('0xDA44478',28),
                  29:[1,0,1,1,0,0,1,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,1,1,0,0,1,1,1],
                  30:[1,0,0,0,1,1,0,0,0,1,0,1,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,1,1,1],
                  31:[0,1,0,1,0,1,0,0,1,0,0,1,0,0,1,1,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1],
                  32:[0,0,0,1,1,1,0,1,1,0,0,0,1,1,0,1,0,0,1,1,1,0,1,1,1,1,1,1,0,0,1,1],
                  33:[0,1,1,0,0,1,1,0,0,1,0,1,0,1,0,1,0,0,1,0,1,1,0,0,0,0,1,1,1,1,1,1,1],
                  34:[1,1,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,0,1,0,0,1,0,1,0,1,0,1],
                  49:hex2vec('0x012ABEC79E46F',49),
                  50:hex2vec('0x025863ABC266F',50),
                  51:hex2vec('0x71C077376ADB4',51),
                  52:hex2vec('0x0945AE0F3246F',52),
                  53:hex2vec('0x0132AA7F8D2C6F',53),
                  54:hex2vec('0x0266A2814B3C6F',54),
                  55:hex2vec('0x04C26AA1E3246F',55),
                  56:hex2vec('0x099BAACB47BC6F',56),
                  57:hex2vec('0x01268A8ED623C6F',57),
                  58:hex2vec('0x023CE545C9ED66F',58),
                  59:hex2vec('0x049D38128A1DC6F',59),
                  60:hex2vec('0x0AB8DF0C973252F',60),
                  61:hex2vec('0x00459CC5F4694BAF',61),
                  62:hex2vec('0x008B1B5318BE4BAF',62),
                  63:hex2vec('0x04CF5A2471657C6F',63),
                  64:hex2vec('0x122B21E9978C2BBF',64),
                  65:hex2vec('0x1045A6A6270AC4BBF',65),
                  66:hex2vec('0x088B8CF1325A50BBF',66),
                  67:hex2vec('0x01169F29AC67C4B9F',67),
                  68:hex2vec('0x122B43963E8662BBF',68),
                  69:hex2vec('0x1D9024F657C5EE71EA',69),
                  70:hex2vec('0x1A133B4E3093EDD57E',70),
                  71:hex2vec('0x63383AB6B452ED93FE',71),
                  72:hex2vec('0xE4CD5AF0D054433D82',72),
                  73:hex2vec('0x1B66B26359C3E2BC00A',73),
                  74:hex2vec('0x36DDBED681F98C70EAE',74),
                  75:hex2vec('0x6399C983D03EFDB556D',75),
                  76:hex2vec('0xDB69891118E2C2A1FA0',76),
                  77:hex2vec('0x1961AE251DC950FDDBF4',77),
                  78:hex2vec('0x328B457F0461E4ED7B73',78),
                  79:hex2vec('0x76CF68F327438AC6FA80',79),
                  80:hex2vec('0xCE43C8D986ED429F7D75',80),
                  81:hex2vec('0x0E3C32FA1FEFD2519AB32',81),
                  82:hex2vec('0x3CB25D380CE3B7765695F',82),
                  83:hex2vec('0x711763AE7DBB8482D3A5A',83),
                  84:hex2vec('0xCE79CCCDB6003C1E95AAA',84),
                  85:hex2vec('0x19900199463E51E8B4B574',85),
                  86:hex2vec('0x3603FB659181A2A52A38C7',86),
                  87:hex2vec('0x7F7184F04F4E5E4D9B56AA',87),
                  88:hex2vec('0x9076589AF5702502CE2CE2',88),
                  89:hex2vec('0x180E09434E1BBC44ACDAC8A',89),
                  90:hex2vec('0x3326D87C3A91DA8AFA84211',90),
                  91:hex2vec('0x77F80E632661C3459492A55',91),
                  92:hex2vec('0xCC6181859D9244A5EAA87F0',92),
                  93:hex2vec('0x187B2ECB802FB4F56BCCECE5',93),
                  94:hex2vec('0x319D9676CAFEADD68825F878',94),
                  95:hex2vec('0x69566B2ACCC8BC3CE0DE0005',95),
                  96:hex2vec('0xCF963FD09B1381657A8A098E',96),
                  97:hex2vec('0x1A843DC410898B2D3AE8FC362',97),
                  98:hex2vec('0x30E05C18A1525596DCCE600DF',98),
                  99:hex2vec('0x72E6DB6A75E6A9E81F0846777',99),
                  100:hex2vec('0xDF490FFB1F8390A54E3CD9AAE',100),
                  101:hex2vec('0x1A5048216CCF18F83E910DD4C5',101),
                  102:hex2vec('0x2945A4F11CE44FF664850D182A',102),
                  103:hex2vec('0x77FAAB2C6E065AC4BE18F274CB',103),
                  104:hex2vec('0xE568ED4982F9660EBA2F611184',104),
                  105:hex2vec('0x1C6387FF5DA4FA325C895958DC5',105),
                  }

    print("msl filter data with code ", msl_code_length)
    code = code_table[msl_code_length]
    # note that the codes are time reversed as defined compared to what we send
    x_msl = numpy.zeros(msl_baud_length*len(code),dtype=numpy.complex64)
    idx = 0
    for c in code:
        block = len(x_msl)/len(code)
        if c == 1:
            x_msl[idx * block:idx*block+block] = numpy.ones(block,dtype=numpy.float64) + 0j
        else:
            x_msl[idx * block:idx*block+block] = -1.0*numpy.ones(block,dtype=numpy.float64) + 0j
        idx += 1

    #print("convolving filter ", x_msl)
    #print(len(x_msl))
    dc = numpy.convolve(d,x_msl,'same')

    return dc


# load data and plot it!

def usage():
    print("usage : %s" % sys.argv[0])
    print("        -i <input file> -p <type> [-c <chan>] [-f <ch>:<size>] [-r <range>] [-b <bins>] [-t <title>] [-l] [-d] [-m <code info>]")
    print(" -i <input file>        The name of the file to load and display, may be a wildcard for multiple files.")
    print(" -o <offset freq>       Center frequency offset.")
    print(" -p <type>              The type of plot to make : power, iq, voltage, histogram, spectrum, specgram, rti.")
    print(" -c <chan>:<subchan>    Channel and subchannel to plot for mutliple channel data, default is first channel. (e.g. test:0 or test:1, string:integer)")
    print(" -a <time>              Absolute start time for 'zero' sample (seconds). <year>-<month>-<day>T<hour>:<minute>:<second> in UT")
    print(" -r <start>:<stop>:[<skip>]:[<modulus>]:[<integration>]  Data range to display in samples. start:stop:modulus:integration, full no modulus by default.")
    print(" -z <low>:<high>        Dynamic range setting for log plots e.g. -50:0 in dB.")
    print(" -b <bins>              Number of bins for histogram, spectral, specgram, and rti modes.")
    print(" -t <title>             A string \'Data from 2009-01-21 on 100.0 MHz FM\' for a title. The single quotes are important.")
    print(" -l                     Enable log scale for plots.")
    print(" -d                     Enable detrending for spectral estimates.")
    print(" -m <code_length>:<baud_length>  MSL code length (bauds) and samples per baud for voltage level convolution.")
    print("\nWhen using wildcards in the input file name it is important to use small quotes \'example*.bin\'")

if __name__ == "__main__":

    # default values
    input_files = []
    sfreq = 0.0
    cfreq = 0.0
    plot_type = None
    channel = ''
    subchan = 0 # sub channel to plot
    atime = 0
    start_sample = 0L
    stop_sample  = -1L
    modulus = None
    integration = 1

    zscale = (0,0)

    bins = 256

    title = ""
    log_scale = False
    detrend = False
    show_plots = True
    plot_file = ""

    msl_code_length = 0
    msl_baud_length = 0


    # parse the command line arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hi:p:c:o:r:b:z:m:t:a:ld:s:')
    except:
        usage()
        sys.exit()

    for opt, val in opts:

        if opt in ('-h'):
            usage()
            sys.exit()
        elif opt in ('-i'):
            input_files.append(val)
        elif opt in ('-s'):
            plot_file = val
        elif opt in ('-p'):
            plot_type = val
        elif opt in ('-c'):
            sp = string.split(val,':')
            if len(sp) > 0:
                channel = sp[0]
            if len(sp) > 1:
                subchan = int(sp[1])
        elif opt in ('-o'):
            cfreq = float(val)
        elif opt in ('-a'):
            tuple_time = time.strptime(val, "%Y-%m-%dT%H:%M:%S")
            print tuple_time
            atime = calendar.timegm(tuple_time)
            print atime

        elif opt in ('-r'):
            sp = string.split(val,':')

            if len(sp) > 0:
                start_sample = long(sp[0])

            if len(sp) > 1:
                stop_sample = long(sp[1])

            if len(sp) > 2:
                modulus = long(sp[2])

            if len(sp) > 3:
                integration = long(sp[3])

            if len(sp) > 4:
                print("Unknown range format.")
                usage()
                sys.exit()

        elif opt in ('-z'):
            zl,zh = string.split(val,':')
            zscale = (float(zl), float(zh))

        elif opt in ('-b'):
            bins = int(val)
        elif opt in ('-t'):
            title = val
        elif opt in ('-l'):
            log_scale = True
        elif opt in ('-d'):
            detrend = True
        elif opt in ('-m'):
            cl,bl = string.split(val,':')
            msl_code_length = int(cl)
            msl_baud_length = int(bl)

    for f in input_files:
        print("file %s" % f)

        try:
            print "loading data"

            drf = digital_rf_hdf5.read_hdf5(f)

            chans = drf.get_channels()
            if channel == '':
                chidx = 0
            else:
                chidx = chans.index(channel)

            ustart, ustop = drf.get_bounds(chans[chidx])
            print ustart, ustop

            print "loading metadata"
	    
            mdf = dmd.read_digital_metadata(f + '/' + chans[chidx] + '/metadata')
            mdstart = mdf.get_bounds()

            # Note this could be more intelligent. Assumes no changes.
	    #mdt = mdf.read_latest()
            # print mdt
            #md = mdt[mdt.keys()[0]]
            sfreq = drf.get_rf_file_metadata(chans[chidx])["samples_per_second"].tolist()[0]
            toffset = start_sample

            print toffset

            if atime == 0:
                atime = ustart
            else:
                atime = int(atime * sfreq)

            sstart = atime + long(toffset)
            dlen = stop_sample - start_sample + 1

            print sstart, dlen
            md_dict = mdf.read(mdstart[0],mdstart[1])
            #print md_dict
            md_df = pd.DataFrame.from_dict(md_dict, orient='index')
            #print md_df
            df_idx = md_df.index
            df_loc = df_idx.get_loc(sstart,method='ffill')
            metadata = md_df.iloc[df_loc]
            #print metadata
            cfreq = metadata['center_frequencies'][0]

            d = drf.read_vector(sstart, dlen, chans[chidx])

            d = d[:,subchan]
            print d.shape

            print "d", d[0:10]
            #edata = EchotekBinaryData.EchotekBinaryData(f, frame_mode, frame_size, start_sample, stop_sample, verify=False)
            #d = edata.getComplexData()

            if len(d) < (stop_sample - start_sample):
                print "Probable end of file, the data size is less than expected value."
                sys.exit()

            if msl_code_length > 0:
                d = apply_msl_filter(d, msl_code_length, msl_baud_length)

        except:
            print("problem loading file %s" % f)
            traceback.print_exc(file=sys.stdout)
            sys.exit()

        save_plot = False
        if plot_file != "":
            path_head,_ = os.path.split(plot_file)
            if not os.path.isdir(path_head):
                os.makedirs(path_head)
            save_plot = True
        print "generating plot"

        try:
            if plot_type == 'power':
                power_process(d, sfreq, toffset, modulus, integration, log_scale, zscale, title)

            elif plot_type == 'iq':
                iq_process(d, sfreq, toffset, modulus, integration, log_scale, title)

            elif plot_type == 'phase':
                phase_process(d, sfreq, toffset, modulus, integration, log_scale, title)

            elif plot_type == 'voltage':
                voltage_process(d, sfreq,  toffset, modulus, integration, log_scale, title)

            elif plot_type == 'histogram':
                histogram_process(d, sfreq, toffset, modulus, integration, bins, log_scale, title)

            elif plot_type == 'spectrum':
                spectrum_process(d, sfreq, cfreq, toffset, modulus, integration, bins, log_scale, zscale, detrend, title, 'b')

            elif plot_type == 'specgram':
                specgram_process(d, sfreq, cfreq, toffset, modulus, integration, bins, detrend, log_scale, zscale, title)

            elif plot_type == 'rti':
                rti_process(d, sfreq, toffset, modulus, integration, detrend, log_scale, zscale, title)

            elif plot_type == 'sti':
                sti_process(d, sfreq, cfreq, toffset, modulus, integration, bins, detrend, log_scale, zscale, title)

            else:
                print("Unknown plot type %s" % plot_type)
                save_plot = False
            
            if save_plot:
                print 'saving plot'
                matplotlib.pylab.savefig(plot_file)
        except:
            traceback.print_exc(file=sys.stdout)
            sys.exit()
