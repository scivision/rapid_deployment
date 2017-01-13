#!/usr/bin/env python
#
# pythonic phase difference curve
#
# Read 150 and 400 MHz signals
#
# resample 150 MHz signal to 15 kHz
#
#
# Use this to determine peak frequency
#
#
import matplotlib as mpl
mpl.use('Agg')
import os
import sys
import numpy
import scipy.constants
import scipy.interpolate
import scipy.linalg
import csv
import matplotlib.pyplot as plt
import stuffr
import h5py
import glob
from mpl_toolkits.basemap import Basemap

import beacon_conf as c

#numpy.seterr(all='raise')

debug=True
verbose_debug=True
debug_plot = False

## read list of recorded flights
def read_flights(dname):
    f = open("%s/flights.log"%(dname))
    c = csv.reader(f,delimiter=" ")
    results = []
    idx = []
    for row in c:
        results.append({"name":row[12],
                        "idx":int(row[0]),
                        "t0":float(row[1]),
                        "t1":float(row[2]),
                        "fname0":row[3],
                        "fname1":row[4],
                        "freq0":float(row[5]),
                        "freq1":float(row[6]),
                        "dec0":int(row[7]),
                        "dec1":int(row[8]),
                        "sr0":float(row[9]),
                        "sr1":float(row[10])})
        idx.append(int(row[0]))
    f.close()
    return({"idx":numpy.array(idx,dtype=numpy.int),"table":results})

## find correct row
def get_flight_info(dname,i):
    f = read_flights(dname)
    idx = numpy.where(f["idx"] == i)[0]
    return(f["table"][idx])

## read ephemeris file
def read_ephem(fname):
    if debug:
        print("reading %s"%(fname))
    f = numpy.genfromtxt(fname)
    t = f[:,0]
    el = f[:,1]
    az = f[:,2]
    r = f[:,3]
    vel2 = numpy.diff(f[:,3])/numpy.diff(t)
    vel = f[:,4]
    alt= f[:,5]
    sublat= f[:,6]
    sublon= f[:,7]
    t[len(t)-1]=t[len(t)-1]+600
    t[len(t)-2]=t[len(t)-1]+500
    tdop = (t[0:(len(t)-1)]+t[1:len(t)])/2.0
    tdop[0] = tdop[0]-600.0
    return({"t":t,
            "el":scipy.interpolate.interp1d(t,el),
            "az":scipy.interpolate.interp1d(t,az),
            "r":scipy.interpolate.interp1d(t,r,kind="cubic"),
            "vel":scipy.interpolate.interp1d(t,vel,kind="cubic"),
            "vel2":scipy.interpolate.interp1d(tdop,vel2,kind="cubic"),
            "alt":scipy.interpolate.interp1d(t,alt),
            "sublat":scipy.interpolate.interp1d(t,sublat),
            "sublon":scipy.interpolate.interp1d(t,sublon),
            "ell":el,
            "azl":az,
            "rl":r,
            "altl":alt,
            "vell":vel,
            "sublatl":sublat,
            "sublonl":sublon})

## read datafile
def read_data(dname,flight=42):
    f = get_flight_info(dname,flight)
    e = read_ephem("%s/flight-%06d.pass"%(dname,flight))
    z0 = numpy.fromfile("%s/flight-%06d.001.bin"%(dname,flight),dtype=numpy.complex64)
    z1 = numpy.fromfile("%s/flight-%06d.002.bin"%(dname,flight),dtype=numpy.complex64)
    return({"f":f,"e":e,"z0":z0,"z1":z1,"dname":dname,"flight":flight})

def outlier_removed_fit(m,n_iter=10,polyord=7):
    m2 = numpy.copy(m)
    tv = numpy.linspace(-1,1,num=len(m))
    A = numpy.zeros([len(m),polyord])
    for j in range(polyord):
        A[:,j] = tv**(float(j))
    A2 = numpy.copy(A)
    fit = None
    for i in range(n_iter):
        xhat = scipy.linalg.lstsq(A2,m2)[0]
        fit = numpy.dot(A,xhat)
        resid = (fit - m2)
        std = numpy.std(resid)
        bidx = numpy.where(numpy.abs(resid) > 2.0*std)[0]
        for bi in bidx:
            A2[bi,:]=0.0
            m2[bi]=0.0
    if debug_plot:
        plt.plot(m2,label="outlier removed")
        plt.plot(m,label="original")
        plt.plot(fit,label="fit")
        plt.legend()
        plt.ylim([numpy.min(fit)-std*3.0,numpy.max(fit)+std*3.0])
        plt.show()
    return(fit)

## estimate quicklook spectrum
def get_quicklook_spectrum(d,window=4096*2,n_window=500,bandwidth=100):
    step = int(numpy.floor((len(d["z0"])-2.0*window)/n_window))
    res0 = numpy.zeros([n_window,window],dtype=numpy.complex64)
    res1 = numpy.zeros([n_window,window],dtype=numpy.complex64)

    freqm = numpy.zeros([n_window])

    tvec = numpy.zeros([n_window])
    fvec = numpy.linspace(-d["f"]["sr0"]/2.0,d["f"]["sr0"]/2.0,num=window)

    wfun = stuffr.hanning(window)
    idx = numpy.arange(window)
    sr = d["f"]["sr0"]

    fidx = numpy.arange(-bandwidth/2,bandwidth/2,dtype=numpy.int) + window/2

    for i in range(n_window):
        tvec[i]=float(i*step)/sr + d["f"]["t0"]
        z00 = d["z0"][i*step + idx]
        z01 = d["z0"][i*step + idx + window]
        z10 = d["z1"][i*step + idx]
        z11 = d["z1"][i*step + idx + window]
        doppler0 = 1e6*d["f"]["freq0"]*(d["e"]["vel2"](float(i*step)/sr + d["f"]["t0"]))/scipy.constants.c
        doppler1 = 1e6*d["f"]["freq1"]*(d["e"]["vel2"](float(i*step)/sr + d["f"]["t0"]))/scipy.constants.c

        F0 = numpy.fft.fftshift(numpy.fft.fft(wfun*z00*numpy.exp(1.0j*2.0*numpy.pi*doppler0*(idx/sr))))
        F1 = numpy.fft.fftshift(numpy.fft.fft(wfun*z01*numpy.exp(1.0j*2.0*numpy.pi*doppler0*(idx/sr + float(window)/sr))))
        res0[i,:]=F0*numpy.conj(F1)

        F0 = numpy.fft.fftshift(numpy.fft.fft(wfun*z10*numpy.exp(1.0j*2.0*numpy.pi*doppler1*(idx/sr))))
        F1 = numpy.fft.fftshift(numpy.fft.fft(wfun*z11*numpy.exp(1.0j*2.0*numpy.pi*doppler1*(idx/sr+float(window)/sr))))

        res1[i,:]=F0*numpy.conj(F1)
        freqm[i] = fvec[fidx[numpy.argmax(numpy.abs(res1[i,fidx]))]]

        res0[i,:]=res0[i,:]/numpy.median(numpy.abs(res0[i,:]))
        res1[i,:]=res1[i,:]/numpy.median(numpy.abs(res1[i,:]))

    tvec[0]=tvec[0]-100
    tvec[len(tvec)-1]=tvec[len(tvec)-1]+100

    dopfit = outlier_removed_fit(freqm)
    doppler_residual = scipy.interpolate.interp1d(tvec,dopfit)
    cspec = numpy.mean(numpy.abs(res0*numpy.conj(res1)),axis=0)
    return({"cspec":cspec,"max_bin":numpy.argmax(cspec),
            "doppler_residual":doppler_residual,
            "tvec":tvec,"fvec":fvec,"res1":res1,"res0":res0})

#
# combination of coherent and incoherent integration
#
def estimate_phase_curve(d,fft_window=512*2*2,incoh_int=160,sfactor=4):
    n_points = (len(d["z0"])-2*fft_window)/((fft_window/sfactor)*incoh_int)
    ql = get_quicklook_spectrum(d,window=8192)
    fi = fft_window/2
    #    fi = ql["max_bin"]
    #print(fi)

    wfun = stuffr.hanning(fft_window)

    tvec = numpy.zeros([n_points],dtype=numpy.float64)

    spec0 = numpy.zeros([fft_window],dtype=numpy.float32)
    spec1 = numpy.zeros([fft_window],dtype=numpy.float32)

    phase0 = numpy.zeros([n_points],dtype=numpy.complex64)
    phase1 = numpy.zeros([n_points],dtype=numpy.complex64)

    std0 = numpy.zeros([n_points],dtype=numpy.complex64)
    std1 = numpy.zeros([n_points],dtype=numpy.complex64)

    snr0 = numpy.zeros([n_points],dtype=numpy.float32)
    snr1 = numpy.zeros([n_points],dtype=numpy.float32)

    pwin0 = numpy.zeros([incoh_int],dtype=numpy.complex64)
    pwin1 = numpy.zeros([incoh_int],dtype=numpy.complex64)

    phase_amp_0 = numpy.zeros([incoh_int*n_points],dtype=numpy.complex64)
    phase_amp_1 = numpy.zeros([incoh_int*n_points],dtype=numpy.complex64)

    idx = numpy.arange(fft_window)

    sr = d["f"]["sr0"]
    wi = 0

    phase_00 = numpy.exp(1.0j*0.0)
    phase_01 = numpy.exp(1.0j*0.0)
    phase_10 = numpy.exp(1.0j*0.0)
    phase_11 = numpy.exp(1.0j*0.0)

    for i in range(n_points):
        if debug:
            sys.stdout.write(".")
            sys.stdout.flush()
#            print("%d/%d"%(i,n_points))
        tvec[i] = i*incoh_int*(fft_window/sfactor)/d["f"]["sr0"] + d["f"]["t0"]
        spec0[:]=0.0
        spec1[:]=0.0

        for k in range(incoh_int):
            z00 = d["z0"][wi*fft_window/sfactor + idx]
            z01 = d["z0"][wi*fft_window/sfactor + idx + fft_window/sfactor]
            z10 = d["z1"][wi*fft_window/sfactor + idx]
            z11 = d["z1"][wi*fft_window/sfactor + idx + fft_window/sfactor]
            tdop = float(wi*fft_window/sfactor)/sr + d["f"]["t0"]
            doppler0 = -1.0*(150.0/400.0)*ql["doppler_residual"](tdop) + 1e6*d["f"]["freq0"]*d["e"]["vel2"](tdop)/scipy.constants.c
            doppler1 = -1.0*ql["doppler_residual"](tdop) + 1e6*d["f"]["freq1"]*d["e"]["vel2"](tdop)/scipy.constants.c

            # channel 0
            csin00 = phase_00*numpy.array(numpy.exp(1.0j*2.0*numpy.pi*doppler0*(idx/sr)),dtype=numpy.complex64)
            F = numpy.fft.fftshift(numpy.fft.fft(wfun*z00*csin00))
            spec0 += numpy.abs(F)**2.0
            F0 = F[fi]
            phase_amp_0[wi]=F0

            csin01 = phase_00*numpy.array(numpy.exp(1.0j*2.0*numpy.pi*doppler0*(idx/sr + float(fft_window/sfactor)/sr)),dtype=numpy.complex64)
            phase_00 = phase_00*numpy.exp(1.0j*2.0*numpy.pi*doppler0*((fft_window/sfactor+1.0)/sr))
            F = numpy.fft.fftshift(numpy.fft.fft(wfun*z01*csin01))
            spec0 += numpy.abs(F)**2.0
            F1 = F[fi]

            pwin0[k] = (F0*numpy.conj(F1))#/(numpy.abs(F0)*numpy.abs(F1))
            phase0[i] += pwin0[k]

            # channel 1
            csin10 = phase_10*numpy.array(numpy.exp(1.0j*2.0*numpy.pi*doppler1*(idx/sr)),dtype=numpy.complex64)
            F = numpy.fft.fftshift(numpy.fft.fft(wfun*z10*csin10))
            spec1 += numpy.abs(F)**2.0
            F0 = F[fi]
            phase_amp_1[wi]=F0

            csin11 = phase_10*numpy.array(numpy.exp(1.0j*2.0*numpy.pi*doppler1*(idx/sr + float(fft_window/sfactor)/sr)),dtype=numpy.complex64)
            phase_10 = phase_10*numpy.exp(1.0j*2.0*numpy.pi*doppler1*((fft_window/sfactor+1.0)/sr))
            F = numpy.fft.fftshift(numpy.fft.fft(wfun*z11*csin11))
            spec1 += numpy.abs(F)**2.0
            F1 = F[fi]

            pwin1[k] = (F0*numpy.conj(F1))#/(numpy.abs(F0)*numpy.abs(F1))
            phase1[i]+= pwin1[k]
            wi += 1

        std0[i] = numpy.std(numpy.angle(pwin0))
        std1[i] = numpy.std(numpy.angle(pwin1))
        snr0[i] = spec0[fi]/numpy.median(spec0)
        snr1[i] = spec1[fi]/numpy.median(spec1)

    phasecurve = float(incoh_int)*numpy.cumsum(numpy.angle(phase1)*(d["f"]["freq0"]/d["f"]["freq1"])-numpy.angle(phase0))

    diff_1 = numpy.angle(phase_amp_1[1:(len(phase_amp_1))]*numpy.conj(phase_amp_1[0:(len(phase_amp_1)-1)]))*(d["f"]["freq0"]/d["f"]["freq1"])
    diff_0 = numpy.angle(phase_amp_0[1:(len(phase_amp_0))]*numpy.conj(phase_amp_0[0:(len(phase_amp_0)-1)]))

    phasecurve_amp = stuffr.decimate(numpy.array(-1.0*numpy.cumsum(diff_1 - diff_0),dtype=numpy.float32),dec=sfactor)
    stdcurve = numpy.sqrt(numpy.cumsum(float(sfactor)*incoh_int*(std0**2.0 + std1**2.0)))

    return({"phasecurve":phasecurve,
            "phasecurve_amp":phasecurve_amp,
            "stdcurve":stdcurve,
            "phase_diff0":phase0,"phase_diff1":phase1,
            "std0":std0,"std1":std1,
            "snr0":snr0,"snr1":snr1,
            "t":tvec,
            "t0_amp":d["f"]["t0"],
            "sr_amp":numpy.float(d["f"]["sr0"])/numpy.float(fft_window),
            "ql":ql})

## overflight
def analyze_flight(dname,flight):
    f = get_flight_info(dname,flight)
    dur = f["t1"]-f["t0"]

    if dur > c.min_duration:
        d = read_data(dname=dname,flight=flight)
        d["pc"] = estimate_phase_curve(d)
        save_result(d)
        plot_overview(d)
     # mark as analyzed
        f = file("%s/nanalyzed.log"%(dname),"a")
        f.write("%d\n"%(flight))
        f.close()
    else:
        if debug:
            print("too short flight. not analyzing. deleting binary files.")
            os.system("rm -f %s/flight-%06d.001.bin"%(dname,flight))
            os.system("rm -f %s/flight-%06d.002.bin"%(dname,flight))
    return(d)

 ## use top directory
def analyze_all(maxn=10):
    last_analyzed=""
    dl = glob.glob("%s/????.??.??"%(c.datadir))
    ## make sure that we don't use too much memory
    dl.sort()
    n = 0
    # go through all days
    for d in dl:
        if verbose_debug:
            print("%s"%(d))
        last_analyzed=d
        # go through all days when recorder has been active
        if os.path.exists("%s/flights.log"%(d)):
            f = read_flights(d)
            an_idx = numpy.array([],dtype=numpy.int)
            if not os.path.exists("%s/nanalyzed.log"%(d)):
                os.system("touch %s/nanalyzed.log"%(d))
            tmpf = file("%s/nanalyzed.log"%(d),"r")
            l = tmpf.readlines()
            tmpf.close()
            if(len(l) > 0):
                an_idx = numpy.genfromtxt("%s/nanalyzed.log"%(d))

            for i in f["idx"]:
                # if analyzed, skip
                if not numpy.in1d(numpy.array([i],dtype=numpy.int),an_idx):
                    if verbose_debug:
                        print("%d not analyzed yet "%(i))
                    if os.path.exists("%s/flight-%06d.001.bin"%(d,i)) and os.path.exists("%s/flight-%06d.002.bin"%(d,i)):
                        analyze_flight(d,flight=i)
                        n += 1
                        if n>maxn:
                            return
                    else:
                        if verbose_debug:
                            print("binary files missing %s flight %d. skipping."%(d,i))

        else:
            if debug:
                print("no data")
    return(last_analyzed)





def save_result(d):
    # remove old file from the way if exists
    os.system("echo \"%s/%d\" > %s/last_an.txt"%(d["dname"],d["flight"],c.datadir))

    fname = "%s/%d@%s@%s.h5"%(d["dname"],numpy.int(d["f"]["t0"]),d["f"]["name"],c.station)
    os.system("rm -f %s"%(fname))
    h = h5py.File(fname)

    # phase curve derived from incoherently averaged signal
    h["phasecurve_150MHz_rad"]=d["pc"]["phasecurve"]
    h["phasecurve_std"]=d["pc"]["stdcurve"]
    h["doppler_residual_ch1"]=d["pc"]["ql"]["doppler_residual"](d["pc"]["t"])
    h["t"]=d["pc"]["t"]

    # phase curve derived from amplitude domain signal
    h["phasecurve_amp"]=d["pc"]["phasecurve_amp"]
    h["t0_amp"]=d["pc"]["t0_amp"]
    h["sr_amp"]=d["pc"]["sr_amp"]

    h["freq0_MHz"]=d["f"]["freq0"]
    h["freq1_MHz"]=d["f"]["freq1"]

    h["snr0"]=d["pc"]["snr0"]
    h["snr1"]=d["pc"]["snr1"]

    tle = file("%s/beacon.tle"%(d["dname"]))
    h["tle"]=tle.readlines()
    tle.close()

    h["ephem_t"] = d["e"]["t"]
    h["ephem_el"] = d["e"]["ell"]
    h["ephem_az"] = d["e"]["azl"]
    h["ephem_range_m"] = d["e"]["rl"]
    h["ephem_altitude_m"] = d["e"]["altl"]
    h["ephem_vel_ms"] = d["e"]["vell"]
    h["ephem_sublat_deg"] = d["e"]["sublatl"]
    h["ephem_sublon_deg"] = d["e"]["sublonl"]

    h["station"] = c.station
    h["station_lat"] = c.station_latitude
    h["station_lon"] = c.station_longitude
    h["station_alt"] = c.station_elevation
    h.close()

def plot_overview(d,flim=[-0.5e3,0.5e3]):
    plt.figure(figsize=(15,15))
    plt.subplot(321)
    # Lambert Conformal Conic map.
    m = Basemap()
    t = d["e"]["t"]
    t = t[0:(len(t)-2)]
    tlim = [t[0],t[len(t)-1]]

    m.drawcoastlines(color="gray")
    m.plot(d["e"]["sublon"](t),d["e"]["sublat"](t),linewidth=3)
    m.plot(c.station_longitude, c.station_latitude, "o")
    plt.title("%s %s\n%s %s"%(c.station,
                              d["f"]["name"],
                              stuffr.unix2datestr(d["f"]["t0"]),
                              stuffr.unix2datestr(d["f"]["t1"])))


    plt.subplot(322)
    plt.plot(d["e"]["sublat"](d["pc"]["t"]),d["pc"]["phasecurve"],color="black")
    plt.xlabel("Latitude (deg)")
    plt.ylabel("Phase difference (3/8)*400-150 MHz (rad)")
    plt.title("Phase difference - Sublatitude")

    # quicklook spectra
    plt.subplot(323)
    tvec=d["pc"]["ql"]["tvec"]
    fvec=d["pc"]["ql"]["fvec"]
    res0=d["pc"]["ql"]["res0"]
    res1=d["pc"]["ql"]["res1"]
    plt.pcolormesh(tvec,fvec,numpy.transpose(stuffr.comprz_dB(numpy.abs(res0))))
    plt.plot(tvec,(150.0/400.0)*d["pc"]["ql"]["doppler_residual"](tvec),"k--",label="doppler resid")
    plt.ylim(flim)
    plt.xlabel("Unix time (s)")
    plt.ylabel("Frequency (Hz)")
    plt.title("Power ch0 (dB) %1.2f MHz"%(d["f"]["freq0"]))
    plt.legend()
    plt.colorbar(orientation="horizontal")

    plt.xlim(tlim)

    plt.subplot(324)
    plt.pcolormesh(tvec,fvec,numpy.transpose(stuffr.comprz_dB(numpy.abs(res1))))
    plt.plot(tvec,d["pc"]["ql"]["doppler_residual"](tvec),"k--",label="doppler resid")
    plt.ylim(flim)
    plt.xlabel("Unix time (s)")
    plt.ylabel("Frequency (Hz)")
    plt.title("Power ch1 (dB), %1.2f MHz"%(d["f"]["freq1"]))
    plt.colorbar(orientation="horizontal")
    plt.legend()
    plt.xlim(tlim)

    plt.subplot(325)
    plt.plot(d["pc"]["t"],d["pc"]["phasecurve"],color="black",label="incoh")
    plt.plot(d["pc"]["t"],d["pc"]["phasecurve"]+d["pc"]["stdcurve"],color="gray")
    plt.plot(d["pc"]["t"],d["pc"]["phasecurve"]-d["pc"]["stdcurve"],color="gray")
    plt.plot(float(d["pc"]["t0_amp"])+numpy.arange(float(len(d["pc"]["phasecurve_amp"])))/float(d["pc"]["sr_amp"]),d["pc"]["phasecurve_amp"],label="amp")

    plt.xlabel("Unix time (s)")
    plt.ylabel("Phase (rad)")
    plt.legend()
    plt.title("Phase curve (150 MHz scale")
    plt.xlim(tlim)
    plt.ylim([numpy.min(d["pc"]["phasecurve"]-2.0*d["pc"]["stdcurve"]),numpy.max(d["pc"]["phasecurve"]+2.0*d["pc"]["stdcurve"]) ])

    plt.subplot(326)
    plt.plot(d["pc"]["t"],10.0*numpy.log10(d["pc"]["snr0"]),label="150")
    plt.plot(d["pc"]["t"],10.0*numpy.log10(d["pc"]["snr1"]),label="400")
    plt.xlim(tlim)
    plt.legend()
    plt.title("SNR")
    plt.xlabel("Unix time (s)")
    plt.ylabel("SNR (dB)")

    pngname = "%s/%d@%s@%s.png"%(d["dname"],int(d["f"]["t0"]),d["f"]["name"],c.station)
    plt.savefig(pngname,dpi=72)
    plt.clf()

if __name__ == "__main__":
    analyze_all()

