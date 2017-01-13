#!/usr/bin/env python
"""@package predictpasses

  Read in beacon.tle and predict the passes of all satellites in that
  file. Create passes.txt, which contains transit times, center frequencies
  peak elevations and satellite names. The file beacon.tle can be generated using
  get_beacontle.py, which downloads the necessary files from celestrak combines
  TLE lines with satellite frequency information.

  Passes.txt is then read by the beacon satellite recorder program, to determine
  what and when to record.

  We rely on the proven pyephem package for ephemeris calculations.

  Also write a flight-000???.pass file, which contains an az,el,range
  timeseries of the pass. 
                            date                 directory                    lat        long       elev   sat_elev_cutoff
  Usage: ./predictpasses.py "2010/4/30 10:00:00" /data/beacon/data/2010.04.30 67:21:51.0 26:37:45.0 197.03 10

  (c) 2009 Juha Vierinen
  
"""
import ephem, sys, math
import os, errno, time 

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise

def make_dir(ddir,startd):
    # the start date in pyephem format: 
    dtuple = startd.tuple()
    dirname = "%s/%04d.%02d.%02d" % (ddir,dtuple[0],dtuple[1],dtuple[2])
    mkdir_p(dirname)


#
# encapsulate all beacon.tle file operations in this
#
class satellitedata:
    
    # load satellite ephemeris
    def __init__(self, fname, ddir, startd, obs, elcutoff):
        self.ut0 = ephem.Date('1970/1/1 00:00:00')
        self.fname = fname
        self.obs = obs
        self.startd = startd
        self.ddir = ddir
        self.readData()
        self.elcutoff = elcutoff
        self.getPasses(elcutoff)
        
        
    def readData(self):
        f = open(self.fname, 'r')
        data = f.readlines()
        f.close()
        
        self.nsatellites = int(math.floor(len(data)/4))
        self.name = []
        self.tle2 = []
        self.tle3 = []
        self.ppm = []

        print "number of satellites: ",self.nsatellites," starting search at ",self.startd
        for i in range(0,self.nsatellites):

            self.name.append(data[i*4].strip())
            self.tle2.append(data[i*4+1].strip())
            self.tle3.append(data[i*4+2].strip())
            self.ppm.append(float(data[i*4+3].strip()))

            print i," ",self.name[i]," ",self.ppm[i]," ppm"
        print "done"

    def getPasses(self,elcutoff):

        self.risings = []
        self.risingsUTC = []
        self.settings = []
        self.peakel = []
        f = 0
        passidx = 0
        for si in range(0,self.nsatellites):
            sat = ephem.readtle(self.name[si],self.tle2[si],self.tle3[si])
            self.risings.append([])
            self.risingsUTC.append([])
            self.settings.append([])
            self.peakel.append([])
            satUp = 0
            peak = 0
            
            # look a little over one day in advance
            for i in range(1,int(60*60*24.5),10):
                self.obs.date = self.startd
                self.obs.date = self.obs.date + float(i)*ephem.second
                sat.compute(self.obs)
                el = 360.0*float(sat.alt)/2.0/math.pi
                az = 360.0*float(-sat.az)/2.0/math.pi + 90
                r = float(sat.range)
                alt = float(sat.elevation)
                vel = float(sat.range_velocity)
                sublat = 360.0*float(sat.sublat)/2.0/math.pi
                sublong = 360.0*float(sat.sublong)/2.0/math.pi
                
                # make sure that we only log satellites that go up after t0, not ones that already are up.
                if el > elcutoff and i < 60.0*60.0*24.0 and satUp == 0:
                    f = file("%s/flight-%06d.pass"%(self.ddir,passidx),"w")
                    risetime = ((self.obs.date - self.ut0)*3600.0*24.0 + i*ephem.second)
                    print passidx," rise ",self.name[si]," ",risetime
#                    print time.ctime(risetime)
                    self.risings[si].append(risetime)
                    self.risingsUTC[si].append(time.ctime(risetime))
                    satUp = 1
                    
                if el < elcutoff and satUp == 1:
                    settime = ((self.obs.date - self.ut0)*3600.0*24.0 + i*ephem.second)
                    self.settings[si].append(settime)
                    print "set  ",self.name[si]," ",settime," peak elevation ",peak,"\n"
                    f.close()
                    satUp = 0
                    self.peakel[si].append(peak)
                    peak = 0
                    passidx = passidx + 1
                if satUp == 1:
                    f.write("%1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f\n"%( ((self.obs.date - self.ut0)*3600.0*24.0 + i*ephem.second),el,az,r,vel,alt,sublat, sublong))
                if el > peak and satUp == 1:
                    peak = el

    ## write a file containing passes.
    def writeFile(self,fname):
        f = open(fname, 'w')

        # b=time.strftime("%Y.%m.%d_%H.%M.%S",time.strptime(time.ctime(a)))
        # stmp="%s.%09d" % (b,int((a-math.floor(a))*1000000000))
        f.write("# start_time end_time freq1 (MHz) freq2 (MHz) bw1 (kHz) bw2 (kHz) peak_elevation name\n")
        for si in range(0,self.nsatellites):
            for pi in range(0,min(len(self.risings[si]),len(self.settings[si]))):
                f.write("%1.10f %1.10f %1.6f %1.6f %1.3f %1.3f %1.3f %s\n" %
                        (self.risings[si][pi],
                         self.settings[si][pi],
                         150 + 150*self.ppm[si]/1e6,
                         400 + 400*self.ppm[si]/1e6,
                         40, 40, self.peakel[si][pi],
                         self.name[si]))
        f.close()

#                            date                 directory                    lat        long       elev   sat_elev_cutoff
#  Usage: ./predictpasses.py "2010/4/30 10:00:00" /data/beacon/data/2010.04.30 67:21:51.0 26:37:45.0 197.03 10

def write_pass_files(datestr="2010/4/30 10:00:00",
                     data_dir="/data/beacon/data/2010.04.30",
                     latitude="67:21:51.0",
                     longitude="26:37:45.0",
                     altitude=197.03,
                     sat_elev_cutoff=10.0):

    # location of the receiver station
    obs = ephem.Observer()
    startd = ephem.Date(datestr)
    ddir = data_dir
    obs.lat = latitude
    obs.long = longitude
    obs.elevation = altitude
    elevationCutoff = int(sat_elev_cutoff)
    
    mkdir_p(ddir)
    sats = satellitedata("%s/beacon.tle"%(ddir),ddir,startd,obs,elevationCutoff)
    sats.writeFile("%s/passes.txt"%(ddir))
