"""The beacon_plot module is used to create plots using data that is normally
introduced from the beacon receiver HDF5 file.

$Id: BeaconPlot.py 142 2007-08-10 19:54:09Z damiana $
"""

import math
import datetime
import string

import matplotlib.pylab
import matplotlib.toolkits.basemap
import matplotlib.font_manager

import numpy

# use beacon_analysis to get phases for mag/phase plot
# (path is relative to the path of the driver script which is
# currently in src/control)
import sys
sys.path.append("../modules/highlevel")
import beacon_analysis

class BeaconPlot:
    """The class BeaconPlot is used to create engineering and science plots
    for the software radio beacon receiver. Several different plots are
    possible, as listed in the usage example.

    Usage example:
        import BeaconPlot
        
        bp = BeaconPlot.BeaconPlot()
        bp.timePlot([vS4List, uS4List], interval, ['VHF','UHF'],
                    'S4', startTime, 'COSMIC 3', 's4plot.png')
        bp.azElPlot(azVals, elVals, interval,
                    startTime, 'COSMIC 3', 'azel.png')
        bp.mapPlot(eLats, eLons, fLats, fLons, 42.6, -71.5,
                   1000000000, startTime, 'COSMIC 3', tecVals,
                   'TEC', 'map.png', units='e- / m^2')
        bp.magPhasePlot(uhfVoltages, 10000000, uhfFreqs, ranges,
                        startTime, 'COSMIC 3', 'UHF', 'magphase.png')
        bp.voltageSpectrogram(uhfVoltages, 10000000, uhfFreqs, ranges,
                              startTime, 'COSMIC 3', 'UHF', 'spec.png')

    Public attributes:
        dpi: controls the output resolution of images
        thumbDpi: controls the output resolution of thumbnails

    Non-standard Python modules used: numpy, matplotlib,
                                      matplotlib.toolkits.basemap

    Written by "Damian Ancukiewicz":
    mailto:damiana@haystack.mit.edu Jul 11, 2007
    """

    def __init__(self, dpi=100, thumbDpi=20):
        """__init__ is a public method that initializes the beaconPlot object.

        Inputs:
            dpi: Sets the resolution of all output images.
            thumbDpi: Sets the resolution of the output thumbnails

        Returns: none

        Affects: dpi, thumbDpi, __nan

        Exceptions: none
        """
        # allow NaN values to be used/returned in methods
        self.__nan = 1e300*1e300-1e300*1e300
        # set dpi
        self.dpi = dpi
        self.thumbDpi = thumbDpi

    def __isNan(self, num):
        """is_nan checks the inputted value to see if it contains
        the IEEE-754 standard NaN value. If it does, True
        is returned; otherwise, False is returned. Algorithm by
        Bill Rideout.
        """
    
        if ( (num < 10.0) == False and (num > -10.0) == False ):
            return True
        else:
            return False

    def timePlot(self, dataArray, interval, dataNameList,
                 valueName, startTime, satName, outFilename,
                 axisType='linear',units='', yLimits=(), legend=True):
        """timePlot is a public method that creates a plot of a value
           vs.  time for a satellite pass. It can be used for both S4
           and sigma-phi.  The vertical axis can be eiter linear or
           logarithmic.

        Inputs:    
            dataArray: list or array of values spaced at even time
                        intervals; must be two-dimensional; size of
                        first dimension determines how many lines are
                        drawn; only up to 5 lines are allowed              
            interval: (ns) integer representing the interval between
                      successive values      
            dataNameList: list of names (and units) of data that is
                          being plotted; must have as many elements as
                          size of first dimension of dataArray
            valueName: string with the name of the value that will be
                       shown in the title and y-axis
            startTime: datetime object representing the time of the
                       first measurement
            satName: string with the name of satellite
            outFilename: string with the filename to which the PNG
                         image of the graph will be written
            axisType: either 'linear' or 'log'; determines the type of
                      vertical axis used
            units: String represeinting the units used for the y-axis
            yLimits: a two-tuple of floats containing the minimum
                     and maximum values for the y-axis; will scale
                     automatically if not provided
            legend: boolean type that determines whether the legend is
                    shown or not
         
         Returns: nothing
         
         Affects: nothing
         
         Exceptions: throws SyntaxError if axisType is formated
                     incorrectly; throws ValueError if there is a
                     mismatch between the first dimension of dataArray
                     and the length of dataNameList or if the size
                     than 5
        """
        
        # see if we are dealing with a linear or logarithmic plog
        if (axisType == 'linear'):
            log = False
        elif (axisType == 'log'):
            log = True
        else:
            raise SyntaxError,'Axis type of "%s" is incorrect' % (axisType)
        
        if (len(dataArray) != len(dataNameList)):
            raise ValueError, 'Mismatch between size of dataArray and dataNameList'
        if (len(dataArray) > 5):
            raise ValueError, 'Size of dataArray is greater than 5'
        
        # make an array of seconds passed since beginning of pass for
        # each value
        timeArray = numpy.arange(0,len(dataArray[0])*interval/1000000000.0,
                                 interval/1000000000.0)
        # figure out end time
        endTime = startTime + datetime.timedelta(microseconds = len(dataArray[0])*interval/1000)
        # create array of logs of values if axis is logarithmic
        if (log == True):
            logArray = numpy.zeros( (len(dataArray),len(dataArray[0])) )
            for i in range(len(dataArray)):
                for j in range(len(dataArray[0])):
                    logArray[i][j] = math.log(dataArray[i][j])
            # replace values with their logs
            dataArray = logArray
            
        # list of colors to choose from
        colors = 'brgky'
        # create plot
        matplotlib.pylab.clf()
        for i in range(len(dataArray)):
            matplotlib.pylab.plot(timeArray,dataArray[i],colors[i]+'-')
        if (len(yLimits) == 2):
            matplotlib.pylab.ylim(yLimits)
        matplotlib.pylab.xlabel('time since beginning of pass (s)')
        ylbl = ''
        if (log == True):
            ylbl += 'log '+valueName
        else:
            ylbl += valueName
        if (len(units) != 0):
            ylbl += ' ('+units+')'
        matplotlib.pylab.ylabel(ylbl)
        matplotlib.pylab.title(valueName+' vs. time: %s\n%s to %s' % 
                            (satName, startTime.strftime("%Y-%m-%d %H:%M:%S"), 
                            endTime.strftime("%Y-%m-%d %H:%M:%S")))
        if (legend == True):
            matplotlib.pylab.legend(dataNameList)
        matplotlib.pylab.grid(True)
        # save plot
        matplotlib.pylab.savefig(outFilename, dpi=self.dpi)
        # save thumbnail
        parts = outFilename.split('.')
        thumbFilename = string.join(parts[0:-1],'.')+'-small.'+parts[-1]
        matplotlib.pylab.savefig(thumbFilename, dpi=self.thumbDpi)
    
    def azElPlot(self, satAzList, satElList, interval,
              startTime, satName, outFilename):
        """azElPlot is a public method that creates a polar plot 
           of the azimuth and elevation of a satellite for a pass.
           
           Inputs:
               satAzList: list/array of the azimuth of the satellite (deg) at
                          regularly spaced intervals over time.
               satElList: list/array of the elevation of the satellite (deg) at
                          regularly spaced intervals over time
               interval: (ns) an integer representing the time interval
                         between samples
               startTime: a datetime object representing the start time
                          of the pass
               satName:  string representing the satellite's name
               outFilename: string with the output filename of the PNG
                            image of the plot
           
           Returns: nothing
           
           Affects: nothing
           
           Exceptions: none
        """
        # figure out end time
        endTime = startTime + datetime.timedelta(microseconds = len(satAzList)*interval/1000)
        # make angles: rotate by 90 degrees to make 0 degrees straight up
        theta = []
        for i in satAzList:
            theta.append(-math.radians(i) + math.pi/2.0)
        # make radii in the proper range (90 degrees: r=0; 0 degrees: r=90)  
        r = []
        for i in satElList:
            r.append(90.0 - i)
        # make plot
        matplotlib.pylab.clf()
        matplotlib.pylab.polar(theta,r,'b-')
        #  draw a point with zero area (invisible) and r=90
        # to set scale properly (dirty, but effective)
        matplotlib.pylab.scatter([math.pi/2.0],[90],s=0)
        # draw beginning (green) and end (red) points
        matplotlib.pylab.scatter(theta[0:1],r[0:1],s=20, c='g', zorder=10)
        matplotlib.pylab.scatter([theta[-1]],[r[-1]], s=20, c='r', zorder=10)
        # make azimuth and elevation labels
        thetalabels = [str(i % 360) for i in range(90,-270,-45)]
        matplotlib.pylab.thetagrids(range(0,360,45),
                                    thetalabels, frac=1.06, size=10)
        matplotlib.pylab.title('Azimuth and elevation relative to receiver: %s\n%s to %s' % 
                           (satName, startTime.strftime("%Y-%m-%d %H:%M:%S"), 
                            endTime.strftime("%Y-%m-%d %H:%M:%S")), size=12)
        rlabels = [str(i) for i in range(80,-10,-10)]
        matplotlib.pylab.rgrids(range(10,100,10),rlabels, size=10)
        # save plot 
        matplotlib.pylab.savefig(outFilename, dpi=self.dpi)
        # save thumbnail
        parts = outFilename.split('.')
        thumbFilename = string.join(parts[0:-1],'.')+'-small.'+parts[-1]
        matplotlib.pylab.savefig(thumbFilename, dpi=self.thumbDpi)

    def mapPlot(self, eLatList, eLonList, fLatList, fLonList,
                recLat, recLon, interval, startTime, satName,
                valueList, valueName, outFilename, units=''):
        
        """ mapPlot is a public method that draws the intersections
            between the line of sight from the receiver to the
            satellite and the E and F regions of the
            ionosphere. Additionally, it takes in a list of values
            that it then draws on the F region, with different colors
            corresponding to different values. Thus, it can be used to
            plot S4, sigma phi, TEC, etc. over distance.

            Inputs:
                eLatList, eLonList  - lists/arrays of the latitudes and
                                      longitudes of the E region intersections
                fLatList, fLonList - lists/arrays of the latitudes and
                                     longitudes of the F region intersections
                recLat, recLon - floats representing the latitude and
                                 longitude of the receiver
                interval - (ns) an integer representing the interval
                           between successive data/ephemeris values
                startTime - a datetime object representing the start
                            time of the pass
                satName - a string with the name of the satellite
                valueList - a list containing the values that will
                            be plotted by color on the F region line
                valueName - a string with the name of the vlaue plotted
                outFilename - a string containing the file to write to
                units - a string with the units used for plotting

                All lists must have the same length.

            Returns: nothing

            Affects: nothing

            Exceptions: none
            
            (note: would be good to show legend as a picture)
        """
        # figure out where the corners of the map should be
        extremeLats = [min(eLatList), max(eLatList), min(fLatList),
                       max(fLatList), recLat]
        extremeLons = [min(eLonList), max(eLonList), min(fLonList),
                       max(fLonList), recLon]
        minLat, maxLat = min(extremeLats), max(extremeLats)
        minLon, maxLon = min(extremeLons), max(extremeLons)

        latRange = maxLat - minLat
        lonRange = maxLon - minLon
        maxRange = max([latRange, lonRange])
        Dlat = (maxRange - latRange) / 2.
        Dlon = (maxRange - lonRange) / 2.

        # clear plot if there is any
        matplotlib.pylab.clf()
        # will be overwritten by the map; is used to later draw
        # a colorbar
        matplotlib.pylab.scatter(eLatList, valueList, c=valueList)
        
        # create a new map
        map = matplotlib.toolkits.basemap.Basemap(projection='mill',
                                                  lat_0 = recLat,
                                                  lat_ts = recLat,
                                                  llcrnrlat = minLat-Dlat-2.,
                                                  llcrnrlon = minLon-Dlon-2.,
                                                  urcrnrlat = maxLat+Dlat+2.,
                                                  urcrnrlon = maxLon+Dlon+2.,
                                                  lon_0=recLon,
                                                  resolution='i',
                                                  area_thresh=500.)
        map.drawcoastlines(color='0.75')
        map.drawcountries(color='0.75')
        map.fillcontinents(color='0.98')
        map.drawmapboundary()

        # figure out how to draw parallels and meridians
        # find scale (multiple of 5 degrees)
        mapScale = round(maxRange / (4.*5.))*5.
        if (mapScale == 0.):
            mapScale = 2.

        # find point nearest center of map that is on the grid
        nearestLat = round((minLat+maxLat)/(2.*mapScale))*mapScale
        nearestLon = round((minLon+maxLon)/(2.*mapScale))*mapScale

        # draw the parallels/meridians
        map.drawmeridians(numpy.arange(nearestLon-2*mapScale,
                                       nearestLon+3*mapScale,mapScale),
                          labels=[0,0,0,1],linewidth=0.5,size=10)
        map.drawparallels(numpy.arange(nearestLat-2*mapScale,
                                       nearestLat+3*mapScale,mapScale),
                          labels=[0,1,0,0],linewidth=0.5,size=10)

        # turn latitude/longitude into map's X-Y coordinates
        exList,eyList = map(eLonList,eLatList)
        fxList,fyList = map(fLonList,fLatList)
        recX,recY = map(recLon,recLat)

        # make sure to not draw points on F region where there are
        # missing values
        for i in range(len(fxList)):
            if self.__isNan(valueList[i]):
                fxList[i], fyList[i] = self.__nan, self.__nan

        # draw points
        recPoint = map.plot([recX],[recY], 'b+')
        eLine = map.plot(exList,eyList, zorder=10)
        map.scatter(fxList,fyList, c=valueList, linewidth=0, zorder=10, s=5)

        endTime = startTime + datetime.timedelta(microseconds = len(eLatList)*interval/1000)
        matplotlib.pylab.title('E & F region penetration, %s spatial distribution: %s\n%s to %s' % 
               (valueName, satName, startTime.strftime("%Y-%m-%d %H:%M:%S"), 
                endTime.strftime("%Y-%m-%d %H:%M:%S")), size=12)

        # draw the colorbar (based off of the now-invisible scatter plot,
        # which has the same values)
        matplotlib.pylab.colorbar(pad=0.11)

        unitsX, unitsY = map(maxLon+Dlon+9.5,(minLat+maxLat)/2.)
        matplotlib.pylab.text(unitsX,unitsY,units)
        matplotlib.pylab.legend((recPoint,eLine),
                                ('Receiver', 'E region'),
                                loc=(-0.35,0.9),
                         prop=matplotlib.font_manager.FontProperties(size=10))
        # save map
        matplotlib.pylab.savefig(outFilename, dpi=self.dpi)
        # save thumbnail
        parts = outFilename.split('.')
        thumbFilename = string.join(parts[0:-1],'.')+'-small.'+parts[-1]
        matplotlib.pylab.savefig(thumbFilename, dpi=self.thumbDpi)

    def magPhasePlot(self, voltages, interval, frequencies,
                     ranges, startTime, satName, freqStr, outFilename):
        """ magPhasePlot is a public method that plots the magnitude and phase
            of a complex voltage signal, one above the other.

            Inputs:
                voltages: (V) list or array of complex voltages
                interval:  (ns) integer representing the interval
                           between successive voltage values
                frequencies: (Hz) the frequencies at which the voltages
                             were found at (same length as voltages)
                ranges: (m) the ranges from the satellite to receiver
                             (same length as voltages)
                startTime: datetime object containing the start time of pass
                satName: string with the satellite name
                freqStr: a string describing the frequency at which
                              the signal was found at (ex. 'VHF')
                outFilename: string of the fiename to which to write the file
                
            Returns: nothing

            Affects: nothing

            Exceptions: none
        """
        # find magniutde 
        magnitudes = numpy.zeros(len(voltages))
        for i in range(len(voltages)):
            magnitudes[i] = abs(voltages[i])
        # find phases
        phases = beacon_analysis.find_phases(voltages, frequencies, ranges)
        phases = beacon_analysis.detrend_phases(phases)

        # normalize phases to -pi..pi 
        for i in range(len(phases)):
            phases[i] = (phases[i] + math.pi) % (2*math.pi) - math.pi
        
        # make an array of seconds passed since beginning of pass for
        # each value
        timeArray = numpy.arange(0,len(voltages)*interval/1000000000.0,
                                 interval/1000000000.0)
        # figure out end time
        endTime = startTime + datetime.timedelta(microseconds = len(voltages)*interval/1000)

        # draw plot
        matplotlib.pylab.clf()
        matplotlib.pylab.subplot(211)
        matplotlib.pylab.plot(timeArray,magnitudes)
        matplotlib.pylab.ylim(ymin=0)
        matplotlib.pylab.title('Voltage magnitude and phase (%s): %s\n%s to %s' % 
               (freqStr, satName, startTime.strftime("%Y-%m-%d %H:%M:%S"), 
                endTime.strftime("%Y-%m-%d %H:%M:%S")), size=12)
        matplotlib.pylab.ylabel('magnitude (V)')
        
        matplotlib.pylab.subplot(212)
        matplotlib.pylab.plot(timeArray,phases)
        # draw line at y=0
        matplotlib.pylab.plot(timeArray,numpy.zeros(len(voltages)),'k:')
        matplotlib.pylab.ylim((-math.pi,math.pi))
        matplotlib.pylab.ylabel('phase (rad)')
        matplotlib.pylab.xlabel('time since beginning of pass (s)')

        # save plot
        matplotlib.pylab.savefig(outFilename, dpi=self.dpi)
        # save thumbnail
        parts = outFilename.split('.')
        thumbFilename = string.join(parts[0:-1],'.')+'-small.'+parts[-1]
        matplotlib.pylab.savefig(thumbFilename, dpi=self.thumbDpi)

    def voltageSpectrogram(self, voltageList, interval, frequencies,
                           ranges, startTime, satName, freqStr, outFilename):
        """ powerSpectrogram draws a spectrogram of a list of complex
            signal voltages.

            Inputs:
                voltageList: (V) a list of complex voltages
                interval: (ns) an integer representing the sampling
                          interval
                frequencies: (Hz) a list of floats representing the
                             frequencies at which the signal was found
                ranges:   (m) a list of floats with the ranges to the
                              satellite
                startTime: a datetime object representing the
                           start time of the pass
                satName: a string with the name of the satellite
                freqStr: a string describing the frequency at which
                              the signal was found at (ex. 'VHF')
                outFilename: a string containing the name of the image
                             file to write to
                             
            Returns: nothing

            Affects: nothing

            Exceptions: none
            """
        # find phase, detrend phase, find magnitude, and reconstruct
        # voltages from those to compensate for difference in range
        magnitudes = numpy.zeros(len(voltageList))
        for i in range(len(voltageList)):
            magnitudes[i] = abs(voltageList[i])
        phases = beacon_analysis.find_phases(voltageList,frequencies,ranges)
        phases = beacon_analysis.detrend_phases(phases)
        newVoltages = numpy.zeros(len(voltageList), dtype='complex')
        for i in range(len(voltageList)):
            newVoltages[i] = magnitudes[i]*(math.cos(phases[i])+1j*math.sin(phases[i]))
            
        # figure out sampling frequency
        samplingFreq = 1e9 / interval

        # figure out end time
        endTime = startTime + datetime.timedelta(microseconds = len(voltageList)*interval/1000)

        # draw spectrogram
        matplotlib.pylab.clf()
        matplotlib.pylab.specgram(newVoltages, Fs=samplingFreq)
        matplotlib.pylab.colorbar()
        matplotlib.pylab.xlabel('time since beginning of pass (s)')
        matplotlib.pylab.ylabel('frequency (Hz)')
        matplotlib.pylab.title('Power spectral density (%s): %s\n%s to %s' % 
               (freqStr, satName, startTime.strftime("%Y-%m-%d %H:%M:%S"), 
                endTime.strftime("%Y-%m-%d %H:%M:%S")), size=12)
    
        # save plot
        matplotlib.pylab.savefig(outFilename, dpi=self.dpi)
        # save thumbnail
        parts = outFilename.split('.')
        thumbFilename = string.join(parts[0:-1],'.')+'-small.'+parts[-1]
        matplotlib.pylab.savefig(thumbFilename, dpi=self.thumbDpi)
                                                    

            
        

    
