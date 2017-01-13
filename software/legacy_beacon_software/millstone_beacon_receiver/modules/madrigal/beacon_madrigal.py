"""
The MadrigalLoad module is used to load an HDF5 file, formatted
according to the current beacon receiver HDF5 spec.

The module creates a Madrigal experiment for each HDF5 file given and
automatically creates the experiment header.

$Id: MadrigalLoad.py 142 2007-08-10 19:54:09Z damiana $
"""

import os,sys
import datetime
import math

import tables
import numpy

import madrigal.metadata
import madrigal.cedar
import madrigal.admin
import madrigal._Madrec

# NOTE THIS IS THE DEPRECATED MADRIGAL API
# USE THE MADRIGAL 3.0 API!

# path is relative to the path of the driver script (which currently
# is in src/control)
sys.path.append('../modules/plot')
import beacon_plot

def find_pierce_point(recLat, recLon, satLat, satLon, satAlt, ppAlt):
    """find_pierce_point is a public function that finds the
    latitude and longitude at which a line of sight from a given
    latitude and longitude with a given azimuth and elevation will
    pierce a given altitude.

    Inputs:
        recLat (deg) - latitude of the receiver 
        recLon (deg) - longitude of the receiver
        satLat (deg) - latitude of the satellite 
        satLon (deg) - longitude of the satellite 
        satAlt (m) - altitude of the satellite 
        ppAlt (m) - altitude of the pierce point 

    Returns: a 2-tuple containing the latitude and longitude
             (in degrees) of the pierce point.

    Algorithm: the receiver's and satellite's coordinates are
               converted to Cartesian coordinates based on the
               simplification of the Earth as a sphere.  The
               intersection is then found between a parametrized
               line from the receiver to the satellite and a
               sphere containing the pierce altitude. This is
               converted back to latitude and longitude.

    Affects: nothing

    Exceptions: throws ValueError if altitude is less than 0 m
    
    (note: might want to have Earth modeled as an oblate spheroid;
     should not be too hard to accomplish)
    """
    # check for proper values
    if (ppAlt < 0.0):
        raise ValueError,'Pierce altitude of %s is incorrect' % (ppAlt)
    
    # define radius of Earth in meters
    R = 6372800.
    
    # get spherical coordinates
    recPhi = math.radians(90 - recLat)
    recTheta = math.radians(recLon)
    satR = R + satAlt
    satPhi = math.radians(90 - satLat)
    satTheta = math.radians(satLon)
    
    # convert receiver coordinates to Cartesian
    recX = R * math.sin(recPhi) * math.cos(recTheta)
    recY = R * math.sin(recPhi) * math.sin(recTheta)
    recZ = R * math.cos(recPhi)
    satX = satR * math.sin(satPhi) * math.cos(satTheta)
    satY = satR * math.sin(satPhi) * math.sin(satTheta)
    satZ = satR * math.cos(satPhi)
    
    # make quadratic equation of the form at^2 + bt + c = 0
    a = math.pow(satX-recX,2)+math.pow(satY-recY,2)+math.pow(satZ-recZ,2)
    b = 2*(recX*(satX-recX)+recY*(satY-recY)+recZ*(satZ-recZ))
    c = math.pow(recX,2)+math.pow(recY,2)+math.pow(recZ,2)-math.pow(R+ppAlt,2)
    
    # solve for t - only the positive result has the actual point
    t = (-b + math.sqrt(math.pow(b,2)-4*a*c))/(2*a)
    
    # find the pierce point in Cartesian coords
    ppX = recX + (satX-recX)*t
    ppY = recY + (satY-recY)*t
    ppZ = recZ + (satZ-recZ)*t
    
    # convert back to latitude/longitude
    ppLat = 90 - math.degrees(math.acos(ppZ/(R+ppAlt)))
    ppLon = math.degrees(math.atan2(ppY,ppX))

    # return answer
    return (ppLat,ppLon)
    
class MadrigalLoad:
    """The class MadrigalLoad is used to create a Madrigal experiment
       from an HDF5 file formatted according to the current beacon receiver
       HDF5 spec.

       Usage example:
            import tables
            import MadrigalLoad
            
            mbl = MadrigalLoad.MadrigalLoad(log, 8888, 8888, 'final', 
                               'Scintillation parameters and TEC',
                               'MIDAS-M Beacon Receiver', 100, 250)
            h5file = tables.openFile('receiver-2007081001.h5')
            mbl.load(h5file)
            h5file.close()

       Written by "Damian Ancukiewicz":
       mailto:damiana@haystack.mit.edu Jul 05, 2007

       Requires pytables, numpy, madrigal.* , BeaconPlot
    """

    def __init__(self, log, kinst, kindat, kindatType,
                 kindatDesc, expTitle, eAltitude, fAltitude):
        """__init__ is a public method that creates a
        madrigalBeaconLoad object and initializes all necessary
        attributes including the hard-coded kinst and kindat values.

           Inputs:
               log - a logging object to log status of the loading.

           Returns: nothing

           Affects: The following is initialized:
               __log: Instance of the logger object
               __kinst: Identifies to Madrigal the type of instrument
                        used to create the experiment.
               __kindat: Identifies to Madrigal the algorithm that was
                         used to process the receiver data
               __kindatType: Madrigal file description corresponding
                             to the kindat
               __kindatDesc: Description of the kindat
               __expTitle: title of the experiment
               __eAltitude: altitude chosen for E region penetration
               __fAltitude: altitude chosen for F region penetration
               __madrigalDB: instance of a madrigal.metadata.madrigalDB object

           Exceptions: none
        """


        # replace this all with a config script later
        self.__log = log
        self.__madrigalDB = madrigal.metadata.MadrigalDB()

        # hard-coded; set by convention (completely bogus for now)
        self.__kinst = kinst
        self.__kindat = kindat
        # __kindatType can't have commas
        self.__kindatType = kindatType
        self.__kindatDesc = kindatDesc
        self.__expTitle = expTitle
        # altitudes (in m) set for E & F regions for finding pierce points
        self.__eAltitude = eAltitude
        self.__fAltitude = fAltitude

        # set a nan value so that it can be returned by functions
        self.__nan = 1e300*1e300 - 1e300*1e300

    def __cNan(self, num):
        """__cNan is a private method that checks a number to
        see if it is equal to the IEEE-754 standard NaN value. If it
        isn't, it returns the orignal number. If it is, it returns the
        string 'missing' so that Madrigal can recognize a missing value.

        Inputs:
            num - the number to be checked

        Returns: num if num is not NaN, and 'missing' if it is

         Algorithm: if both (num < 10.0) and (num > -10.0)
         are False, num is nan. Algorithm written by Bill Rideout.
        """
        if ( (num < 10.0) == False and (num > -10.0) == False ):
            return 'missing'
        else:
            return num
        
    def __createExpCatalog(self, startTime, endTime):
        """ __createExpCatalog is a private method that creates a
            properly formatted experiment catalog so that it can then
            be put into a Madrigal experiment.

            Inputs:
                 startTime - start time of the experiment
                 endTime - end time of the experiment

            Returns: a Madrigal experiment catalog record object

            Affects: nothing

            Exceptions: none
            
            (note: this should contain much more information;
             it should be expanded to contain more pass data)
        """
        
        # make the string that will form the meat of the record
        s = 'KRECC       2001 Catalogue Record, Version 1\n'
        s += 'KINSTE      %s'% self.__kinst +' MIDAS_M Beacon Receiver\n'
        s += 'MODEXP         0 Usual mode of operation\n'
        s += 'IBYRE       %s'% startTime.year +' Beginning year\n'
        s += 'IBDTE       %s'% startTime.strftime("%m%d") +' Beginning month and day\n'
        s += 'IBHME       %s'% startTime.strftime("%H%M") +' Beginning UT hour and minute\n'
        s += 'IBCSE       %s%s'% (startTime.second,startTime.microsecond/10000) +' Beginning centisecond\n'
        s += 'IBYRE       %s'% endTime.year +' Ending year\n'
        s += 'IBDTE       %s'% endTime.strftime("%m%d") +' Ending month and day\n'
        s += 'IBHME       %s'% endTime.strftime("%H%M") +' Ending UT hour and minute\n'
        s += 'IBCSE       %s%s'% (endTime.second,endTime.microsecond/10000) +' Ending centisecond\n'
        s += 'C\n'
        s += 'CPURP    Beacon receiver satellite pass\n'
        s += 'C\n'
        s += 'CIREM    See http://www.openradar.org\n'
        s += 'CSREM    See http://www.openradar.org\n'
        s += 'C\n'
        s += 'CPI      <FIXME>\n'
        # strip newlines, make every line 80 columns
        s = self.__padto80(s)

        return madrigal.cedar.MadrigalCatalogRecord(self.__kinst,
                                                    1, # modexp of 0
                                                    startTime.year, startTime.month,
                                                    startTime.day, startTime.hour,
                                                    startTime.minute, startTime.second,
                                                    startTime.microsecond / 10000,
                                                    endTime.year, endTime.month,
                                                    endTime.day, endTime.hour,
                                                    endTime.minute, endTime.second,
                                                    endTime.microsecond / 10000, s)
    

    def __createExpHeader(self, startTime, endTime):
        """ __createExpHeader is a private method that creates a
            properly formatted experiment header record so that it can
            then be put into a Madrigal experiment.

            Inputs:
                 startTime - start time of the experiment
                 endTime - end time of the experiment

            Returns: a Madrigal experiment header record object

            Affects: nothing

            Exceptions: none
            
            (note: this should contain much more information;
             it should be expanded to contain more pass data)
        """
        s = 'KRECH       3002  Header record, Version 1\n'
        s += 'KINST       %s'% self.__kinst +' MIDAS_M Beacon Receiver\n'
        s += 'KINDAT      %s'% self.__kindat +' '+self.__kindatDesc + '\n'
        # strip newlines, make every line 80 columns
        s = self.__padto80(s)

        return madrigal.cedar.MadrigalHeaderRecord(self.__kinst, self.__kindat,
                                                   startTime.year, startTime.month,
                                                   startTime.day, startTime.hour,
                                                   startTime.minute, startTime.second,
                                                   startTime.microsecond / 10000,
                                                   endTime.year, endTime.month,
                                                   endTime.day, endTime.hour,
                                                   endTime.minute, endTime.second,
                                                   endTime.microsecond / 10000,
                                                   18, # num. of 1D parameters
                                                   5, # num. of 2D parameters
                                                   s)

    def __padto80(self, text):
        """__padto80 is a private method that takes input text, and
        formats it for a catalog or header record.  Each line must
        contain with 80 characters or less.  All line feeds will be
        removed, and each line will be padded to 80 characters with
        spaces if needed. Taken from example code in the Madrigal
        documentation.

        Inputs:
            text - the text to be padded, where lines are split by newlines

        Returns: the text as a single line, with each virtual line padded
                 to 80 characters

        Affects: nothing

        Exceptions: none
        """
        lines = text.split('\n')
        retStr = ''
        for line in lines:
            if len(line) > 80:
                raise ValueError, 'line <%s> greater than 80 charaters' % (line)
            retStr += line + ' ' * (80 - len(line))
        return retStr

    def __createPlotPage(self, expPath, filePath, hasTEC):
        """__createPlotPage is a private method that generates
        an HTML page in the Madrigal experiment page with links
        for plots. This allows the user to view the plots that are in
        the experiment directory.

        Inputs:
            expPath - a string containing the path to the current
                      Madrigal experiment
            filePath - a string containing the filename that will be
                       written. 
            hasTEC - a boolean type that is true if TEC values were
                     calculated and TEC plots were made

        Returns: nothing

        Affects: adds a file named plots/index.html with links to plots.
                 A link to this file will show up in the web interface.

        Exceptions: none
        """

        # all plot names and descriptions
        plotNames = { 's4':'S4 plot', 'sigmaphi':'Sigma-phi plot',
                      'azel':'Azimuth and elevation plot' }
        # only add TEC plots if TEC exists
        if (hasTEC == True):
            plotNames['slanttec'] = 'Slant-path TEC'
            plotNames['verticaltec'] = 'Vertical TEC'
            plotNames['vtecmap'] = 'Vertical TEC map plot' 

        # create HTML file
        s = "<html>\n<head>\n"
        s += "<title>Plots for experiment %s</title>\n"%expPath
        s += "</head>\n<body BGCOLOR=#FFFFFF LINK=#008000 VLINK=#003366>\n"
        s += "<center>\n<h2>Plots for experiment %s</h2>\n"%expPath
        s += "<p><table border='0' cellpadding='2'>\n"

        # add links and thumbs of all plots
        for name in plotNames.keys():
            s += "<tr><td><a href='%s'>" % (name+'.png')
            s += "<img src='%s'></td>" % (name+'-small.png')
            s += "<td><a href='%s'>%s</a></td></tr>\n" % (name+'.png', plotNames[name])

        s += "</table>\n</center>\n</body>\n</html>\n"

        # write file
        htmlFile = open(filePath,'w')
        htmlFile.write(s)
        htmlFile.close()
        
         
    def load(self, h5file, plots=True, plotDpi=100, thumbDpi=20):
        """ load parses the HDF5 file and loads it into Madrigal.

        Inputs:
            hdf5File - a PyTables HDF5 file object that will be
                       read in
            plots - a boolean type that determines whether to
                    draw and insert plots into the experiment path
            plotDpi - an integer that defines the resolution of the
                      generated plots that will be put in the experiment
                      directory
            thumbDpi - an integer that defines the resolution of the
                       generated thumbnails that will also be put
                       in the experiment directory 

        Returns: the path of the created experiment file

        Affects: the local Madrigal database; no public or private
                 attributes will be affected

        Exceptions: none

        Note: $MADROOT/bin/updateMaster must be run after loading the
        experiment in order to have it show up in the database
        """

        # read in all relevant arrays

        print 'Reading HDF5 file'
        # metadata
        fsatArray = numpy.array(h5file.root.Metadata.Fsat.read())
        # fields are commented out if not used yet
        #bandwidthArray = numpy.array(h5file.root.Metadata.Bandwidth.read())
        startUTCStr = h5file.root.Metadata.StartUTC.read()[0]
        sampInterval = int(h5file.root.Metadata.SampInterval.read()[0])
        outInterval = int(h5file.root.Metadata.OutInterval.read()[0])
        s4Interval = int(h5file.root.Metadata.S4Interval.read()[0])
        sigmaPhiInterval = int(h5file.root.Metadata.SigmaPhiInterval.read()[0])
        tecInterval = int(h5file.root.Metadata.TECInterval.read()[0])
        tleList = h5file.root.Metadata.TLE.read()

        # mid-level data
        frequencyArray = numpy.array(h5file.root.Data.Frequency.read())
        # get sample number and number of frequencies
        numFreq = int(h5file.root.Data.Voltage.shape[0])
        samples = int(h5file.root.Data.Voltage.shape[1])
        
        # high-level data
        s4Array = numpy.array(h5file.root.HighLevelData.S4.read())
        sigmaPhiArray = numpy.array(h5file.root.HighLevelData.SigmaPhi.read())
        slantTECArray = numpy.array(h5file.root.HighLevelData.SlantTEC.read())
        verticalTECArray = numpy.array(h5file.root.HighLevelData.VerticalTEC.read())

        # ephemeris
        azimuthArray = numpy.array(h5file.root.Ephemeris.Azimuth.read())
        elevationArray = numpy.array(h5file.root.Ephemeris.Elevation.read())
        altitudeArray = numpy.array(h5file.root.Ephemeris.Altitude.read())
        rangeArray = numpy.array(h5file.root.Ephemeris.Range.read())
        satLatitudeArray = numpy.array(h5file.root.Ephemeris.Latitude.read())
        satLongitudeArray = numpy.array(h5file.root.Ephemeris.Longitude.read())

        # receiver configuration
        #antDescStr = h5file.root.ReceiverConfig.AntennaDescription.read()[0]
        receiverLatitude = h5file.root.ReceiverConfig.Latitude.read()[0]
        receiverLongitude = h5file.root.ReceiverConfig.Longitude.read()[0]
        #receiverAltitude = h5file.root.ReceiverConfig.Altitude.read()[0]
        #receiverAzimuth = h5file.root.ReceiverConfig.Azimuth.read()[0]
        #receiverElevation = h5file.root.ReceiverConfig.Elevation.read()[0]

        # close file
        h5file.close()

        print 'Calculating pierce points'
        # find E & F region pierce points for all ephemeris
        ePLatitudeArray = numpy.zeros(samples)
        ePLongitudeArray = numpy.zeros(samples)
        fPLatitudeArray = numpy.zeros(samples)
        fPLongitudeArray = numpy.zeros(samples)
        for i in range(samples):
            (ePLat,ePLon) = find_pierce_point(receiverLatitude,
                                              receiverLongitude,
                                              satLatitudeArray[i],
                                              satLongitudeArray[i],
                                              altitudeArray[i],
                                              self.__eAltitude)
            (fPLat,fPLon) = find_pierce_point(receiverLatitude,
                                              receiverLongitude,
                                              satLatitudeArray[i],
                                              satLongitudeArray[i],
                                              altitudeArray[i],
                                              self.__fAltitude)
            ePLatitudeArray[i] = ePLat
            ePLongitudeArray[i] = ePLon
            fPLatitudeArray[i] = fPLat
            fPLongitudeArray[i] = fPLon
            
        # decide which TEC parameter to use from the four provided
        # for now, a simple algorithm: whichever has the most data points
        # if none have data, hasTEC is False and no TEC is
        # inserted or plotted
        numTECValues = [0,0,0,0]
        for i in range(len(slantTECArray)):
            for j in range(4):
                if (self.__cNan(slantTECArray[i][j][0]) != 'missing'):
                    numTECValues[j] += 1
        if (max(numTECValues) == 0):
            hasTEC = False
        else:
            hasTEC = True
            tecParameter = numTECValues.index(max(numTECValues))
        
        # make  time object for start time (which is formatted
        # according to ISO 8601) and strings to represent it
        startTime = datetime.datetime(int(startUTCStr[0:4]),          # year
                                      int(startUTCStr[4:6]),          # month
                                      int(startUTCStr[6:8]),          # day
                                      int(startUTCStr[9:11]),         # hour
                                      int(startUTCStr[11:13]),        # min
                                      int(startUTCStr[13:15]),        # sec
                                      int(startUTCStr[16:-1]) / 1000) # us
        dateStr = startTime.strftime('%y%m%d')
        dateStrMadrigal = startTime.strftime('%d%b%y').lower()
        
        # figure out end time based on the number of samples and the 
        # sampling rate
        duration = datetime.timedelta(
            microseconds = samples * sampInterval / 1000)
        endTime = startTime + duration

        # find instrument mnemonic
        madInstrument = madrigal.metadata.MadrigalInstrument(self.__madrigalDB)
        instMnemonic = madInstrument.getInstrumentMnemonic(self.__kinst)

        # put down the experiment filename
        expBasename = '%s%sa.000' % (instMnemonic, dateStr)
        expFilename = '/tmp/'+expBasename

        # create a Madrigal file; will be later loaded into Madrigal and
        # then deleted
        madFile = madrigal.cedar.MadrigalCedarFile(expFilename, True)

        print expFilename + ' created.'

        # add catalog and header files
        madFile.append(self.__createExpCatalog(startTime, endTime))
        madFile.append(self.__createExpHeader(startTime, endTime))

        # create a list of 2-tuples, each containing the start time
        # and data type of a parameter. The type is designated as
        # 1: S4, 2: sigma-phi and 3: TEC. Once sorted, this list will
        # determine in what order the parameters are inserted as records
        # into the Madrigal file.
        paramList = []
        for i in range(s4Array.shape[1]):
            td = datetime.timedelta(microseconds = i * s4Interval / 1000)
            paramList.append((startTime+td,1))
        for i in range(sigmaPhiArray.shape[1]):
            td = datetime.timedelta(microseconds = i * sigmaPhiInterval / 1000)
            paramList.append((startTime+td,2))
        # only do this if TEC exists
        if (hasTEC == True):
            for i in range(verticalTECArray.shape[0]):
                td = datetime.timedelta(microseconds = i * tecInterval / 1000)
                paramList.append((startTime+td,3))
        paramList.sort()

        # initialize indices for data; each is incremented once a record
        # of it is pulled
        s4Index, sigmaPhiIndex, tecIndex = 0,0,0

        # iterate through this list and add records to Madrigal
        print 'Adding records:'
        for i in range(len(paramList)):
            if (i % 250 == 0):
                print str(i+1)+' of '+str(len(paramList))
            # shared 1D and 2D parameters
            oneDList = ['gdalt','range','azm','elm','gdlat',
                        'glon','eplat','eplon','fplat','fplon']
            twoDList = ['tfreq']
            # figure out the start and end times for the record, its interval
            # and index, and what parameters will be included
            if (paramList[i][1] == 1):
                interval = s4Interval
                index = s4Index
                twoDList.extend(['s4','ds4'])
            elif (paramList[i][1] == 2):
                interval = sigmaPhiInterval
                index = sigmaPhiIndex
                twoDList.extend(['sgmph','dsgmph'])
            elif (paramList[i][1] == 3):
                interval = tecInterval
                index = tecIndex
                oneDList.extend(['rtec','drtec','rneli','drneli'])
                
            recBeginTime = startTime + datetime.timedelta(microseconds = index * interval / 1000)
            recEndTime = recBeginTime + datetime.timedelta(microseconds = interval / 1000)
            # figure out which indices to use for ephemeris and frequency data
            # (if more than one index for each sample duration, use
            # the middle index)
            ephMultiplier = interval / sampInterval
            ephIndex = index * ephMultiplier + int(ephMultiplier / 2)
            # interval might not necessarily be an integer multiple
            # of the frequency/SNR interval
            freqMultiplier = float(interval) / float(outInterval)
            freqIndex = int( (float(index) + 0.5) * freqMultiplier )
            # make sure indices aren't too big (this can happen only
            # at the end)
            if (ephIndex >= samples):
                ephIndex = samples - 1
            if (freqIndex >= frequencyArray.shape[1]):
                freqIndex = frequencyArray.shape[1] - 1
            
            # create the data record
            dataRec = madrigal.cedar.MadrigalDataRecord(self.__kinst,
                                                        self.__kindat,
                                                        recBeginTime.year,
                                                        recBeginTime.month,
                                                        recBeginTime.day,
                                                        recBeginTime.hour,
                                                        recBeginTime.minute,
                                                        recBeginTime.second,
                                                        recBeginTime.microsecond / 10000,
                                                        recEndTime.year,
                                                        recEndTime.month,
                                                        recEndTime.day,
                                                        recEndTime.hour,
                                                        recEndTime.minute,
                                                        recEndTime.second,
                                                        recEndTime.microsecond / 10000,
                                                        oneDList,
                                                        twoDList,
                                                        numFreq)

            # add all common parameters to record
            # (if nan, __cNan() will let Madrigal recognize 'missing')
            dataRec.set1D('gdalt',self.__cNan(altitudeArray[ephIndex] / 1000))
            dataRec.set1D('range',self.__cNan(rangeArray[ephIndex] / 1000))
            dataRec.set1D('azm',self.__cNan(azimuthArray[ephIndex]))
            dataRec.set1D('elm',self.__cNan(elevationArray[ephIndex]))
            dataRec.set1D('gdlat',self.__cNan(satLatitudeArray[ephIndex]))
            dataRec.set1D('glon',self.__cNan(satLongitudeArray[ephIndex]))
            dataRec.set1D('eplat',self.__cNan(ePLatitudeArray[ephIndex]))
            dataRec.set1D('eplon',self.__cNan(ePLongitudeArray[ephIndex]))
            dataRec.set1D('fplat',self.__cNan(fPLatitudeArray[ephIndex]))
            dataRec.set1D('fplon',self.__cNan(fPLongitudeArray[ephIndex]))
                
            # set 2D parameters
            for f in range(numFreq):
                # set frequency 
                dataRec.set2D('tfreq',f,self.__cNan(fsatArray[f]))
                # give S4 or sigma-phi if appropriate
                if (paramList[i][1] == 1):
                    dataRec.set2D('s4',f,self.__cNan(s4Array[f,index,0]))
                    dataRec.set2D('ds4',f,self.__cNan(1.0))
                if (paramList[i][1] == 2):
                    dataRec.set2D('sgmph',f,self.__cNan(sigmaPhiArray[f,index,0]))
                    dataRec.set2D('dsgmph',f,self.__cNan(1.0))
            
            # put in TEC parameters if appropriate
            if (paramList[i][1] == 3):
                dataRec.set1D('rtec',self.__cNan(verticalTECArray[index,tecParameter,0]))
                dataRec.set1D('drtec',self.__cNan(verticalTECArray[index,tecParameter,1]))
                dataRec.set1D('rneli',self.__cNan(slantTECArray[index,tecParameter,0]))
                dataRec.set1D('drneli',self.__cNan(slantTECArray[index,tecParameter,1]))

            # append the record
            madFile.append(dataRec)
            # dump to the file every 1000 records so that the memory
            # isn't overloaded
            if (i % 1000 == 0 and i != 0):
                madFile.dump()
                
            # increment the proper index since the record was pulled
            if (paramList[i][1] == 1):
                s4Index += 1
            elif (paramList[i][1] == 2):
                sigmaPhiIndex += 1
            elif (paramList[i][1] == 3):
                tecIndex += 1

        # dump the last bit
        madFile.dump()
              
        # create a new Madrigal experiment
        dbAdmin = madrigal.admin.MadrigalDBAdmin()

        # see how many experiments already exist for the day so that
        # an optional character can be appended (this code makes
        # the assumption that the extra character is always appended serially)
        index = 0
        count = 0 
        madrigalExp = madrigal.metadata.MadrigalExperiment(self.__madrigalDB)
        while 1:
            url = madrigalExp.getExpUrlByPosition(index)
            if url == None:
                break
            if (url.find(instMnemonic+'/'+dateStrMadrigal) != -1):
                count += 1
            index += 1

        # too many experiments
        if (count > 26):
            errStr = 'Limit of 27 experiments reached for '+dateStrMadrigal
            self.__log.error(errStr)
            raise IOError, errStr
        
        if (count == 0):
            extraChar = ''
        else:
            extraChar = chr(96+count)

        # make sure all metadata strings are without commas
        self.__expTitle = self.__expTitle.replace(',', ';')

        # create the experiment
        expPath = dbAdmin.createMadrigalExperiment(expFilename,
                                                   self.__expTitle,
                                                   0, # public
                                                   self.__kindatType,
                                                   self.__kinst,
                                                   1, # default
                                                   extraChar)
    
        print 'Experiment loaded into ' + expPath
        os.remove(expFilename)
        print expFilename + ' deleted.'

        # no plots: stop
        if (plots == False):
            return expPath
        
        # make plots and put into Madrigal
        bp = beacon_plot.BeaconPlot(dpi=plotDpi, thumbDpi=thumbDpi)
        plotPath = expPath+'/plots'
        mkdirStr = 'mkdir '+expPath+'/plots'
        if (os.system(mkdirStr) != 0):
            raise IOError, 'Could not create /plots subdirectory in ' + expPath

        freqNames = ['%s MHz'%(int(round(i/1e6))) for i in fsatArray]

        print 'Generating S4 plot'
        # take only the values and not the errors
        s4Values = s4Array[:,:,0]
        bp.timePlot(s4Values, s4Interval, freqNames, 'S4', startTime,
                    tleList[0], plotPath+'/s4.png', yLimits=(0,1))
        
        print 'Generating Sigma-phi plot'
        sphValues = sigmaPhiArray[:,:,0]
        bp.timePlot(sphValues, sigmaPhiInterval, freqNames, 'Sigma-phi',
                    startTime, tleList[0], plotPath+'/sigmaphi.png',
                    units='rad', yLimits=(0,math.pi))

        if (hasTEC == True):
            print 'Generating TEC plots'
            slantTECValues = slantTECArray[:,tecParameter,0]
            verticalTECValues = verticalTECArray[:,tecParameter,0]
            bp.timePlot([slantTECValues], tecInterval, [''],
                        'Relative slant-path TEC',
                        startTime, tleList[0], plotPath+'/slanttec.png',
                        units='e- / m^2', legend=False)
            bp.timePlot([verticalTECValues], tecInterval, [''],
                        'Relative vertical TEC',
                        startTime, tleList[0], plotPath+'/verticaltec.png',
                        units='TECU',legend=False)

        print 'Generating az/el plot'
        bp.azElPlot(azimuthArray, elevationArray, sampInterval,
                    startTime, tleList[0], plotPath+'/azel.png')

        if (hasTEC == True):
            print 'Generating vertical TEC map plot'
            # see if vertical TEC values need to be scaled because
            # there are less of them than ephemeris values
            if (tecInterval == sampInterval):
                vScaledTECValues = verticalTECValues
            else:
                multiple = tecInterval / sampInterval
                vScaledTECValues = numpy.zeros(samples)
                for i in range(samples / multiple):
                    for j in range(len(multiple)):
                        vScaledTECValues[i*multiple+j] = verticalTECValues[i]

            bp.mapPlot(ePLatitudeArray, ePLongitudeArray,
                       fPLatitudeArray, fPLongitudeArray,
                       receiverLatitude, receiverLongitude,
                       sampInterval, startTime,
                       tleList[0], vScaledTECValues,
                       'relative vertical TEC',
                       plotPath+'/vtecmap.png', units='TECU')

        # create page to link to plots
        self.__createPlotPage(expPath, '/tmp/index.html', hasTEC)
        dbAdmin.addWebFile(expPath, '/tmp/index.html', 'plots/')
        
        # required to have Madrigal recognize the new experiment
        print 'Updating master metadata'
        dbAdmin.updateMaster()

        print 'Done.'
        return expPath
