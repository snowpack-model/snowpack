"""
test.py: This script is an example of how to use the python-wrapper for meteoIO. This script does the same as
        the data_converter (which is implemented in c++ and can be found in \doc\examples). It reads in, filters
        and writes out a data file according to the ini-file ('io.ini').

Author: Thiemo Theile
Date created: 4/1/2019
Python Version: 2.7 or 3.x
"""


import meteoioWrapper as mio

#print what is inside the meteoioWrapper (classes and functions)
print(dir(mio))

cfg = mio.PyConfig(b"io.ini")
TZ = cfg.get(b"TIME_ZONE", b"Input", b"default")

io = mio.PyIOManager(cfg)

Tstep = 1/24. #sampling rate in days

#define the time range for which the meteo data is processed (d_start until d_end)
d_start = mio.PyDate()
d_start.setDate(2009,1,1,12,0,float(TZ),0)
d_end = mio.PyDate()
d_end.setDate(2009,1,1,15,0,float(TZ),0)

print( "Powered by MeteoIO " + str(mio.PyGetLibVersion()) )
print( "Reading data from " + str(d_start.toString(b"ISO")) + " to " + str(d_end.toString(b"ISO_TZ")) )

dateList = mio.createDateList(d_start,d_end,Tstep)

doResample = True #set to True if you want the data to be resampled to the sampling rate Tstep

vecVecMeteo=[] #so we can keep and output the data that has been read

if doResample:
    mapIDs = dict()
    insert_position=0
    for date in dateList:
        meteoSet = io.getMeteoData(date) #read 1 timestep at once, forcing resampling to the timestep
        for meteoData in meteoSet:
            stationID = meteoData.getStationID()
            if not (stationID in mapIDs):
                mapIDs[stationID]=insert_position
                insert_position=insert_position+1
                vecVecMeteo.append( [] )
            vecVecMeteo[ mapIDs[stationID] ].append(meteoData)
else:
    vecVecMeteo = io.getMeteoDataRange(d_start, d_end) #This would be the call that does NOT resample the data, instead of the above "for" loop

#example of how to access meta data:
stationID = vecVecMeteo[0][0].meta.getStationID()
dateOfFirstDataPoint = vecVecMeteo[0][0].date.toString()
print("station name: "+str(stationID))
#example of how to get certain meteo-data (wind speed, temperature,...) for certain station
temperatures = mio.getTimeSeriesOfCertainParameter(b"TA",stationID, vecVecMeteo)
print(temperatures)

print("Writing output data")
io.writeMeteoData(vecVecMeteo,b"")