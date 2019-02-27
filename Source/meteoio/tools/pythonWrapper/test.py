import meteoioWrapper as mio

print(dir(mio))

test = mio.PyCoords()
test.setLatLon(10,11,1560,0)
print(test.toString("FULL"))

cfg = mio.PyConfig("io.ini".encode('utf-8'))
TZ = cfg.get("TIME_ZONE".encode('utf-8'), "Input".encode('utf-8'),"default".encode('utf-8'))

io = mio.PyIOManager(cfg)

Tstep = 1/24. #sampling rate in days

d_start = mio.PyDate()
d_start.setDate(2009,1,1,12,0,float(TZ),0)
d_end = mio.PyDate()
d_end.setDate(2009,1,1,15,0,float(TZ),0)

print( "Powered by MeteoIO " + str(mio.PyGetLibVersion()) )
print( "Reading data from " + str(d_start.toString("ISO".encode('utf-8'))) + " to " + str(d_end.toString("ISO".encode('utf-8'))) )

dateList = mio.createDateList(d_start,d_end,Tstep)

doResample = False

vecVecMeteo=[] #so we can keep and output the data that has been read

if doResample:
    mapIDs = dict()
    insert_position=0
    for date in dateList:
        meteoSet = io.getMeteoData(date) #read 1 timestep at once, forcing resampling to the timestep
        for meteoData in meteoSet:
            print(meteoData.toString())
            stationID = meteoData.getStationID()
            if not (stationID in mapIDs):
                mapIDs[stationID]=insert_position
                insert_position=insert_position+1
                vecVecMeteo.append( [] )
            vecVecMeteo[ mapIDs[stationID] ].append(meteoData)
else:
    vecVecMeteo = io.getMeteoDataRange(d_start, d_end) #This would be the call that does NOT resample the data, instead of the above "for" loop

print("Writing output data")
io.writeMeteoData(vecVecMeteo,"".encode('utf-8'))