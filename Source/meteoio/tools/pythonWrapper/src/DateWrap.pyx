"""DateWrap.pyx: This file wraps the c++-class Date (frommeteoio/dataClasses/Date.h) to the python-class PyDate.
   Author: Thiemo Theile
   Date created: 4/1/2019
"""

# distutils: language = c++

from Date cimport Date

#this is dublicated code from Date.pxd to make the enums known here (I found no elegant way to do this...)
ctypedef enum DATE_FORMATS:
    ISO, # ISO 8601 extended format combined date: YYYY-MM-DDTHH:mm:SS.sss (fields might be dropped, in the least to the most significant order)
    ISO_TZ, # ISO 8601 format (same as ISO) but with time zone specification
    ISO_Z, # ISO 8601 format, forcing GMT and Zulu (Z) timezone specification
    FULL, # ISO 8601 followed by the julian date (in parenthesis)
    NUM, # ISO 8601 basic format date: YYYYMMDDHHmmSS (fields might be dropped, in the least to the most significant order)
    DIN, #DIN5008 format: DD.MM.YYYY HH:MM:SS.sss
    ISO_WEEK, # ISO 8601 week date: YYYY-Www-D (for example: 2014-W41-1)
    ISO_DATE # ISO 8601 date format without the time (ie YYYY-MM-DD)

cdef class PyDate:
    cdef Date c_date

    def stringToType(self, typeString):
        type = ISO
        if typeString=="ISO":
            type=ISO
        if typeString=="ISO_TZ":
            type=ISO_TZ
        if typeString=="ISO_Z":
            type=ISO_Z
        if typeString=="FULL":
            type=FULL
        if typeString=="NUM":
            type=NUM
        if typeString=="DIN":
            type=DIN
        if typeString=="ISO_WEEK":
            type=ISO_WEEK
        if typeString=="ISO_DATE":
            type=ISO_DATE
        return type

    def __cinit__(self):
        self.c_date = Date()
    #overloading of constructor not possible with python. one way would be to use classmethods
    
    def setDate(self, year, month, day, hour, minute, in_timezone, in_dst):
        #print("setting date...")
        self.c_date.setDate(year, month, day, hour, minute, in_timezone, in_dst)
        #print(self.c_date.toString())

    def toString(self, typeString="ISO", gmt=0):
        type = self.stringToType(typeString)
        return self.c_date.toString(type, gmt)
       
    #def __add__(PyDate self, PyDate anotherDate):
    #    newDate = PyDate()
    #    newDate.c_date = self.c_date + anotherDate.c_date
    #    return newDate
        
    def __add__(PyDate self, double dtDays):
        newDate = PyDate()
        newDate.c_date = self.c_date + dtDays
        return newDate
        
    def __sub__(PyDate self, double dtDays):
        newDate = PyDate()
        newDate.c_date = self.c_date - dtDays
        return newDate    
         
    def __lt__(self, PyDate anotherDate):
        return (self.c_date < anotherDate.c_date)
        
    def __le__(self, PyDate anotherDate):
        return (self.c_date <= anotherDate.c_date)
        
    def __gt__(self, PyDate anotherDate):
        return (self.c_date > anotherDate.c_date)
        
    def __ge__(self, PyDate anotherDate):
        return (self.c_date >= anotherDate.c_date)
        
    def __eq__(self, PyDate anotherDate):
        return (self.c_date == anotherDate.c_date)

         
######### useful helper functions: #####################

def createDateList(dStart,dStop,dtDays):
    """ Create a list of PyDate-objects starting from dStart until dStop with the interval dtDays """
    dateList=[]
    dI = dStart
    while( dI < dStop ):
        dateList.append(dI)
        dI = dI + dtDays
    return dateList