# distutils: language = c++

from Date cimport Date

cdef class PyDate:
    cdef Date c_date

    def stringToType(self, typeString):
        type=0
        if typeString=="ISO":
            type=0
        if typeString=="ISO_TZ":
            type=1
        if typeString=="ISO_Z":
            type=2
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



         
        
def createDateList(dStart,dStop,dtDays):
    dateList=[]
    dI = dStart
    while( dI < dStop ):
        #print(dI.toString())
        dateList.append(dI)
        dI = dI + dtDays
    return dateList