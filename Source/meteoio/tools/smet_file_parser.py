import numpy as np
from datetime import datetime


class SMETFileParser:
    """
    This class is intended to be used to read SMET files.
    """
    
    def __init__(self, fileName = None):
        """
        Constructor. The full name of the file to be parsed should be passed as 
        argument. 
        """
        if fileName is not None:
            self.setFile(fileName)
        
        
    def setFile(self, fileName):
        """
        Parses the SMET file whose name is given as input argument.
        """
        self._fileName = fileName
        self._parseFile()
        
    
    def getHeader(self):
        """
        Returns the key-value pairs defined in the file header as a dictionnary.
        """
        if self._fileName is not None:
            return self._header
        else:
            raise RuntimeError("Cannot return file header: no file parsed. " +
                               "Please specify a file to parse using method setFile()")
    
    
    def getHeaderKeys(self):
        """
        Returns the list of keys present in the file header.
        """
        if self._fileName is not None:
            return self._header.keys()
        else:
            raise RuntimeError("Cannot return file header keys: no file parsed. " +
                               "Please specify a file to parse using method setFile()")
    
    
    def getHeaderValue(self, key):
        """
        Returns the value associated with the header key passed as input 
        argument.
        """
        key = key.lower()
        if self._fileName is not None and self._header.has_key(key):
            return self._header[key]
        elif self._fileName is not None:
            raise RuntimeError("Cannot return value associated with key " +
                               key + " in header of file " + self._fileName +
                               ": key does not exist")
        else:
            raise RuntimeError("Cannot return value associated with header key " +
                               key + ": no file parsed. Please specify a file " +
                               "to parse using method setFile()")
    
    
    def getFieldNames(self):
        """
        Returns the names of the fields contained in the SMET file.
        """
        if self._fileName is not None:
            return self._fields
        else:
            raise RuntimeError("Cannot return field names: no file parsed. " +
                               "Please specify a file to parse using method setFile()")
    
    
    def getFieldData(self, fieldName):
        """
        Returns the values of the field whose name is passed as input argument.
        """
        fieldName = fieldName.upper()
        if self._fileName is not None and fieldName in self._data.dtype.names:
            return self._data[fieldName]
        elif self._fileName is not None:
            raise RuntimeError("Cannot return values of field " + fieldName +
                               " in file " + self._fileName + ": field does not exist")
        else:
            raise RuntimeError("Cannot return values of field " + fieldName +
                               ": no file parsed. Please specify a file " +
                               "to parse using method setFile()")
            
    
    def getData(self):
        """
        Returns all the values contained in section DATA of the file as a 
        numpy.ndarray object.
        """
        if self._fileName is not None:
            return self._data
        else:
            raise RuntimeError("Cannot return field data: no file parsed. " +
                               "Please specify a file to parse using method setFile()")
    
    
    def _parseFile(self):
        iDataStart = self._parseHeaderSection()
        self._parseDataSection(iDataStart + 1)
        
        
    def _parseHeaderSection(self):
        keys   = []
        values = []
        with open(self._fileName, 'r') as f:
            iDataLine = self._reachSectionStart("HEADER", f)
            for iLine, line in enumerate(f):
                line = line.split('#', 1)[0].strip()
                if line.upper() == "[DATA]":
                    iDataLine += iLine + 1
                    break
                elif line == "":
                    continue
                elif line.count('=') != 1:
                    raise RuntimeError("Cannot parse line " + str(iLine + iDataLine + 2) +
                                       " in file " + self._fileName + ": key must be " +
                                       "separated from its corresponding value using sign '='")
            
                key, value = line.split('=')
                keys.append(key.strip().lower())
                try:
                    values.append(float(value))
                except:
                    values.append(value.strip())
            
        self._header = dict(zip(keys, values)) 
        try:
            self._fields = self._header["fields"].upper().split()
            self._nodata = self._header["nodata"]
        except KeyError as err:
            raise RuntimeError("Missing key '" + str(err) + "' in HEADER section "
                               "of file " + self._fileName)
        
        return iDataLine
        
        
    def _parseDataSection(self, numHeaderLines):
        dtypes = ['f' if f != "TIMESTAMP" else 'O' for f in self._fields]
        dtypes = ','.join(dtypes)
        kwargs = dict(comments="#", skip_header=numHeaderLines, dtype=dtypes)
        if "TIMESTAMP" in self._fields:
            convertTime = lambda x: datetime.strptime(x, '%Y-%m-%dT%H:%M:%S')
            iTimeField  = self._fields.index("TIMESTAMP")
            kwargs["converters"] = {iTimeField: convertTime}
        
        try:
            self._data = np.genfromtxt(self._fileName, **kwargs)
        except Exception as e:
            raise RuntimeError("Error while parsing DATA section in file " + 
                               self._fileName + ": " + str(e))
            
        self._data.dtype.names = self._fields
        for field in set(self._fields) - set(["TIMESTAMP"]):
            self._data[field][self._data[field] == self._nodata] = np.NaN
        
        
    def _reachSectionStart(self, sectionName, fileObject):
        fileObject.seek(0)
        sectionName = "[" + sectionName.upper() + "]"
        for iLine, line in enumerate(fileObject):
            if line.split('#', 1)[0].strip().upper() == sectionName:
                return iLine
        
        raise EOFError("Cannot read section " + sectionName.upper() + " in file "
                       + self._fileName + ": section does not exist")
