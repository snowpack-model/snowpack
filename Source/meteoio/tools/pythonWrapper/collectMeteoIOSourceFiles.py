import os
from glob import glob
import os.path as path

setOfFilesToSkip={
    "../../meteoio/IOHandler.cmake.cc",
    "../../meteoio/plugins/CosmoXMLIO.cc",
    "../../meteoio/plugins/PSQLIO.cc",
    "../../meteoio/plugins/NetCDFIO.cc",
    "../../meteoio/plugins/libncpp.cc",
    "../../meteoio/plugins/ImisIO.cc",
    "../../meteoio/plugins/DBO.cc",
    "../../meteoio/plugins/GRIBIO.cc",
    "../../meteoio/plugins/GSNIO.cc",
    "../../meteoio/plugins/OshdIO.cc",
    "../../meteoio/plugins/template.cc",
    "../../meteoio/plugins/BormaIO.cc",
    "../../meteoio/plugins/SASEIO.cc",
    "../../meteoio/plugins/PNGIO.cc",
    "../../meteoio/plugins/libMatioWrapper.cc",
    "../../meteoio/meteoStats/RandomNumberGenerator.cc",
    "../../meteoio/jnative/DEMLoader.cc"
}

def getListOfSourcefiles(filesToSkip=setOfFilesToSkip):
    files = []
    current_dir = os.getcwd()
    start_dir = path.abspath(path.join(current_dir,"../../meteoio"))
    pattern = "*.cc"

    for dir, _, _ in os.walk(start_dir):
        files.extend(glob(os.path.join(dir, pattern)))

    listOfFiles=[]
    for file in files:
        splitted = file.split("\\")
        newFile = "../../" + splitted[3]
        for i in range(4, len(splitted)):
            split = splitted[i]
            newFile = newFile + "/" + split
        newFile = newFile + ""
        if(filesToSkip.isdisjoint({newFile})):
            listOfFiles.append(newFile)
        else:
            print("we don't want this file:"+newFile)
    return listOfFiles
