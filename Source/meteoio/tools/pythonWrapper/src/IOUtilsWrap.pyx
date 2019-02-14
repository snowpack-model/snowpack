# distutils: language = c++

cimport IOUtils

def PyGetLibVersion(short_version=0):
    return IOUtils.getLibVersion(short_version)