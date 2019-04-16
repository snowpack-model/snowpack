"""IOUtilsWrap.pyx: This file wraps the c-function getLibVersion to the python-function PyGetLibVersion.
   Author: Thiemo Theile
   Date created: 4/1/2019
"""

# distutils: language = c++

cimport IOUtils

def PyGetLibVersion(short_version=0):
    return IOUtils.getLibVersion(short_version)