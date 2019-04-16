Date: 13.02.2019
Author: Thiemo Theile


This python-wrapper uses Cython (https://cython.org/) to make the meteoIO-library accessible from python-scripts. 

The wrapper works on Linux, Windows and Mac.

HOW TO USE THE WRAPPER:
After successful installation (see below) a good starting point is the test-script "test.py". This script does the same as the data_converter (which can be found in \doc\examples). It reads in, filters and writes out a data file according to the ini-file ('io.ini'). Also it is shown how to get a time series of a certain meteo parameter (for example temperature). If you want to get information about further functionality of the wrapper you have to look at the pyx-files in the subfoler 'src'. 
Not all functionality of meteoio is included in the wrapper. If you need a certain function, it can easily be added to the wrapper by adjusting the pxd- and pyx files. 
 

WINDOWS AND CYTHON 
Building the wrapper on windows is a tough task. I only managed to get this running with one combination of c++-compiler (msvc 8.0 for python 2.7) and python (python 2.7). Unfortunately it is not possible to build a meteoio-dll with msvc 8.0, so all the source-code has to be compiled when building the wrapper. The more elegant way is to link the wrapper to the meteoio-dll (which is done on LINUX and MAC).  
 
Other compilers under WINDOWS:
With MinGW there are known compatibility problems. With Visual Studio 2017 and python 3.7 it should work, but it did not. I managed to build the wrapper, but I got strange crashes when 
applying the wrapper in python. Hopefully in future there is a compiler and python 3.x combination, which builds the meteoio-library and wrapper without problems under WINDOWS.  


GETTING STARTED

WINDOWS:
  
1) Required Tools:
 -CMake (https://cmake.org/)
 -Microsoft Visual C++ Compiler for Python 2.7 (https://www.microsoft.com/en-us/download/details.aspx?id=44266)
 -python 2.7 32-bit with cython:
	-Anaconda3 32-bit (!!!) (download Anaconda3-2018.12-Windows-x86.exe from https://repo.continuum.io/archive/)
	-open "Anaconda Promt" and type:
		- conda create -n py27 python=2.7
		- activate py27
		- conda install cython
		- deactivate
		
2) Building the wrapper:
 -run CMake and configure + generate with the default settings (this step is necessary to generate the file IOHandler.cc from IOHandler.cmake.cc)
 -run setupWin.py with the argument "install_lib"
 "python.exe setupWin.py install_lib"
 (if you have several python versions installed make sure that "python 2.7 32-bit with cython" is used)
  => now the file "meteoioWrapper.*.pyd" should have been created and copied to the subfolder of your python-installation ..\Lib\site-packages. The pyd-file is the wrapper-module which can be imported in python (e.g. "import meteoioWrapper as mio", see test.py).

 
3) Using the wrapper

run the test-script:
"python.exe test.py"
(if you have several python versions installed make sure that "python 2.7 32-bit" is used)


LINUX (Ubuntu):

1) Required Tools:
-gcc:
	"$sudo apt install build-essential"
	(maybe you first need: "$sudo apt update" and "$sudo apt upgrade")

-svn:
	"$sudo apt install subversion"
-cmake:
	"$sudo apt install cmake-curses-gui"

- meteoio: (see here: https://models.slf.ch/p/meteoio/page/Compiling-MeteoIO )
make sure that meteoio is installed correctly ("make install")!!!

- python: python 2.7 or 3.x should be fine

- cython (version > 0.27): install with anaconda or pip
	"pip install cython" or "pip3 install cython"

 
2) Building the wrapper
 go to the subfolder where the setupPosix.py-file is located (../tools/pythonWrapper/) and type:
 "python setupPosix.py build_ext --inplace" 
 
 Alternatively, if you have sudo-rights, you can build and install the wrapper (this is not tested under Linux yet):
 "sudo python setupPosix.py install_lib" 

3) Using the wrapper
run the test-script:
"python test.py"
 

MAC
1) Required Tools:

- meteoio and c-compiler: (see here: https://models.slf.ch/p/meteoio/page/Compiling-MeteoIO )
make sure that meteoio is installed correctly ("make install")!!!

- python: download and install latest release from https://www.python.org/downloads/mac-osx/

-cython: "sudo pip3 install cython"

2) Building the wrapper
 go to the subfolder where the setupPosix.py-file is located (../tools/pythonWrapper/) and type:
 "python3 setupPosix.py build_ext --inplace" 
 This will build the wrapper in the current directoy, alternatively you can build and install the wrapper (this is not tested under MAC yet):
 "python setupPosix.py install_lib" 

3) Using the wrapper
run the test-script:
"python test.py"


KNOWN PROBLEMS 

Strings (only python 3.x):
There is some issue with different handling of strings between cython, c++ and python. Therefore a python string which is passed to a c++-function has to be converted to utf-8: (e.g. cfg = mio.PyConfig("io.ini".encode('utf-8'), see test.py, alternatively you can write: 'cfg = mio.PyConfig(b"io.ini")'

Memory (only on WINDOWS, maybe because on windows only the 32bit-version works):
-Memory problems with large input files:
  -in IOManager in function getMeteoDataRange() when reading in large files (> 500000 values). 
  -in  IOManager in function writeMeteoData() for large files (> 300000 values). Error: "MemoryError: bad allocation"

