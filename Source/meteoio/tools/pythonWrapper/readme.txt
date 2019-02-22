Date: 13.02.2019
Author: Thiemo Theile


This python-wrapper uses Cython (https://cython.org/) to make the meteoIO-library accessible from python-scripts. 

The wrapper works on Linux, Windows and Mac (???to be tested!!!).

WINDOWS AND CYTHON 
Building the wrapper on windows is a tough task. I only managed to get this running with one combination of c++-compiler (msvc 8.0 for python 2.7) and python (python 2.7). Unfortunately it is not possible to build a meteoio-dll with msvc 8.0, so all the source-code has to be compiled when building the wrapper. The more elegant way is to link the wrapper to the meteoio-dll (which is done on LINUX).  
 
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
 -run setupWin.py with the arguments "build_ext --inplace"
 "python.exe setupWin.py build_ext --inplace"
 (if you have several python versions installed make sure that "python 2.7 32-bit with cython" is used)
 
 => now the file "meteoioWrapper.*.pyd" should have been created. This is the wrapper-module which can be imported in python (e.g. "import meteoioWrapper as mio", see test.py).
 
3) Using the wrapper

run the test-script:
"python.exe test.py"
(if you have several python versions installed make sure that "python 2.7 32-bit with cython" is used)


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

- python:

- cython:
 
2) Building the wrapper
 go to the subfolder where the setup.py-file is located (../tools/pythonWrapper/) and type:
 "python.exe setupLinux.py build_ext --inplace" 

3) Using the wrapper
run the test-script:
"python.exe test.py"
 