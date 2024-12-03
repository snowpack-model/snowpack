# Updating third party libraries

Some third party libraries are packaged with MeteoIO. This choice is made based on the following conditions:
 * they provide some important functionality (unlike plugins that might be only useful for some special cases)
 * they are small enough to package within MeteoIO
 * they are easy to integrate / compile


## ctsa
CTSA is a C software package for univariate time series analysis. It is available on [github](https://github.com/rafat/ctsa). In order to update it to a new version, here are the steps to follow:
 * copy the content of src/*.c and src/*.h from the ctsa repository into thirdParty/ctsa
 * copy the file *ctsa.h* at the root of the ctsa repository into thirdParty
 * re-add the following line after the `typedef struct auto_arima_set* auto_arima_object;`line (around line 16) in ctsa.h: `auto_arima_object auto_arima_copy(auto_arima_object original);`
 * check that no new file has been added into the ctsa repository. If necessary, add any new c files into thirdParty/CMakeLists.txt
 * make sure that our file *thirdParty/ctsa/ctsa_wrapper.c* is still there...


## tinyexpr
TinyExpr is a very small parser and evaluation library for evaluating math expressions from C. It is available on it [home page](https://codeplea.com/tinyexpr) and its repository is on [github](https://github.com/codeplea/tinyexpr). In order to update it:
 * copy tinyexpr.c and tinyexpr.h into *thirdParty/*


## picojson
Picojson is a header-file-only, JSON parser serializer in C++. It is available on [github](https://github.com/kazuho/picojson). As this is a header-file-only library, updating it is only a matter of copying the file *picojson.h* from its repository into *thirdParty/*.
