#############################
## Global Configuration    ##
#############################

# set Memory check command
set(CTEST_MEMORYCHECK_COMMAND valgrind)

#get all include flags needed to make executables
SET(SHAREDNAME ${PROJECT_NAME})
SET(LIBRARIES ${SHAREDNAME} ${CMAKE_DL_LIBS})

# go back to the source dir to have all .h files as would be installed
INCLUDE_DIRECTORIES(..)

###################
## Add Tests     ##
###################
ADD_SUBDIRECTORY(simple)
ADD_SUBDIRECTORY(basics)
