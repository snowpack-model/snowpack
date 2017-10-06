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
ADD_SUBDIRECTORY(res1exp)
ADD_SUBDIRECTORY(res5exp)
ADD_SUBDIRECTORY(basics)
ADD_SUBDIRECTORY(mass_and_energy_balance)
ADD_SUBDIRECTORY(linearsolver)
ADD_SUBDIRECTORY(implicitsolver)
ADD_SUBDIRECTORY(albedo)

