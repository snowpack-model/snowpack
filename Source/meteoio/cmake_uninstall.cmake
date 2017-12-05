IF(NOT EXISTS "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/meteoio/install_manifest.txt")
  MESSAGE(FATAL_ERROR "Cannot find install manifest: \"/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/meteoio/install_manifest.txt\"")
ENDIF(NOT EXISTS "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/meteoio/install_manifest.txt")

FILE(READ "/media/nander/46d45ce3-2115-47b4-90e0-f01b0b010844/src/polarsnowpack/snowpack/Source/meteoio/install_manifest.txt" files)
STRING(REGEX REPLACE "\n" ";" files "${files}")
FOREACH(file ${files})
  MESSAGE(STATUS "Uninstalling \"$ENV{DESTDIR}${file}\"")
  IF(EXISTS "$ENV{DESTDIR}${file}")
    EXEC_PROGRAM(
      "/usr/bin/cmake" ARGS "-E remove \"$ENV{DESTDIR}${file}\""
      OUTPUT_VARIABLE rm_out
      RETURN_VALUE rm_retval
      )
    IF(NOT "${rm_retval}" STREQUAL 0)
      MESSAGE(FATAL_ERROR "Problem when removing \"$ENV{DESTDIR}${file}\"")
    ENDIF(NOT "${rm_retval}" STREQUAL 0)
  ELSE(EXISTS "$ENV{DESTDIR}${file}")
    MESSAGE(STATUS "File \"$ENV{DESTDIR}${file}\" does not exist.")
  ENDIF(EXISTS "$ENV{DESTDIR}${file}")
ENDFOREACH(file)
