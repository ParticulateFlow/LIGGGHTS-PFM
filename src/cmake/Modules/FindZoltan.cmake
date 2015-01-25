# - Find Zoltan
# Find the native Zoltan includes and library
# This module defines
#  ZOLTAN_INCLUDE_DIR, where to find zoltan.h, etc.
#  ZOLTAN_LIBRARIES, the libraries needed to use Zoltan.
#  ZOLTAN_FOUND, If false, do not try to use Zoltan.

find_path(ZOLTAN_INCLUDE_DIR zoltan.h)

set(ZOLTAN_NAMES ${ZOLTAN_NAMES} zoltan)
find_library(ZOLTAN_LIBRARY NAMES ${ZOLTAN_NAMES} )

# handle the QUIETLY and REQUIRED arguments and set ZOLTAN_FOUND to TRUE if
# all listed variables are TRUE
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ZOLTAN DEFAULT_MSG ZOLTAN_LIBRARY ZOLTAN_INCLUDE_DIR)

if(ZOLTAN_FOUND)
  set(ZOLTAN_LIBRARIES ${ZOLTAN_LIBRARY})
endif()
