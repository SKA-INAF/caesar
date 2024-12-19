# Check if LOG4CXX_DIR is defined by the user
set(LOG4CXX_DIR $ENV{LOG4CXX_DIR})

MESSAGE(STATUS "LOG4CXX_DIR: ${LOG4CXX_DIR}")

# Search for Log4Cxx headers
find_path(LOG4CXX_INCLUDE_DIR
  NAMES log4cxx/log4cxx.h
  HINTS
    ${LOG4CXX_DIR}/include
  PATHS
    /usr/include
    /usr/local/include
    /opt/local/include
)

# Search for libraries
find_library(LOG4CXX_LIBRARY
  NAMES log4cxx
  HINTS
    ${LOG4CXX_DIR}/lib
    ${LOG4CXX_DIR}/lib64
  PATHS
    /usr/lib
    /usr/local/lib
    /usr/lib64
    /usr/local/lib64
)

# Check if both include and library are found
if(LOG4CXX_INCLUDE_DIR AND LOG4CXX_LIBRARY)
	set(LOG4CXX_FOUND TRUE)
	set(LOG4CXX_LIBRARIES ${LOG4CXX_LIBRARY})
	set(LOG4CXX_INCLUDE_DIRS ${LOG4CXX_INCLUDE_DIR})
else()
	set(LOG4CXX_FOUND FALSE)
endif()

# Provide results to the calling script
if(LOG4CXX_FOUND)
	message(STATUS "Found Log4Cxx: ${LOG4CXX_LIBRARY}, includes at ${LOG4CXX_INCLUDE_DIR}")
else()
	message(WARNING "Log4Cxx not found. Please set LOG4CXX_DIR to the installation path.")
endif()

# Make variables available to parent scope
INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Log4Cxx DEFAULT_MSG LOG4CXX_LIBRARIES LOG4CXX_INCLUDE_DIRS)
MARK_AS_ADVANCED(LOG4CXX_LIBRARIES LOG4CXX_INCLUDE_DIRS)
