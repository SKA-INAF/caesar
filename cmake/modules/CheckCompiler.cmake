#---------------------------------------------------------------------------------------------------
#  CheckCompiler.cmake
#---------------------------------------------------------------------------------------------------

#----Get the compiler file name (to ensure re-location)---------------------------------------------
get_filename_component(_compiler_name ${CMAKE_CXX_COMPILER} NAME)
get_filename_component(_compiler_path ${CMAKE_CXX_COMPILER} PATH)
if("$ENV{PATH}" MATCHES ${_compiler_path})
  set(CXX ${_compiler_name})
else()
  set(CXX ${CMAKE_CXX_COMPILER})
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" uppercase_CMAKE_BUILD_TYPE)

#----Test if clang setup works----------------------------------------------------------------------
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  exec_program(${CMAKE_CXX_COMPILER} ARGS "--version 2>&1 | grep version" OUTPUT_VARIABLE _clang_version_info)
  string(REGEX REPLACE "^.*[ ]version[ ]([0-9]+)\\.[0-9]+.*" "\\1" CLANG_MAJOR "${_clang_version_info}")
  string(REGEX REPLACE "^.*[ ]version[ ][0-9]+\\.([0-9]+).*" "\\1" CLANG_MINOR "${_clang_version_info}")
  message(STATUS "Found Clang. Major version ${CLANG_MAJOR}, minor version ${CLANG_MINOR}")
  set(COMPILER_VERSION clang${CLANG_MAJOR}${CLANG_MINOR})
else()
  set(CLANG_MAJOR 0)
  set(CLANG_MINOR 0)
endif()

#---Obtain the major and minor version of the GNU compiler-------------------------------------------
if (CMAKE_COMPILER_IS_GNUCXX)
  exec_program(${CMAKE_C_COMPILER} ARGS "-dumpversion" OUTPUT_VARIABLE _gcc_version_info)
  string(REGEX REPLACE "^([0-9]+).*$"                   "\\1" GCC_MAJOR ${_gcc_version_info})
  string(REGEX REPLACE "^[0-9]+\\.([0-9]+).*$"          "\\1" GCC_MINOR ${_gcc_version_info})
  string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1" GCC_PATCH ${_gcc_version_info})

  if(GCC_PATCH MATCHES "\\.+")
    set(GCC_PATCH "")
  endif()
  if(GCC_MINOR MATCHES "\\.+")
    set(GCC_MINOR "")
  endif()
  if(GCC_MAJOR MATCHES "\\.+")
    set(GCC_MAJOR "")
  endif()
  message(STATUS "Found GCC. Major version ${GCC_MAJOR}, minor version ${GCC_MINOR}")
  set(COMPILER_VERSION gcc${GCC_MAJOR}${GCC_MINOR}${GCC_PATCH})
else()
  set(GCC_MAJOR 0)
  set(GCC_MINOR 0)
endif()

#---Set a default build type for single-configuration CMake generators if no build type is set------
set(CMAKE_CONFIGURATION_TYPES Release MinSizeRel Debug RelWithDebInfo)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are: Release, MinSizeRel, Debug, RelWithDebInfo." FORCE)
endif()
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

include(CheckCXXCompilerFlag)
include(CheckCCompilerFlag)

#---Check for cxx11 option------------------------------------------------------------
if(cxx11 AND cxx14)
  message(STATUS "c++11 mode requested but superseded by request for c++14 mode")
  set(cxx11 OFF CACHE BOOL "" FORCE)
endif()
if(cxx11)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" HAS_CXX11)
  if(NOT HAS_CXX11)
    message(STATUS "Current compiler does not suppport -std=c++11 option. Switching OFF cxx11 option")
    set(cxx11 OFF CACHE BOOL "" FORCE)
  endif()
endif()
if(cxx14)
  CHECK_CXX_COMPILER_FLAG("-std=c++14" HAS_CXX14)
  if(NOT HAS_CXX14)
    message(STATUS "Current compiler does not suppport -std=c++14 option. Switching OFF cxx14 option")
    set(cxx14 OFF CACHE BOOL "" FORCE)
  endif()
endif()

#---Check for libcxx option------------------------------------------------------------
if(libcxx)
  CHECK_CXX_COMPILER_FLAG("-stdlib=libc++" HAS_LIBCXX11)
  if(NOT HAS_LIBCXX11)
    message(STATUS "Current compiler does not suppport -stdlib=libc++ option. Switching OFF libcxx option")
    set(libcxx OFF CACHE BOOL "" FORCE)
  endif()
endif()


#---Need to locate thead libraries and options to set properly some compilation flags----------------
find_package(Threads)
if(CMAKE_USE_PTHREADS_INIT)
  set(CMAKE_THREAD_FLAG -pthread)
else()
  set(CMAKE_THREAD_FLAG)
endif()

#---Setup details depending opn the major platform type----------------------------------------------
if(CMAKE_SYSTEM_NAME MATCHES Linux)
  include(SetUpLinux)
elseif(APPLE)
  include(SetUpMacOS)
elseif(WIN32)
  include(SetupWindows)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_THREAD_FLAG}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_THREAD_FLAG}")

if(cxx11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
  endif()
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
  endif()
endif()

if(cxx14)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
  if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
  endif()
endif()

if(libcxx)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

if(gcctoolchain)
  CHECK_CXX_COMPILER_FLAG("--gcc-toolchain=${gcctoolchain}" HAS_GCCTOOLCHAIN)
  if(HAS_GCCTOOLCHAIN)
     set(CMAKE_CXX_FLAGS "--gcc-toolchain=${gcctoolchain} ${CMAKE_CXX_FLAGS}")
  endif()
endif()

if(gnuinstall)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DR__HAVE_CONFIG")
endif()

#---Print the final compiler flags--------------------------------------------------------------------
message(STATUS "SFINDER Platform: ${ROOT_PLATFORM}")
message(STATUS "SFINDER Architecture: ${ROOT_ARCHITECTURE}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler Flags: ${CMAKE_CXX_FLAGS} ${ALL_CXX_FLAGS_${uppercase_CMAKE_BUILD_TYPE}}")
