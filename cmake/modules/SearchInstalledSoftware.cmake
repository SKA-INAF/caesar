#---Check for installed packages depending on the build options/components eamnbled -
#include(ExternalProject)
#include(FindPackageHandleStandardArgs)


#===============================
#==    Check for BOOST       ===
#===============================
find_package(Boost REQUIRED COMPONENTS filesystem system regex thread log log_setup)
add_definitions(-DBOOST_LOG_USE_NATIVE_SYSLOG -DBOOST_LOG_DYN_LINK)
message (STATUS "BOOST HEADERS: ${Boost_INCLUDE_DIRS}, LIBS: ${Boost_LIBRARIES}")
#-------------------------

#=================================
#==    Check for OPENCV        ===
#=================================
# NB: C API was deprecated in OpenCV >4, so first determine version and skip cmake find accordingly
set(OpenCV_DIR "$ENV{OPENCV_DIR}")
message(STATUS "OpenCV_DIR: " "${OpenCV_DIR}")

if(NOT EXISTS ${OPENCV_VERSION_FILE_DIR}) #may alreay be in cache
	find_path(OPENCV_VERSION_FILE_DIR "version.hpp" 
		PATHS "${OpenCV_DIR}" 
		PATH_SUFFIXES "include" "include/opencv" "include/opencv2" "include/opencv2/core" "include/opencv4/opencv2/core" 
		DOC ""
		NO_DEFAULT_PATH
	)
	if(NOT EXISTS ${OPENCV_VERSION_FILE_DIR})
		message(FATAL_ERROR "OpenCV version file not found")
	else()
		set(OPENCV_VERSION_FILE ${OPENCV_VERSION_FILE_DIR}/version.hpp)
	endif()
else()	
	set(OPENCV_VERSION_FILE ${OPENCV_VERSION_FILE_DIR}/version.hpp)
endif()

message(STATUS "OPENCV_VERSION_FILE: " ${OPENCV_VERSION_FILE})

#file(STRINGS ${OPENCV_VERSION_FILE} OpenCV_VERSIONS_TMP REGEX "^#define CV_[A-Z]+_VERSION[ \t]+[0-9]+$")
file(STRINGS ${OPENCV_VERSION_FILE} OpenCV_VERSIONS_TMP REGEX "^#define CV_VERSION_[A-Z]+[ \t]+[0-9]+$")
message(STATUS "OpenCV_VERSIONS_TMP: " ${OpenCV_VERSIONS_TMP})
string(REGEX REPLACE ".*#define CV_VERSION_MAJOR[ \t]+([0-9]+).*" "\\1" OpenCV_VERSION_MAJOR "${OpenCV_VERSIONS_TMP}")
string(REGEX REPLACE ".*#define CV_VERSION_MINOR[ \t]+([0-9]+).*" "\\1" OpenCV_VERSION_MINOR "${OpenCV_VERSIONS_TMP}")
string(REGEX REPLACE ".*#define CV_VERSION_REVISION[ \t]+([0-9]+).*" "\\1" OpenCV_VERSION_REVISION "${OpenCV_VERSIONS_TMP}")

message(STATUS "OpenCV_VERSION_MAJOR: " ${OpenCV_VERSION_MAJOR})
message(STATUS "OpenCV_VERSION_MINOR: " ${OpenCV_VERSION_MINOR})
message(STATUS "OpenCV_VERSION_REVISION: " ${OpenCV_VERSION_REVISION})

set(OpenCV_VERSION ${OpenCV_VERSION_MAJOR}.${OpenCV_VERSION_MINOR}.${OpenCV_VERSION_REVISION} CACHE STRING "" FORCE)

message(STATUS "OpenCV_VERSION: " ${OpenCV_VERSION})

#if(${OpenCV_VERSION} VERSION_LESS 4.0.0)
#	find_package(OpenCV REQUIRED)
#	message (STATUS "OPENCV HEADERS: ${OpenCV_INCLUDE_DIRS}, LIBS: ${OpenCV_LIBS}, ${OpenCV_LIB_COMPONENTS}")
#endif()

find_package(OpenCV2 REQUIRED)
message (STATUS "OPENCV2 HEADERS: ${OpenCV2_INCLUDE_DIRS}, LIBS: ${OpenCV2_LIBRARIES}")
#-------------------------

#=================================
#==   Check for R PROJECT      ===
#=================================

option(ENABLE_R "Enable building of R components" OFF)
if(ENABLE_R)
	MESSAGE(STATUS "Looking for R")
	##find_package (R REQUIRED COMPONENTS base RInside Rcpp rrcovHD truncnorm FNN akima)
	find_package (R REQUIRED COMPONENTS base RInside Rcpp)	
	if(R_FOUND)
		MESSAGE(STATUS "R found, defining preprocessor flag R_ENABLED")
		add_definitions(-DR_ENABLED=1)
		message (STATUS "R_INCLUDE_DIR: ${R_INCLUDE_DIR}")
		message (STATUS "R_LIBRARIES: ${R_LIBRARIES}")
	else()
		MESSAGE(SEND_ERROR "R not found!")
	endif()
endif()


#message (STATUS "Looking for R")
#find_package (R REQUIRED COMPONENTS base RInside Rcpp rrcovHD truncnorm FNN akima)
#message (STATUS "R_INCLUDE_DIR: ${R_INCLUDE_DIR}")
#message (STATUS "R_LIBRARIES: ${R_LIBRARIES}")
#-------------------------


#==================================
#==    Check for ROOT           ===
#==================================
message(STATUS "Looking for ROOT")

if(NOT DEFINED ENV{ROOTSYS})
	message(SEND_ERROR "ROOTSYS variable not defined!")
endif()

#list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}/cmake)
#list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
#list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/cmake)
#list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/cmake/modules)


# - Include ROOTConfig.cmake
SET (ROOT_CMAKE_CONFIG_FILE $ENV{ROOTSYS}/cmake/ROOTConfig.cmake)
if(EXISTS ${ROOT_CMAKE_CONFIG_FILE})
	MESSAGE(STATUS "Including ROOT config file ${ROOT_CMAKE_CONFIG_FILE}")
  #include(${ROOT_CMAKE_CONFIG_FILE})
	list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

else()
	MESSAGE(STATUS "ROOTConfig.cmake was not found, trying to use FindROOT.cmake ...")
	list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
endif()

	
# - Determine ROOT python version and pyroot components
find_program(ROOT_CONFIG_EXECUTABLE NAMES root-config)

#execute_process( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --python-version
execute_process( COMMAND bash "-c" "${ROOT_CONFIG_EXECUTABLE} --python-version"
                 OUTPUT_VARIABLE PYVERSION_ROOT
 	               OUTPUT_STRIP_TRAILING_WHITESPACE )


if ("${PYVERSION_ROOT}" STREQUAL "") # root-config --python is empty, try a workaround (parse --config python lib)
	execute_process( COMMAND bash "-c" "${ROOT_CONFIG_EXECUTABLE} --config | grep -Po '^.*?\\K(?<=PYTHON_LIBRARY_RELEASE=).*?(?=.so).so'"
                   OUTPUT_VARIABLE ROOT_PYTHON_LIB_PATH
                   OUTPUT_STRIP_TRAILING_WHITESPACE ) 	

	MESSAGE(STATUS "ROOT_PYTHON_LIB_PATH: ${ROOT_PYTHON_LIB_PATH}")
	string(REGEX MATCH "([0-9]+)\\.([0-9]+)" ROOT_PYTHON_LIB_VERSION ${ROOT_PYTHON_LIB_PATH})
	MESSAGE(STATUS "ROOT_PYTHON_LIB_VERSION: ${ROOT_PYTHON_LIB_VERSION}")
	if ("${ROOT_PYTHON_LIB_VERSION}" STREQUAL "") 
		MESSAGE(SEND_ERROR "Cannot parse python version used in ROOT building!")
	endif()
	string(REPLACE "." ";" VERSION_LIST ${ROOT_PYTHON_LIB_VERSION})
else()
	MESSAGE(STATUS "PYVERSION_ROOT: " ${PYVERSION_ROOT})
	string(REPLACE "." ";" VERSION_LIST ${PYVERSION_ROOT})
endif()

list(GET VERSION_LIST 0 PYVERSION_ROOT_MAJOR)
list(GET VERSION_LIST 1 PYVERSION_ROOT_MINOR)
#list(GET VERSION_LIST 2 PYVERSION_ROOT_PATCH)
MESSAGE(STATUS "PYVERSION_ROOT MAJOR: " ${PYVERSION_ROOT_MAJOR})
MESSAGE(STATUS "PYVERSION_ROOT MINOR: " ${PYVERSION_ROOT_MINOR})
#MESSAGE(STATUS "PYVERSION_ROOT PATCH: " ${PYVERSION_ROOT_PATCH})
SET(PYROOT_COMPONENTS "ROOTPythonizations${PYVERSION_ROOT_MAJOR}_${PYVERSION_ROOT_MINOR}")
#list(APPEND PYROOT_COMPONENTS "ROOTPythonizations${PYVERSION_ROOT_MAJOR}_${PYVERSION_ROOT_MINOR}" "cppyy${PYVERSION_ROOT_MAJOR}_${PYVERSION_ROOT_MINOR}" "cppyy_backend${PYVERSION_ROOT_MAJOR}_${PYVERSION_ROOT_MINOR}")
MESSAGE(STATUS "PYROOT_COMPONENTS: " ${PYROOT_COMPONENTS})


#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
SET (ROOT_REQUIRED_MODULES MathCore RIO Hist Tree Net Minuit MathMore TMVA)
SET (ROOT_OPTIONAL_MODULES Minuit2 FITSIO ${PYROOT_COMPONENTS})
#SET (ROOT_OPTIONAL_MODULES Minuit2 FITSIO PyROOT)
if(ENABLE_R)
	list(APPEND ROOT_REQUIRED_MODULES RInterface Rtools)
else()
	#list(APPEND ROOT_OPTIONAL_MODULES RInterface Rtools) ## Do not link with ROOT-R if disabled
endif()
 
find_package(ROOT CONFIG REQUIRED)
#find_package(ROOT REQUIRED MODULE COMPONENTS ${ROOT_REQUIRED_MODULES} OPTIONAL_COMPONENTS ${ROOT_OPTIONAL_MODULES})
find_package(ROOT REQUIRED COMPONENTS ${ROOT_REQUIRED_MODULES} OPTIONAL_COMPONENTS ${ROOT_OPTIONAL_MODULES})

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE}) 

message (STATUS "ROOT HEADERS: ${ROOT_INCLUDE_DIRS}, LIBS: ${ROOT_LIBRARIES}, REQUIRED_MODULES: ${ROOT_REQUIRED_MODULES}, OPTIONAL_MODULES: ${ROOT_OPTIONAL_MODULES} ")
message (STATUS "MINUIT2 FOUND? ${ROOT_minuit2_FOUND}")
if(${ROOT_minuit2_FOUND})
	MESSAGE(STATUS "MINUIT2 found, defining preprocessor flag MINUIT2_ENABLED=true and adding library ${ROOT_minuit2_LIBRARY} to linked lib list")
	add_definitions(-DMINUIT2_ENABLED=1)
else()
	MESSAGE(STATUS "MINUIT2 not found")
endif()


if(ENABLE_R)
  message (STATUS "ROOT-R FOUND? ${ROOT_r_FOUND}")
  if(${ROOT_r_FOUND})
		MESSAGE(STATUS "ROOT-R module found, defining preprocessor flag ROOTR_ENABLED=true and adding library ${ROOT_r_LIBRARY} to linked lib list")
		add_definitions(-DROOTR_ENABLED=1)
  else()
    MESSAGE(STATUS "ROOT-R module not found")
	endif()
else()
	MESSAGE(STATUS "ROOT-R module disabled")
endif()
#--------------------------------


#==================================
#==    Check for GSL library    ===
#==================================
##message (STATUS "Looking for GSL library")
##find_package (GSL REQUIRED)
##find_package (GSL)
##message (STATUS "GSL_INCLUDE_DIRS: ${GSL_INCLUDE_DIRS}")
##message (STATUS "GSL_LIBRARIES: ${GSL_LIBRARIES}")

#==================================
#==    Check for Python         ===
#==================================
##message (STATUS "Looking for Python")
##find_package (PythonLibs REQUIRED)
##find_package (PythonLibs)
##message (STATUS "PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}, PYTHON_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}")

#==================================
#== Check for Python Modules    ===
#==================================
#include(FindPythonModule)

# pyfits
#message (STATUS "Looking for Python pyfits module")
##find_python_module(pyfits REQUIRED)
#find_python_module(pyfits)
# astropy
#message (STATUS "Looking for Python astropy module")
##find_python_module(astropy REQUIRED)
#find_python_module(astropy)


#==================================
#==   Check for Log4Cxx         ===
#==================================
option(ENABLE_LOGGING "Enable Logging" ON)
if(ENABLE_LOGGING)
	MESSAGE(STATUS "Looking for Log4Cxx")
	FIND_PACKAGE(Log4Cxx REQUIRED)
	add_definitions(-DLOGGING_ENABLED=1)
	MESSAGE(STATUS "LOG4CXX_INCLUDE_DIR: ${LOG4CXX_INCLUDE_DIRS}")
	MESSAGE(STATUS "LOG4CXX_LIBRARIES: ${LOG4CXX_LIBRARIES}")
endif()

#==================================
#==   Check for CFITSIO         ===
#==================================
MESSAGE(STATUS "Looking for CFITSIO")
FIND_PACKAGE(CFITSIO REQUIRED)
if (NOT CFITSIO_FOUND)
	MESSAGE(SEND_ERROR "CFITSIO library not found!")
endif()
MESSAGE(STATUS "CFITSIO_INCLUDE_DIR: ${CFITSIO_INCLUDE_DIR}")
MESSAGE(STATUS "CFITSIO_LIBRARIES: ${CFITSIO_LIBRARIES}")

#==================================
#==   Check for JSONCPP         ===
#==================================
MESSAGE(STATUS "Looking for JSONCPP")

if(NOT DEFINED ENV{JSONCPP_ROOT})
	MESSAGE(SEND_ERROR "JSONCPP_ROOT variable not defined!")
endif()

SET (JSONCPP_ROOT $ENV{JSONCPP_ROOT})
MESSAGE(STATUS "JSONCPP_ROOT: ${JSONCPP_ROOT}")

FIND_PATH (JSONCPP_INCLUDE_DIR
	NAMES json/json.h
 	HINTS
	${JSONCPP_ROOT}/include
)

FIND_LIBRARY (JSONCPP_LIBRARIES NAMES jsoncpp HINTS ${JSONCPP_ROOT}/lib)
MARK_AS_ADVANCED (JSONCPP_INCLUDE_DIR JSONCPP_LIBRARIES)
#MESSAGE(STATUS "JSONCPP_INCLUDE_DIR: ${JSONCPP_INCLUDE_DIR}")
#MESSAGE(STATUS "JSONCPP_LIBRARIES: ${JSONCPP_LIBRARIES}")


#=================================
#==   Check for PROTOBUF       ===
#=================================

MESSAGE(STATUS "Looking for Protobuf lib")
FIND_PACKAGE(Protobuf REQUIRED COMPONENTS protobuf protoc)
IF (NOT PROTOBUF_FOUND)
	MESSAGE(SEND_ERROR "Protobuf not found!")
endif()
MESSAGE(STATUS "PROTOBUF_INCLUDE_DIR: ${PROTOBUF_INCLUDE_DIRS}")
MESSAGE(STATUS "PROTOBUF_LIBRARIES: ${PROTOBUF_LIBRARIES}")


#=======================================
#==   SET BUILD APPS OPTION          ===
#=======================================
option(BUILD_APPS "Enable building of Caesar applications in apps dir" ON)
if(BUILD_APPS)
	add_definitions(-DBUILD_CAESAR_APPS)
endif()

#==============================================
#==   ENABLE OPENMP OPTION?                 ===
#==============================================
##option(BUILD_WITH_OPENMP "Enable building of Caesar with OpenMP" OFF)
##if(BUILD_WITH_OPENMP)
option(ENABLE_OPENMP "Enable building of Caesar with OpenMP" OFF)
if(ENABLE_OPENMP)
	MESSAGE(STATUS "Looking for OpenMP")
	include(FindOpenMP)
	if(OPENMP_FOUND)
		MESSAGE(STATUS "OpenMP found, defining preprocessor flag OPENMP_ENABLED")
		add_definitions(-DOPENMP_ENABLED=1)
	else()
		MESSAGE(SEND_ERROR "OpenMP not found!")
	endif()
endif()


#==============================================
#==   ENABLE MPI OPTION?                    ===
#==============================================
option(ENABLE_MPI "Enable building of Caesar MPI components" OFF)
if(ENABLE_MPI)
	MESSAGE(STATUS "Looking for MPI")

	find_package(MPI REQUIRED)
	if(MPI_FOUND)
		MESSAGE(STATUS "MPI found, defining preprocessor flag MPI_ENABLED")
		add_definitions(-DMPI_ENABLED=1)
	else()
		MESSAGE(SEND_ERROR "MPI not found!")
	endif()
		
	#MESSAGE(STATUS "MPI_CXX_INCLUDE_PATH: ${MPI_CXX_INCLUDE_PATH}")
	#MESSAGE(STATUS "MPI_LIBRARIES: ${MPI_LIBRARIES}")
	#MESSAGE(STATUS "MPI_CXX_LIBRARIES: ${MPI_CXX_LIBRARIES}")
endif()

#======================================
#==   Check for Google Test         ===
#======================================
option(ENABLE_TEST "Enable unit test building" ON)
if (ENABLE_TEST)
	MESSAGE(STATUS "Looking for GoogleTest")
	add_definitions(-DENABLE_TEST)
	enable_testing()
	FIND_PACKAGE(GTest REQUIRED)
	#MESSAGE(STATUS "GTEST_INCLUDE_DIRS: ${GTEST_INCLUDE_DIRS}")
	#MESSAGE(STATUS "GTEST_LIBRARIES: ${GTEST_BOTH_LIBRARIES}")
endif()

#==================================
#==   Check for Doxygen         ===
#==================================
MESSAGE(STATUS "Looking for Doxygen")
find_package(Doxygen)
if (NOT DOXYGEN_FOUND)
	MESSAGE(STATUS "Doxygen not found, cannot generate project documentation!")
endif()
option(BUILD_DOC "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

#==================================
#==   Check for Sphinx          ===
#==================================
MESSAGE(STATUS "Looking for Sphinx")
find_package(Sphinx)
if (NOT SPHINX_FOUND)
	MESSAGE(STATUS "Sphinx not found, cannot generate sphinx project documentation!")
endif()
option(BUILD_WIKI_DOC "Create and install the HTML based API & wiki documentation (requires Sphinx)" ${SPHINX_FOUND})

#======================================
#==   Check for VTK                 ===
#======================================
option(ENABLE_VTK "Enable VTK" OFF)
if (ENABLE_VTK)
	MESSAGE(STATUS "Looking for VTK library")	
	FIND_PACKAGE(VTK REQUIRED)
	if(VTK_FOUND)
		MESSAGE(STATUS "VTK found, defining preprocessor flag VTK_ENABLED")
		add_definitions(-DVTK_ENABLED=1)
		MESSAGE(STATUS "VTK_INCLUDE_DIRS= ${VTK_INCLUDE_DIRS}")
		MESSAGE(STATUS "VTK_LIBRARY_DIRS= ${VTK_LIBRARY_DIRS}")
		MESSAGE(STATUS "VTK_LIBRARIES= ${VTK_LIBRARIES}")

		#SET(VTK_LIBRARIES ${VTK_LIBRARIES}
		#	vtkHybrid
  	#	vtkVolumeRendering
  	#	vtkRendering
  	#	QVTK
		#)

	else()
		MESSAGE(SEND_ERROR "VTK not found!")
	endif()	
endif()

#======================================
#==   Check for WCSLIB              ===
#======================================
#if (NOT "$ENV{WCSLIB_ROOT_DIR}" STREQUAL "")
#  set(WCSLIB_ROOT_DIR $ENV{WCSLIB_ROOT_DIR})
#endif()
#MESSAGE(STATUS WCSLIB_ROOT_DIR ${WCSLIB_ROOT_DIR})
#MESSAGE(STATUS "Looking for WCS library")	
#FIND_PACKAGE(WCSLIB REQUIRED)

#MESSAGE(STATUS "WCSLIB_INCLUDE_DIRS= ${WCSLIB_INCLUDE_DIRS}")
#MESSAGE(STATUS "WCSLIB_LIBRARIES= ${WCSLIB_LIBRARIES}")

#======================================
#==   Check for CASA CORE           ===
#======================================
option(ENABLE_CASACORE "Enable CASACORE" OFF)
if (ENABLE_CASACORE)
	if (NOT "$ENV{CASACORE_ROOT_DIR}" STREQUAL "")
  	set(CASACORE_ROOT_DIR $ENV{CASACORE_ROOT_DIR})
	endif()

	MESSAGE(STATUS "Looking for CASACORE library")	
	FIND_PACKAGE(Casacore REQUIRED)
	if(CASACORE_FOUND)
		MESSAGE(STATUS "CASACORE found, defining preprocessor flag CASACORE_ENABLED")
		add_definitions(-DCASACORE_ENABLED=1)
		MESSAGE(STATUS "CASACORE_INCLUDE_DIRS= ${CASACORE_INCLUDE_DIRS}")
		MESSAGE(STATUS "CASACORE_LIBRARIES= ${CASACORE_LIBRARIES}")

	else()
		MESSAGE(SEND_ERROR "CASACORE not found!")
	endif()	
endif()

#======================================
#==   Check for GIT                 ===
#======================================
MESSAGE(STATUS "Looking for Git")
find_package(Git)
if (NOT Git_FOUND)
	MESSAGE(STATUS "Git not found, cannot generate software version number!")
endif()

