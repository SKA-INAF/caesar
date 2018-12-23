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
find_package(OpenCV REQUIRED)
message (STATUS "OPENCV HEADERS: ${OpenCV_INCLUDE_DIRS}, LIBS: ${OpenCV_LIBS}, ${OpenCV_LIB_COMPONENTS}")
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

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/cmake)
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/cmake/modules)

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
SET (ROOT_REQUIRED_MODULES MathCore RIO Hist Tree Net PyROOT Minuit MathMore)
SET (ROOT_OPTIONAL_MODULES Minuit2 FITSIO)
if(ENABLE_R)
	list(APPEND ROOT_REQUIRED_MODULES RInterface Rtools)
else()
	#list(APPEND ROOT_OPTIONAL_MODULES RInterface Rtools) ## Do not link with ROOT-R if disabled
endif()
 
find_package(ROOT REQUIRED MODULE COMPONENTS ${ROOT_REQUIRED_MODULES} OPTIONAL_COMPONENTS ${ROOT_OPTIONAL_MODULES})

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
MESSAGE(STATUS "Looking for Log4Cxx")
FIND_PACKAGE(Log4Cxx REQUIRED)
MESSAGE(STATUS "LOG4CXX_INCLUDE_DIR: ${LOG4CXX_INCLUDE_DIRS}")
MESSAGE(STATUS "LOG4CXX_LIBRARIES: ${LOG4CXX_LIBRARIES}")


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


