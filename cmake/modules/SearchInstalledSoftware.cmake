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
SET (ROOT_REQUIRED_MODULES MathCore RIO Hist Tree Net PyROOT)
if(ENABLE_R)
	list(APPEND ROOT_REQUIRED_MODULES RInterface)
endif()

find_package(ROOT REQUIRED MODULE COMPONENTS ${ROOT_REQUIRED_MODULES})

#find_package(ROOT REQUIRED MODULE COMPONENTS MathCore RIO Hist Tree Net FITSIO PyROOT RInterface)
#find_package(ROOT REQUIRED)
#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE}) ### NOT WORKING!!
#include(ROOTUseFile) 
#include_directories(${ROOT_INCLUDE_DIR})

message (STATUS "ROOT HEADERS: ${ROOT_INCLUDE_DIRS}, LIBS: ${ROOT_LIBRARIES}, REQUIRED_MODULES: ${ROOT_REQUIRED_MODULES}")
#--------------------------------


#==================================
#==    Check for GSL library    ===
#==================================
message (STATUS "Looking for GSL library...")
#include(FindGSL REQUIRED)
find_package (GSL REQUIRED)
message (STATUS "GSL_INCLUDE_DIRS: ${GSL_INCLUDE_DIRS}")
message (STATUS "GSL_LIBRARIES: ${GSL_LIBRARIES}")

#==================================
#==    Check for Python         ===
#==================================
message (STATUS "Looking for Python")
find_package (PythonLibs REQUIRED)
message (STATUS "PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}, PYTHON_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}")

#==================================
#== Check for Python Modules    ===
#==================================
include(FindPythonModule)
# pyfits
message (STATUS "Looking for Python pyfits module")
find_python_module(pyfits REQUIRED)
# astropy
message (STATUS "Looking for Python astropy module")
find_python_module(astropy REQUIRED)


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
option(BUILD_WITH_OPENMP "Enable building of Caesar with OpenMP" OFF)
if(BUILD_WITH_OPENMP)
	MESSAGE(STATUS "Looking for OpenMP")
	include(FindOpenMP)
	if(OPENMP_FOUND)
		MESSAGE(STATUS "OpenMP found, defining preprocessor flag OPENMP_ENABLED")
		add_definitions(-DOPENMP_ENABLED=1)
	else()
		MESSAGE(SEND_ERROR "OPENMP not found!")
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


#===========================================
#==   Check for Tango Framework          ===
#===========================================
option(ENABLE_SERVER "Enable Caesar server in building" OFF)
if(ENABLE_SERVER)
	## Define a preprocessor option
	add_definitions(-DBUILD_CAESAR_SERVER)
	
	

	#================================
	#==   Check for MSGPACK       ===
	#================================
	MESSAGE(STATUS "Looking for MessagePack lib")
	FIND_PACKAGE(MsgPack REQUIRED)
	IF (NOT MSGPACK_FOUND)
		MESSAGE(SEND_ERROR "MessagePack not found!")
	endif()
	MESSAGE(STATUS "MSGPACK_INCLUDE_DIR: ${MSGPACK_INCLUDE_DIRS}")
	MESSAGE(STATUS "MSGPACK_LIBRARIES: ${MSGPACK_LIBRARIES}")
	

	#================================
	#==   Check for OMNIORB       ===
	#================================
	MESSAGE(STATUS "Looking for omniORB")
	if(NOT DEFINED ENV{OMNI_ROOT})
		MESSAGE(SEND_ERROR "OMNI_ROOT variable not defined!")
	endif()

	SET (OMNI_ROOT $ENV{OMNI_ROOT})
	MESSAGE(STATUS "OMNI_ROOT: ${OMNI_ROOT}")

	FIND_PATH (OMNIORB_INCLUDE_DIR
		NAMES omniconfig.h
  	HINTS
  	${OMNI_ROOT}/include
	)

	FIND_LIBRARY (OMNIORB_LIB1 NAMES omniORB4 HINTS ${OMNI_ROOT}/lib)
	FIND_LIBRARY (OMNIORB_LIB2 NAMES COS4 HINTS ${OMNI_ROOT}/lib)
	FIND_LIBRARY (OMNIORB_LIB3 NAMES omniDynamic4 HINTS ${OMNI_ROOT}/lib)
	FIND_LIBRARY (OMNIORB_LIB4 NAMES omnithread HINTS ${OMNI_ROOT}/lib)
	list(APPEND OMNIORB_LIBRARIES ${OMNIORB_LIB1} ${OMNIORB_LIB2} ${OMNIORB_LIB3} ${OMNIORB_LIB4})

	MARK_AS_ADVANCED (OMNIORB_INCLUDE_DIR OMNIORB_LIBRARIES)
	MESSAGE(STATUS "OMNIORB_INCLUDE_DIR: ${OMNIORB_INCLUDE_DIR}")
	MESSAGE(STATUS "OMNIORB_LIBRARIES: ${OMNIORB_LIBRARIES}")

	#================================
	#==   Check for ZMQ           ===
	#================================
	MESSAGE(STATUS "Looking for ZMQ")
	if(NOT DEFINED ENV{ZMQ_ROOT})
		MESSAGE(SEND_ERROR "ZMQ_ROOT variable not defined!")
	endif()

	SET (ZMQ_ROOT $ENV{ZMQ_ROOT})
	MESSAGE(STATUS "ZMQ_ROOT: ${ZMQ_ROOT}")

	FIND_PATH (ZMQ_INCLUDE_DIR
		NAMES zmq.h
  	HINTS
  	${ZMQ_ROOT}/include
	)

	FIND_LIBRARY (ZMQ_LIBRARIES NAMES zmq HINTS ${ZMQ_ROOT}/lib)

	MARK_AS_ADVANCED (ZMQ_INCLUDE_DIR ZMQ_LIBRARIES)
	MESSAGE(STATUS "ZMQ_INCLUDE_DIR: ${ZMQ_INCLUDE_DIR}")
	MESSAGE(STATUS "ZMQ_LIBRARIES: ${ZMQ_LIBRARIES}")


	#==================================	
	#==   Check for TANGO           ===
	#==================================
	MESSAGE(STATUS "Looking for TANGO Framework")
	
	## Define Tango preprocessor flag
  add_definitions(-DUSE_TANGO)
	add_definitions(-DAPPENDERS_HAVE_LEVEL_THRESHOLD=1)

	if(NOT DEFINED ENV{TANGO_ROOT})
		MESSAGE(SEND_ERROR "TANGO_ROOT variable not defined!")
	endif()

	SET (TANGO_ROOT $ENV{TANGO_ROOT})
	MESSAGE(STATUS "TANGO_ROOT: ${TANGO_ROOT}")

	FIND_PATH (TANGO_INCLUDE_DIR
		NAMES tango.h
  	HINTS
  	${TANGO_ROOT}/include/tango
	)

	FIND_LIBRARY (TANGO_LIB1 NAMES tango HINTS ${TANGO_ROOT}/lib)
	FIND_LIBRARY (TANGO_LIB2 NAMES log4tango HINTS ${TANGO_ROOT}/lib)
	list(APPEND TANGO_LIBRARIES ${TANGO_LIB1} ${TANGO_LIB2})

	MARK_AS_ADVANCED (TANGO_INCLUDE_DIR TANGO_LIBRARIES)
	MESSAGE(STATUS "TANGO_INCLUDE_DIR: ${TANGO_INCLUDE_DIR}")
	MESSAGE(STATUS "TANGO_LIBRARIES: ${TANGO_LIBRARIES}")


	#=================================
	#==   Check for YAT4           ===
	#=================================
	MESSAGE(STATUS "Looking for YAT Library")
	if(NOT DEFINED ENV{YAT_ROOT})
		MESSAGE(SEND_ERROR "YAT_ROOT variable not defined!")
	endif()

	FIND_PATH (YAT_INCLUDE_DIR
		NAMES yat/threading/Task.h
  	HINTS
  	$ENV{YAT_ROOT}/include
	)
	FIND_LIBRARY (YAT_LIBRARIES 
		NAMES yat
		HINTS 
		$ENV{YAT_ROOT}/lib
	)
	MARK_AS_ADVANCED (YAT_INCLUDE_DIR YAT_LIBRARIES)
	MESSAGE(STATUS "YAT_INCLUDE_DIR: ${YAT_INCLUDE_DIR}")
	MESSAGE(STATUS "YAT_LIBRARIES: ${YAT_LIBRARIES}")
		


	#======================================
	#==   Check for YAT4TANGO           ===
	#======================================
	MESSAGE(STATUS "Looking for YAT4TANGO Library")
	if(NOT DEFINED ENV{YAT4TANGO_ROOT})
		MESSAGE(SEND_ERROR "YAT4TANGO_ROOT variable not defined!")
	endif()

	FIND_PATH (YAT4TANGO_INCLUDE_DIR
		NAMES yat4tango/DeviceTask.h
  	HINTS
  	$ENV{YAT4TANGO_ROOT}/include
	)
	FIND_LIBRARY (YAT4TANGO_LIBRARIES 
		NAMES yat4tango 
		HINTS 
		$ENV{YAT4TANGO_ROOT}/lib
	)
	MARK_AS_ADVANCED (YAT4TANGO_INCLUDE_DIR YAT4TANGO_LIBRARIES)
	MESSAGE(STATUS "YAT4TANGO_INCLUDE_DIR: ${YAT4TANGO_INCLUDE_DIR}")
	MESSAGE(STATUS "YAT4TANGO_LIBRARIES: ${YAT4TANGO_LIBRARIES}")
		

endif() ### close if enable Tango

