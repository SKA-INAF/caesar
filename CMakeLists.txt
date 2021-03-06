#----------------------------------------------------------------
# - Enforce an out-of-source builds before anything else
# ---------------------------------------------------------------
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(STATUS "Caesar requires an out-of-source build.")
    message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
    message(STATUS "CMakeCache.txt")
    message(STATUS "CMakeFiles")
    message(STATUS "Once these files are removed, create a separate directory")
    message(STATUS "and run CMake from there")
    message(FATAL_ERROR "in-source build detected")
endif()

# -------------------------------------------
# CMAKE VERSIONS
# -------------------------------------------
#cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)
cmake_policy(SET CMP0005 NEW)

#---------------------------------------------
# - Project definition
# --------------------------------------------
project(Caesar)

#---Set pathes where to put the libraries, executables and headers------------------------------
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

##set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
##set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/include)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_INSTALL_PREFIX}/bin)
set(BINARY_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/bin)
set(SHARE_SRC_DIR ${CMAKE_SOURCE_DIR}/share)
set(SHARE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/share)
set(CONFIG_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/config)
set(SCRIPT_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/scripts)
set(MACRO_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/macros)
set(DATA_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/data)

message (STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message (STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")
message (STATUS "BINARY_INSTALL_DIR: ${BINARY_INSTALL_DIR}")

set(SFINDER ${CMAKE_BINARY_DIR})
set(HEADER_OUTPUT_PATH ${CMAKE_BINARY_DIR}/include)
set(SFINDER_INCLUDE_DIR ${HEADER_OUTPUT_PATH})


# ----------------------------------------------------------------------------------------------
#---Where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ is checked-------------
# ----------------------------------------------------------------------------------------------
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules)

# - Versioning. We do this here for now
set(${PROJECT_NAME}_VERSION "1.1.0")
set(${PROJECT_NAME}_VERSION_MAJOR "1")
set(${PROJECT_NAME}_VERSION_MINOR "1")
set(${PROJECT_NAME}_VERSION_PATCH "0")

# ----------------------------------------------------------------------------------------------
#---Load some basic macros which are needed later for the configuration and build----------------
#-----------------------------------------------------------------------------------------------
#include(SFinderBuildOptions)
#include(CheckCompiler) # Set compiler flag
#include(MacroEnsureVersion)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
		message(STATUS "Found C++11 support in ${CMAKE_CXX_COMPILER} compiler.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_REENTRANT -Dlinux -g -std=c++11 -O2 -Wall -Wunused -Wuninitialized -fPIC -pthread -m64") 
elseif(COMPILER_SUPPORTS_CXX0X)
		message(STATUS "Found C++0X support in ${CMAKE_CXX_COMPILER} compiler.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_REENTRANT -Dlinux -g -std=c++0x -O2 -Wall -Wunused -Wuninitialized -fPIC -pthread -m64")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()


# - ADD AddressSanitizer flags?
include(sanitizers)
#option(ENABLE_ADDRESS_SANITIZER "Enable AddressSanitizer for deep debug" OFF)
#if (ENABLE_ADDRESS_SANITIZER)
#	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
#	set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
#	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
#	#set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
#	#set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
#endif()

message(STATUS "${CMAKE_CXX_COMPILER} compiler flags: ${CMAKE_CXX_FLAGS}")

#-----------------------------------------------------------------------------------------------
#---Here we look for installed software and switch on and of the different build options--------
include(SearchInstalledSoftware)
#SFINDER_SHOW_OPTIONS()
#-----------------------------------------------------------------------------------------------


#===================================================
#==   INCLUDE OPENMP FLAGS?
#===================================================
##if(BUILD_WITH_OPENMP AND OPENMP_FOUND)
if(ENABLE_OPENMP AND OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#===================================================
#==   INCLUDE MPI FLAGS?
#===================================================
if(ENABLE_MPI AND MPI_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
endif()

#===================================================
#==   INCLUDE OPENCV FLAGS?
#===================================================
add_definitions(-DOPENCV_MAJOR_VERSION=${OpenCV_VERSION_MAJOR})
add_definitions(-DOPENCV_MINOR_VERSION=${OpenCV_VERSION_MINOR})
add_definitions(-DOPENCV_PATCH_VERSION=${OpenCV_VERSION_REVISION})

#==============================================================================================
#===          SET INCLUDES & LIBS
#==============================================================================================
set(INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
set(LINKED_LIBS "")

# Append subdirs (headers & sources files are appended in subdirs)
add_subdirectory(external)
add_subdirectory(proto)
add_subdirectory(src)
if(BUILD_APPS)
	add_subdirectory(apps)
endif()
add_subdirectory(config)
add_subdirectory(scripts)
add_subdirectory(macros)
add_subdirectory(data)


#-- WCSTools
#set(WCSTOOL_DIR ${CMAKE_SOURCE_DIR}/external/wcstools/v3.9.2/libwcs)
#set(WCSTOOL_DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/external/wcstools/v3.9.2)
#set(WCSTOOL_DOWNLOAD_DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/wcstools)
#set(WCSTOOL_BUILD_DIR ${WCSTOOL_DOWNLOAD_DEST_DIR}/v3.9.2)
#set(WCSTOOL_LIB libwcs.a)
#set(WCSTOOL_INC wcs.h)
#set(WCSTOOL_HEADERS ${WCSTOOL_BUILD_DIR}/libwcs/wcs.h)
#list(APPEND WCSTOOL_HEADERS ${WCSTOOL_BUILD_DIR}/libwcs/wcslib.h)
#list(APPEND WCSTOOL_HEADERS ${WCSTOOL_BUILD_DIR}/libwcs/fitshead.h)
#set(WCSTOOL_INCLUDE_DIRS ${CMAKE_INSTALL_INCLUDEDIR})
#set(WCSTOOL_LIBRARIES ${CMAKE_INSTALL_LIBDIR}/${WCSTOOL_LIB})

#-- linterp
set(LINTERP_DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/external/linterp)
set(LINTERP_DOWNLOAD_DEST_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(LINTERP_BUILD_DIR ${LINTERP_DOWNLOAD_DEST_DIR}/linterp)
set(LINTERP_HEADERS ${LINTERP_BUILD_DIR}/linterp.h)
set(LINTERP_INCLUDE_DIRS ${CMAKE_INSTALL_INCLUDEDIR})

# Append external includes found 
list (APPEND INCLUDE_DIRS ${ROOT_INCLUDE_DIRS})
list (APPEND INCLUDE_DIRS ${OpenCV_INCLUDE_DIRS})
list (APPEND INCLUDE_DIRS ${OpenCV2_INCLUDE_DIRS})
if (ENABLE_R)
	list (APPEND INCLUDE_DIRS ${R_INCLUDE_DIR})
endif()
list (APPEND INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
##list (APPEND INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS})
##list (APPEND INCLUDE_DIRS ${WCSTOOL_INCLUDE_DIRS})
list (APPEND INCLUDE_DIRS ${LINTERP_INCLUDE_DIRS})
if (ENABLE_LOGGING)
	list (APPEND INCLUDE_DIRS ${LOG4CXX_INCLUDE_DIRS})
endif()
list (APPEND INCLUDE_DIRS ${JSONCPP_INCLUDE_DIR})
list (APPEND INCLUDE_DIRS ${CFITSIO_INCLUDE_DIR})
list (APPEND INCLUDE_DIRS ${PROTOBUF_INCLUDE_DIR})
##list (APPEND INCLUDE_DIRS ${GSL_INCLUDE_DIRS})
##list (APPEND INCLUDE_DIRS ${WCSLIB_INCLUDE_DIR})

list (APPEND INCLUDE_DIRS ${LINTERP_DOWNLOAD_DIR})
if (ENABLE_TEST)
	list (APPEND INCLUDE_DIRS ${GTEST_INCLUDE_DIRS})
endif()

if (ENABLE_MPI AND MPI_FOUND)
	list (APPEND INCLUDE_DIRS ${MPI_CXX_INCLUDE_PATH})
endif()

if (ENABLE_VTK)
	list (APPEND INCLUDE_DIRS ${VTK_INCLUDE_DIRS})
endif()

if (ENABLE_CASACORE)
	list (APPEND INCLUDE_DIRS ${CASACORE_INCLUDE_DIR})
  list (APPEND INCLUDE_DIRS ${WCSLIB_INCLUDE_DIR})
endif()

include_directories(${INCLUDE_DIRS})

# -- debug --
##message (STATUS "INCLUDE_DIRS: ${INCLUDE_DIRS}")
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
# -----------

# Append external libs found 
list (APPEND LINKED_LIBS ${ROOT_LIBRARIES})
list (APPEND LINKED_LIBS ${OpenCV2_LIBRARIES})
if (ENABLE_R)
	list (APPEND LINKED_LIBS ${R_LIBRARIES})
endif()
list (APPEND LINKED_LIBS ${Boost_LIBRARIES})
##list (APPEND LINKED_LIBS ${PYTHON_LIBRARIES})
##list (APPEND LINKED_LIBS ${WCSTOOL_LIBRARIES})
if (ENABLE_LOGGING)
	list (APPEND LINKED_LIBS ${LOG4CXX_LIBRARIES})
endif()
list (APPEND LINKED_LIBS ${JSONCPP_LIBRARIES})
list (APPEND LINKED_LIBS ${CFITSIO_LIBRARIES})
list (APPEND LINKED_LIBS ${PROTOBUF_LIBRARIES})
##list (APPEND LINKED_LIBS ${GSL_LIBRARIES})
##list (APPEND LINKED_LIBS ${WCSLIB_LIBRARIES})


if (ENABLE_TEST)
  list (APPEND LINKED_LIBS ${GTEST_BOTH_LIBRARIES})
endif()

if (ENABLE_MPI AND MPI_FOUND)
	list (APPEND LINKED_LIBS ${MPI_CXX_LIBRARIES})
endif()

if (ENABLE_VTK)
	list (APPEND LINKED_LIBS ${VTK_LIBRARIES})
endif()

if (ENABLE_CASACORE)
	list (APPEND LINKED_LIBS ${CASACORE_LIBRARIES})
  list (APPEND LINKED_LIBS ${WCSLIB_LIBRARIES})
endif()

message (STATUS "LINKED_LIBS: ${LINKED_LIBS}")


#==============================================================================================
#===          BUILD EXTERNALS
#==============================================================================================
include(ExternalProject)

#ExternalProject_Add(WCSTools
#	DOWNLOAD_COMMAND cp -rp ${WCSTOOL_DOWNLOAD_DIR} ${WCSTOOL_DOWNLOAD_DEST_DIR}
#	BINARY_DIR ${WCSTOOL_BUILD_DIR}
#	BUILD_IN_SOURCE 0
#	CONFIGURE_COMMAND ""
#	BUILD_COMMAND make
#	INSTALL_COMMAND mkdir --parents ${CMAKE_INSTALL_INCLUDEDIR} COMMAND mkdir --parents ${CMAKE_INSTALL_LIBDIR} COMMAND cp ${WCSTOOL_BUILD_DIR}/libwcs/${WCSTOOL_LIB} ${CMAKE_INSTALL_LIBDIR} COMMAND cp ${WCSTOOL_HEADERS} ${CMAKE_INSTALL_INCLUDEDIR}
#)

ExternalProject_Add(LINTERP
	DOWNLOAD_COMMAND cp -rp ${LINTERP_DOWNLOAD_DIR} ${LINTERP_DOWNLOAD_DEST_DIR}
	BINARY_DIR ${LINTERP_BUILD_DIR}
  BUILD_IN_SOURCE 0
  CONFIGURE_COMMAND ""
	BUILD_COMMAND ""
	INSTALL_COMMAND mkdir --parents ${CMAKE_INSTALL_INCLUDEDIR} COMMAND cp ${LINTERP_HEADERS} ${CMAKE_INSTALL_INCLUDEDIR}
)



#==============================================================================================
#===             GENERATE PROTOBUF SOURCES
#==============================================================================================
set(GEN_PROTO_SRCS "")
set(GEN_PROTO_HEADERS "")
foreach (_protoFile ${PROTO_SRCS})
	list (APPEND GEN_PROTO_SRCS ${CMAKE_BINARY_DIR}/${_protoFile})
endforeach()
foreach (_protoFile ${PROTO_HEADERS})
	list (APPEND GEN_PROTO_HEADERS ${CMAKE_BINARY_DIR}/${_protoFile})
endforeach()
list(REMOVE_DUPLICATES GEN_PROTO_SRCS)
list(REMOVE_DUPLICATES GEN_PROTO_HEADERS)
set(PROTO_GEN_FILES "")
list (APPEND PROTO_GEN_FILES ${GEN_PROTO_SRCS})
list (APPEND PROTO_GEN_FILES ${GEN_PROTO_HEADERS})
#set(PROTO_GEN_FILES_COPY ${PROTO_GEN_FILES})

SET(proto_args 
	PROTOFILES ${PROTO_FILES}
	OUTPATH ${CMAKE_BINARY_DIR}
)
#PROTOBUF_GENERATE_CPP(PROTO_GEN_SRCS PROTO_GEN_HEADERS ${proto_args})	

SET_SOURCE_FILES_PROPERTIES(${GEN_PROTO_SRCS} ${GEN_PROTO_HEADERS} PROPERTIES GENERATED TRUE)
	
# Dummy output which is never actually produced. Anything that depends on this will always be rebuilt.
#add_custom_command(
#	OUTPUT always_rebuild
#  COMMAND cmake -E echo
#)

# custom_output will always be rebuilt because it depends on always_rebuild
add_custom_command(
  ####OUTPUT proto_generation
	OUTPUT ${PROTO_GEN_FILES}
	PRE_BUILD
	COMMAND ${PROTOBUF_PROTOC_EXECUTABLE} --proto_path=${PROTO_DIR} --cpp_out=${CMAKE_BINARY_DIR} ${PROTO_FILES}
	#####DEPENDS always_rebuild  ### Decomment this to always generate protobuf files
	DEPENDS ${PROTO_FILES}
	COMMENT "Generating Protobuf source files from schema"
	VERBATIM
)

#add_custom_target(proto_generation
#	COMMAND ${PROTOBUF_PROTOC_EXECUTABLE} --proto_path=${PROTO_DIR} --cpp_out=${CMAKE_BINARY_DIR} ${PROTO_FILES}
#	DEPENDS ${PROTO_FILES}
#)

add_custom_target(ProtobufSchema ALL 
	#DEPENDS proto_generation
	#DEPENDS ${PROTO_FILES}  ## Will generate source files only when proto files are changed
	DEPENDS ${PROTO_GEN_FILES}
)

#-----------------------------------------------------------------------------------------------



#===============================================================================================
#===                    BUILD PROJECT SHARED LIB
#===============================================================================================
#--- Define Linkdef files
set(DICT_OPTIONS "")
list(APPEND DICT_OPTIONS "-DBOOST_LOG_USE_NATIVE_SYSLOG")
list(APPEND DICT_OPTIONS "-DBOOST_LOG_DYN_LINK")
if(ENABLE_CASACORE)
list(APPEND DICT_OPTIONS "-DCASACORE_ENABLED")
endif()

#--- Generate classes dictionary
##MESSAGE (STATUS "DICT_HEADERS: ${DICT_HEADERS}, $DICT_SRCS: ${DICT_SRCS}")
##ROOT_GENERATE_DICTIONARY(CaesarDict ${DICT_HEADERS} MODULE Caesar LINKDEF LinkDef.h)
ROOT_GENERATE_DICTIONARY(CaesarDict ${DICT_HEADERS} 
  LINKDEF LinkDef.h 
#	OPTIONS -DBOOST_LOG_USE_NATIVE_SYSLOG -DBOOST_LOG_DYN_LINK
  OPTIONS ${DICT_OPTIONS}
)


#---Create a shared library with generated dictionary
add_library(Caesar SHARED ${SRCS} CaesarDict.cxx)
##add_dependencies(Caesar ProtobufSchema LINTERP WCSTools)
add_dependencies(Caesar ProtobufSchema LINTERP)
target_link_libraries(Caesar PUBLIC ${LINKED_LIBS})
#===============================================================================================


#===============================================================================================
#===                     BUILD PROJECT TASKS EXEs
#===============================================================================================
if(BUILD_APPS)
	add_executable(FindSourceExe ${Task_FindSource_SRCS})
	target_link_libraries(FindSourceExe Caesar ${LINKED_LIBS})
	set_target_properties(FindSourceExe PROPERTIES OUTPUT_NAME FindSource)

	add_executable(FindSourceResidualExe ${Task_FindSourceResidual_SRCS})
	target_link_libraries(FindSourceResidualExe Caesar ${LINKED_LIBS})
	set_target_properties(FindSourceResidualExe PROPERTIES OUTPUT_NAME FindSourceResidual)

	add_executable(FindSaliencyExe ${Task_FindSaliency_SRCS})
	target_link_libraries(FindSaliencyExe Caesar ${LINKED_LIBS})
	set_target_properties(FindSaliencyExe PROPERTIES OUTPUT_NAME FindSaliency)

	add_executable(ReadImgExe ${Task_ReadImg_SRCS})
	target_link_libraries(ReadImgExe Caesar ${LINKED_LIBS})
	set_target_properties(ReadImgExe PROPERTIES OUTPUT_NAME ReadImg)

	add_executable(FindBkgExe ${Task_FindBkg_SRCS})
	target_link_libraries(FindBkgExe Caesar ${LINKED_LIBS})
	set_target_properties(FindBkgExe PROPERTIES OUTPUT_NAME FindBkg)

	add_executable(FITS2ROOTExe ${Task_FITS2ROOT_SRCS})
	target_link_libraries(FITS2ROOTExe Caesar ${LINKED_LIBS})
	set_target_properties(FITS2ROOTExe PROPERTIES OUTPUT_NAME FITS2ROOT)

	add_executable(ROOT2FITSExe ${Task_ROOT2FITS_SRCS})
	target_link_libraries(ROOT2FITSExe Caesar ${LINKED_LIBS})
	set_target_properties(ROOT2FITSExe PROPERTIES OUTPUT_NAME ROOT2FITS)

	add_executable(CorrelateSourceCatalogsExe ${Task_CorrelateSourceCatalogs_SRCS})
	target_link_libraries(CorrelateSourceCatalogsExe Caesar ${LINKED_LIBS})
	set_target_properties(CorrelateSourceCatalogsExe PROPERTIES OUTPUT_NAME CorrelateSourceCatalogs)
	
	add_executable(SkyModelConvolverExe ${Task_SkyModelConvolver_SRCS})
	target_link_libraries(SkyModelConvolverExe Caesar ${LINKED_LIBS})
	set_target_properties(SkyModelConvolverExe PROPERTIES OUTPUT_NAME SkyModelConvolver)

	add_executable(SourceMatchFinderExe ${Task_SourceMatchFinder_SRCS})
	target_link_libraries(SourceMatchFinderExe Caesar ${LINKED_LIBS})
	set_target_properties(SourceMatchFinderExe PROPERTIES OUTPUT_NAME SourceMatchFinder)

	add_executable(SourceRegionSelectorExe ${Task_SourceRegionSelector_SRCS})
	target_link_libraries(SourceRegionSelectorExe Caesar ${LINKED_LIBS})
	set_target_properties(SourceRegionSelectorExe PROPERTIES OUTPUT_NAME SourceRegionSelector)

	add_executable(SourceRegionMatchFinderExe ${Task_SourceRegionMatchFinder_SRCS})
	target_link_libraries(SourceRegionMatchFinderExe Caesar ${LINKED_LIBS})
	set_target_properties(SourceRegionMatchFinderExe PROPERTIES OUTPUT_NAME SourceRegionMatchFinder)
	
	add_executable(SourceRegionTaggerExe ${Task_SourceRegionTagger_SRCS})
	target_link_libraries(SourceRegionTaggerExe Caesar ${LINKED_LIBS})
	set_target_properties(SourceRegionTaggerExe PROPERTIES OUTPUT_NAME SourceRegionTagger)
	
	add_executable(SourceSelectorExe ${Task_SourceSelector_SRCS})
	target_link_libraries(SourceSelectorExe Caesar ${LINKED_LIBS})
	set_target_properties(SourceSelectorExe PROPERTIES OUTPUT_NAME SourceSelector)

	add_executable(SourceCountsUnfolderExe ${Task_SourceCountsUnfolder_SRCS})
	target_link_libraries(SourceCountsUnfolderExe Caesar ${LINKED_LIBS})
	set_target_properties(SourceCountsUnfolderExe PROPERTIES OUTPUT_NAME SourceCountsUnfolder)

	add_executable(SourceCatalogsMatchFinderExe ${Task_SourceCatalogsMatchFinder_SRCS})
	target_link_libraries(SourceCatalogsMatchFinderExe Caesar ${LINKED_LIBS})
	set_target_properties(SourceCatalogsMatchFinderExe PROPERTIES OUTPUT_NAME SourceCatalogsMatchFinder)

	add_executable(SourceAstroCatalogsMatchFinderExe ${Task_SourceAstroCatalogsMatchFinder_SRCS})
	target_link_libraries(SourceAstroCatalogsMatchFinderExe Caesar ${LINKED_LIBS})
	set_target_properties(SourceAstroCatalogsMatchFinderExe PROPERTIES OUTPUT_NAME SourceAstroCatalogsMatchFinder)

	add_executable(SourceTreeWriterExe ${Task_SourceTreeWriter_SRCS})
	target_link_libraries(SourceTreeWriterExe Caesar ${LINKED_LIBS})
	set_target_properties(SourceTreeWriterExe PROPERTIES OUTPUT_NAME SourceTreeWriter)

	add_executable(SourceAsciiExporterExe ${Task_SourceAsciiExporter_SRCS})
	target_link_libraries(SourceAsciiExporterExe Caesar ${LINKED_LIBS})
	set_target_properties(SourceAsciiExporterExe PROPERTIES OUTPUT_NAME SourceAsciiExporter)

	add_executable(SourceMVClassifierExe ${Task_SourceMVClassifier_SRCS})
	target_link_libraries(SourceMVClassifierExe Caesar ${LINKED_LIBS})
	set_target_properties(SourceMVClassifierExe PROPERTIES OUTPUT_NAME SourceMVClassifier)

	add_executable(SourceRefitterExe ${Task_SourceRefitter_SRCS})
	target_link_libraries(SourceRefitterExe Caesar ${LINKED_LIBS})
	set_target_properties(SourceRefitterExe PROPERTIES OUTPUT_NAME SourceRefitter)

	if(ENABLE_MPI AND MPI_FOUND)
		include_directories(${MPI_CXX_INCLUDE_PATH})
		add_executable(FindSourceMPIExe ${Task_FindSourceMPI_SRCS})
		target_link_libraries(FindSourceMPIExe Caesar ${LINKED_LIBS} ${MPI_CXX_LIBRARIES})
		if(MPI_COMPILE_FLAGS)
  		set_target_properties(FindSourceMPIExe PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
		endif()
		if(MPI_LINK_FLAGS)
  		set_target_properties(FindSourceMPIExe PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
		endif()
		set_target_properties(FindSourceMPIExe PROPERTIES OUTPUT_NAME FindSourceMPI)
	endif()

endif()
#==============================================================================================



#==============================================================================================
#===                    TESTING
#==============================================================================================
if (ENABLE_TEST)
	# Append subdirs (headers & sources files are appended in subdirs)
	set(TEST_INCLUDE_DIRS ${GTEST_INCLUDE_DIRS})
	add_subdirectory(test)

	enable_testing()

  # Include test & Google Test headers
  include_directories(${TEST_INCLUDE_DIRS})

	# Append external libs found 
	list (APPEND LINKED_LIBS ${GTEST_BOTH_LIBRARIES})
	##MESSAGE (STATUS "LINKED_LIBS: ${LINKED_LIBS}")
	
  ##############
  # Unit Tests
  ##############	
	add_executable(runUnitTests ${UNIT_TEST_MAIN} ${UNIT_TEST_SRCS})
	add_dependencies(runUnitTests Caesar)

  # Standard linking to gtest stuff.
	target_link_libraries(runUnitTests PUBLIC ${LINKED_LIBS} Caesar)

  # This is so you can do 'make test' to see all your tests run, instead of
  # manually running the executable runUnitTests to see those specific tests.
	##add_test(UnitTests ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/runUnitTests --gtest_output="xml:UnitTestResults.xml")
	add_test(UnitTests ${BINARY_INSTALL_DIR}/runUnitTests --gtest_output="xml:UnitTestResults.xml")

  # You can also omit NAME and COMMAND. The second argument could be some other
  # test executable.
	# add_test(that-other-test-I-made runUnitTests)
endif()
#==================================================================================================

#==============================================================================================
#===                    DOCUMENTATION
#==============================================================================================
if(BUILD_DOC)
	if(NOT DOXYGEN_FOUND)
  	MESSAGE(FATAL_ERROR "Doxygen is needed to build the documentation.")
  endif()


	# -- Copy Doxygen doc tree in build dir	
	file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/doc)

	# -- Generate doxygen config file (replace cmake variables)
  set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile.in)
  set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen/Doxyfile)
  configure_file(${doxyfile_in} ${doxyfile} @ONLY)

	# -- Generate doxygen doc
	add_custom_target(generate_doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen
    COMMENT "Generating API documentation with Doxygen"
		VERBATIM
	)

  #add_custom_target(doc
  #  COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
	#	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen
  #  COMMENT "Generating API documentation with Doxygen"
	#	VERBATIM
	#)

	

  # -- Install doc files
  #install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc DESTINATION ${SHARE_INSTALL_DIR})
	#install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen/html DESTINATION ${SHARE_INSTALL_DIR}/doc/doxygen)
	#install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen/man DESTINATION ${SHARE_INSTALL_DIR}/doc/doxygen)
	#install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen/rtf DESTINATION ${SHARE_INSTALL_DIR}/doc/doxygen)
	#install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen/latex DESTINATION ${SHARE_INSTALL_DIR}/doc/doxygen)
	#install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen/xml DESTINATION ${SHARE_INSTALL_DIR}/doc/doxygen)

	add_custom_target(doc
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen/html ${SHARE_INSTALL_DIR}/doc/doxygen/html
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen/man ${SHARE_INSTALL_DIR}/doc/doxygen/man
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen/latex ${SHARE_INSTALL_DIR}/doc/doxygen/latex
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen/rtf ${SHARE_INSTALL_DIR}/doc/doxygen/rtf
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen
    COMMENT "Installing API Doxygen documentation in directory ${SHARE_INSTALL_DIR}/doc/doxygen"
		DEPENDS generate_doc
		VERBATIM
	)

endif()


if(BUILD_WIKI_DOC)
	if(NOT DOXYGEN_FOUND)
  	MESSAGE(FATAL_ERROR "Doxygen is needed to build the Sphinx documentation via breathe.")
  endif()
	if(NOT SPHINX_FOUND)
  	MESSAGE(FATAL_ERROR "Sphinx is needed to build the documentation.")
  endif()

	# -- Copy Sphinx doc tree in build dir	
	file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/doc/sphinx DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/doc)

	# -- Copy LICENSE file in build dir
	file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/COPYRIGHT.md DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

	# -- Generate doxygen config files
  set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/doc/sphinx/Doxyfile.in)
  set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/doc/sphinx/Doxyfile_gen)
  configure_file(${doxyfile_in} ${doxyfile} @ONLY)

  ##add_custom_target(doxydoc
	add_custom_target(generate_doxygen_doc
		COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/sphinx
    COMMENT "Generating XML Doxygen documentation to be used by sphinx breathe"
		VERBATIM
	)

	# -- Generate Sphinx html documentation
	add_custom_target(generate_sphinx_doc
		COMMAND sphinx-build -b html . html
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/sphinx
   	DEPENDS generate_doxygen_doc
   	COMMENT "Build Sphinx html documentation"
		VERBATIM
	)	

	# -- Configure sphinx builder
	#macro(add_sphinx_target TARGET_NAME BUILDER COMMENT_STR)
  #  add_custom_target(${TARGET_NAME}
  #  	COMMAND sphinx-build -b ${BUILDER} . sphinx/${BUILDER}
	#		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/sphinx
  # 		DEPENDS doxydoc
  # 		COMMENT ${COMMENT_STR}
  #	)
	#endmacro(add_sphinx_target)

	##add_sphinx_target(doc-sphinx html "Build Sphinx html documentation")
	
  # -- Install doc files
  #install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc DESTINATION ${SHARE_INSTALL_DIR})
	#install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/sphinx/sphinx/html DESTINATION ${SHARE_INSTALL_DIR}/doc/sphinx)
	
	add_custom_target(doc-sphinx
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/doc/sphinx/html ${SHARE_INSTALL_DIR}/doc/sphinx/html
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Installing sphinx generated html doc tree in ${SHARE_INSTALL_DIR}/doc/sphinx dir"
		DEPENDS generate_sphinx_doc
		VERBATIM
	)
	

endif()
	

#==============================================================================================
#===              INSTALL STUFF
#==============================================================================================
#-- Install targets
set(EXE_FILES "")
if(BUILD_APPS)
	list(APPEND EXE_FILES 
		FindSource
		FindSourceResidual
		FindSaliency
		FindBkg	
		ReadImg
		FITS2ROOT
		ROOT2FITS
		CorrelateSourceCatalogs
		SkyModelConvolver
		SourceMatchFinder
		SourceRegionSelector
		SourceRegionMatchFinder
		SourceRegionTagger
		SourceSelector
		SourceCountsUnfolder
		SourceCatalogsMatchFinder
		SourceAstroCatalogsMatchFinder
		SourceTreeWriter
		SourceAsciiExporter
		SourceMVClassifier
		SourceRefitter
	)
	if(ENABLE_MPI)
		list(APPEND EXE_FILES FindSourceMPI)
	endif()

	foreach (_exeFile ${EXE_FILES})
		install(TARGETS ${_exeFile}Exe DESTINATION ${BINARY_INSTALL_DIR})
	endforeach()
	
endif()

# -- Install library dict map
# NB: EXISTS macro not working because these files are generated only after make, try use the OPTIONAL
#if(EXISTS "${PROJECT_BINARY_DIR}/CaesarDict_rdict.pcm")
#	install(FILES ${PROJECT_BINARY_DIR}/CaesarDict_rdict.pcm DESTINATION ${CMAKE_INSTALL_LIBDIR})
#endif()
#if(EXISTS "${PROJECT_BINARY_DIR}/libCaesarDict_rdict.pcm")
#	install(FILES ${PROJECT_BINARY_DIR}/libCaesarDict_rdict.pcm DESTINATION ${CMAKE_INSTALL_LIBDIR})
#endif()
install(FILES ${PROJECT_BINARY_DIR}/CaesarDict_rdict.pcm DESTINATION ${CMAKE_INSTALL_LIBDIR} OPTIONAL)
install(FILES ${PROJECT_BINARY_DIR}/libCaesarDict_rdict.pcm DESTINATION ${CMAKE_INSTALL_LIBDIR} OPTIONAL)


# -- Install libraries
install(TARGETS Caesar DESTINATION ${CMAKE_INSTALL_LIBDIR})

#-- Install project headers
install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

#-- Install test executable
if (ENABLE_TEST)
	install(TARGETS runUnitTests DESTINATION ${BINARY_INSTALL_DIR})
endif()

#-- Install external headers
## WCS
#install(FILES ${WCSTOOL_BUILD_DIR}/libwcs/${WCSTOOL_LIB} DESTINATION ${CMAKE_INSTALL_LIBDIR})
#install(FILES ${WCSTOOL_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

## LINTERP
install(FILES ${LINTERP_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# -- Install share files
file(GLOB_RECURSE SHARE_FILES ${SHARE_SRC_DIR}/*)
install(FILES ${SHARE_FILES} DESTINATION ${SHARE_INSTALL_DIR})

# -- Install config files
install(FILES ${CONFIG_FILES} DESTINATION ${CONFIG_INSTALL_DIR})

# -- Install script files
install(FILES ${SCRIPT_FILES} 
	PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
	DESTINATION ${SCRIPT_INSTALL_DIR}
)

# -- Install macro files
install(FILES ${MACRO_FILES} 
	PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
	DESTINATION ${MACRO_INSTALL_DIR}
)

# -- Install data files
##MESSAGE(STATUS DATA_FILES: ${DATA_FILES})
install(FILES ${DATA_FILES} DESTINATION ${DATA_INSTALL_DIR})

# -- Install documentation files
if(BUILD_DOC)
	install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc DESTINATION ${SHARE_INSTALL_DIR})
endif()

#==================================================================================================

#===============================================================================================
#==               CLEAN STUFF
#===============================================================================================
set(CLEANED_HEADERS "")
set(CLEANED_LIBS "")
set(CLEANED_EXE "")

#--- Set external WCTOOLS lib & headers to be removed by clean command
list(APPEND CLEANED_LIBS ${WCSTOOL_LIBRARIES})
foreach (_headerFile ${WCSTOOL_HEADERS})
	get_filename_component(_basefilename ${_headerFile} NAME)
	set(_cleanedHeader ${CMAKE_INSTALL_INCLUDEDIR}/${_basefilename})
	list(APPEND CLEANED_HEADERS ${_cleanedHeader})
endforeach()

#-- Set linterp headers to be removed by clean command
foreach (_headerFile ${LINTERP_HEADERS})
	get_filename_component(_basefilename ${_headerFile} NAME)
	set(_cleanedHeader ${CMAKE_INSTALL_INCLUDEDIR}/${_basefilename})
	list(APPEND CLEANED_HEADERS ${_cleanedHeader})
endforeach()

#--- Set CAESAR headers & libs & exe to be removed by clean command
list(APPEND CLEANED_LIBS ${CMAKE_INSTALL_LIBDIR}/CaesarDict_rdict.pcm)
list(APPEND CLEANED_LIBS ${CMAKE_INSTALL_LIBDIR}/libCaesar.so)

foreach (_headerFile ${HEADERS_WPATH})
	set(_cleanedHeader ${CMAKE_INSTALL_INCLUDEDIR}/${_headerFile})
	list(APPEND CLEANED_HEADERS ${_cleanedHeader})
endforeach()

if(BUILD_APPS)
	foreach (_exeFile ${EXE_FILES})
		set(_cleanedExe ${BINARY_INSTALL_DIR}/${_exeFile})
		list(APPEND CLEANED_EXE ${_cleanedExe})
	endforeach()
endif()

if (ENABLE_TEST)	
	list(APPEND CLEANED_EXE ${BINARY_INSTALL_DIR}/runUnitTests)
endif()

# Remove duplicated and add to additional clean targets
list(REMOVE_DUPLICATES CLEANED_HEADERS)
list(REMOVE_DUPLICATES CLEANED_LIBS)
list(REMOVE_DUPLICATES CLEANED_EXE)
list(APPEND ADDITIONAL_CLEAN_TARGETS ${CLEANED_HEADERS})
list(APPEND ADDITIONAL_CLEAN_TARGETS ${CLEANED_LIBS})
list(APPEND ADDITIONAL_CLEAN_TARGETS ${CLEANED_EXE})
##MESSAGE(STATUS "ADDITIONAL_CLEAN_TARGETS: ${ADDITIONAL_CLEAN_TARGETS}")
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${ADDITIONAL_CLEAN_TARGETS}") ## do not forget to put "" in final vars otherwise you'll get a parsing error
#=================================================================================================


#==============================================================================================
#==           PACKAGING 
#==============================================================================================
# all install commands get the same destination. this allows us to use paths
# relative to the executable.
#install(TARGETS example DESTINATION example_destination)
# this is basically a repeat of the file copy instruction that copies the
# resources in the build directory, but here we tell cmake that we want it
# in the package
#install(DIRECTORY resources DESTINATION example_destination)

# now comse everything we need, to create a package
# there are a lot more variables you can set, and some
# you need to set for some package types, but we want to
# be minimal here
set(CPACK_PACKAGE_NAME "CAESAR")
set(CPACK_PACKAGE_VERSION "1.0.0")

# we don't want to split our program up into several things
set(CPACK_MONOLITHIC_INSTALL 1)

# This must be last
include(CPack)
