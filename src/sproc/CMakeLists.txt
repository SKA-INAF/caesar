## INC & SOURCE DIRS
set(INC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

## GET SOURCES and HEADERS
file(GLOB THIS_SRCS ${SRC_DIR}/*.cc)
file(GLOB THIS_HEADERS ${INC_DIR}/*.h)
set(
	THIS_DICT_HEADERS 
	${THIS_DICT_HEADERS} 
	Blob.h
	Source.h
	SourceCube.h
	SourceGroup.h
	SpectralIndexData.h
	SourceMatchData.h
	SourceComponentPars.h
	SourceFitPars.h
	SourceFitter.h
	SourceExporter.h
	SourceSelector.h
	SourceImporter.h
)
set(
	THIS_DICT_SRCS 
	${THIS_DICT_SRCS} 
	${SRC_DIR}/Blob.cc
	${SRC_DIR}/Source.cc
	${SRC_DIR}/SourceCube.cc	
	${SRC_DIR}/SourceGroup.cc
	${SRC_DIR}/SpectralIndexData.cc
	${SRC_DIR}/SourceMatchData.cc	
	${SRC_DIR}/SourceComponentPars.cc
	${SRC_DIR}/SourceFitPars.cc
	${SRC_DIR}/SourceFitter.cc
	${SRC_DIR}/SourceExporter.cc
	${SRC_DIR}/SourceSelector.cc
	${SRC_DIR}/SourceImporter.cc
)

set(THIS_HEADERS_WPATH "")
foreach (_headerFile ${THIS_HEADERS})
	get_filename_component(_basefilename ${_headerFile} NAME)
  list (APPEND THIS_HEADERS_WPATH ${_basefilename})
endforeach()
list(REMOVE_DUPLICATES THIS_HEADERS_WPATH)


## Append srcs and heads to project srcs and heads
list (APPEND SRCS ${THIS_SRCS})
list (APPEND HEADERS ${THIS_HEADERS})
list (APPEND HEADERS_WPATH ${THIS_HEADERS_WPATH})
list (APPEND DICT_HEADERS ${THIS_DICT_HEADERS})
list (APPEND DICT_SRCS ${THIS_DICT_SRCS})

set (SRCS ${SRCS} PARENT_SCOPE)
set (HEADERS ${HEADERS} PARENT_SCOPE)
set (HEADERS_WPATH ${HEADERS_WPATH} PARENT_SCOPE)
set (DICT_HEADERS ${DICT_HEADERS} PARENT_SCOPE)
set (DICT_SRCS ${DICT_SRCS} PARENT_SCOPE)

list (APPEND INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
set (INCLUDE_DIRS ${INCLUDE_DIRS} PARENT_SCOPE)
