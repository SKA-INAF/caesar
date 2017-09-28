## INC & SOURCE DIRS
set(INC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

## GET SOURCES and HEADERS
file(GLOB_RECURSE THIS_SRCS ${SRC_DIR}/*.cc)
file(GLOB_RECURSE THIS_HEADERS ${INC_DIR}/*.h)
set(THIS_DICT_HEADERS ${THIS_DICT_HEADERS} 
	SysUtils.h 
	AstroUtils.h 
	StatsUtils.h 
	CodeUtils.h 
	MathUtils.h 
	GraphicsUtils.h 
	Option.h 
	ConfigParser.h 
	Logger.h
	TaskData.h
	Serializer.h
	Graph.h
)
set(THIS_DICT_SRCS ${THIS_DICT_SRCS} 
	${SRC_DIR}/SysUtils.cc 
	${SRC_DIR}/AstroUtils.cc 
	${SRC_DIR}/StatsUtils.cc 
	${SRC_DIR}/CodeUtils.cc 
	${SRC_DIR}/MathUtils.cc 
	${SRC_DIR}/GraphicsUtils.cc 
	${SRC_DIR}/Option.cc 	
	${SRC_DIR}/ConfigParser.cc 
	${SRC_DIR}/Logger.cc
	${SRC_DIR}/TaskData.cc
	${SRC_DIR}/Serializer.cc
	${SRC_DIR}/Graph.cc
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



