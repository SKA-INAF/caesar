## INC & SOURCE DIRS
set(INC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

## GET SOURCES and HEADERS
##file(GLOB THIS_SRCS ${SRC_DIR}/*.cc)
set(THIS_SRCS ${THIS_SRCS} 
	${SRC_DIR}/SysUtils.cc 
	${SRC_DIR}/AstroUtils.cc 
	${SRC_DIR}/StatsUtils.cc 
	${SRC_DIR}/CodeUtils.cc 
	${SRC_DIR}/MathUtils.cc 
	${SRC_DIR}/GraphicsUtils.cc 
	${SRC_DIR}/EllipseUtils.cc
	${SRC_DIR}/Graph.cc 	 
	${SRC_DIR}/Option.cc 	
	${SRC_DIR}/Cut.cc 	
	${SRC_DIR}/ConfigParser.cc 
	${SRC_DIR}/CutParser.cc 
#	${SRC_DIR}/Logger.cc
	${SRC_DIR}/TaskData.cc
	${SRC_DIR}/Serializer.cc
	${SRC_DIR}/ImgUtils.cc
	${SRC_DIR}/WCSUtils.cc
	${SRC_DIR}/DS9Region.cc
	${SRC_DIR}/DS9RegionParser.cc
	${SRC_DIR}/AstroObject.cc
	${SRC_DIR}/AstroObjectParser.cc
) 


##file(GLOB THIS_HEADERS ${INC_DIR}/*.h)
set(THIS_HEADERS ${THIS_HEADERS} 
	${INC_DIR}/SysUtils.h 
	${INC_DIR}/AstroUtils.h 
	${INC_DIR}/StatsUtils.h 
	${INC_DIR}/CodeUtils.h 
	${INC_DIR}/MathUtils.h 
	${INC_DIR}/GraphicsUtils.h 
	${INC_DIR}/EllipseUtils.h
	${INC_DIR}/Graph.h 
	${INC_DIR}/Option.h 
	${INC_DIR}/Cut.h 	
	${INC_DIR}/ConfigParser.h
	${INC_DIR}/CutParser.h  
#	${INC_DIR}/Logger.h
	${INC_DIR}/TaskData.h
	${INC_DIR}/Serializer.h
	${INC_DIR}/ImgUtils.h
	${INC_DIR}/WCSUtils.h
	${INC_DIR}/DS9Region.h
	${INC_DIR}/DS9RegionParser.h
	${INC_DIR}/AstroObject.h
	${INC_DIR}/AstroObjectParser.h
)


set(THIS_DICT_HEADERS ${THIS_DICT_HEADERS} 
	SysUtils.h 
	AstroUtils.h 
	StatsUtils.h 
	CodeUtils.h 
	MathUtils.h 
	GraphicsUtils.h 
	EllipseUtils.h
	Graph.h
	Option.h
	Cut.h 
	ConfigParser.h 
	CutParser.h 
#	Logger.h
	TaskData.h
	Serializer.h
	ImgUtils.h
	WCSUtils.h
	DS9Region.h
	DS9RegionParser.h
	AstroObject.h
	AstroObjectParser.h
)

set(THIS_DICT_SRCS ${THIS_DICT_SRCS} 
	${SRC_DIR}/SysUtils.cc 
	${SRC_DIR}/AstroUtils.cc 
	${SRC_DIR}/StatsUtils.cc 
	${SRC_DIR}/CodeUtils.cc 
	${SRC_DIR}/MathUtils.cc 
	${SRC_DIR}/GraphicsUtils.cc 
	${SRC_DIR}/EllipseUtils.cc 
	${SRC_DIR}/Graph.cc 
	${SRC_DIR}/Option.cc 	
	${SRC_DIR}/Cut.cc 	
	${SRC_DIR}/ConfigParser.cc 
	${SRC_DIR}/CutParser.cc 
#	${SRC_DIR}/Logger.cc
	${SRC_DIR}/TaskData.cc
	${SRC_DIR}/Serializer.cc
	${SRC_DIR}/ImgUtils.cc
	${SRC_DIR}/WCSUtils.cc
	${SRC_DIR}/DS9Region.cc
	${SRC_DIR}/DS9RegionParser.cc
	${SRC_DIR}/AstroObject.cc
	${SRC_DIR}/AstroObjectParser.cc
) 

if (ENABLE_LOGGING)
	list (APPEND THIS_SRCS ${SRC_DIR}/Logger.cc)
	list (APPEND THIS_HEADERS ${INC_DIR}/Logger.h)
	list (APPEND THIS_DICT_HEADERS Logger.h)
	list (APPEND THIS_DICT_SRCS ${SRC_DIR}/Logger.cc)
endif()


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



