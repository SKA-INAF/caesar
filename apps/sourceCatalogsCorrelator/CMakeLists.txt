set(INC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

## GET SOURCES and HEADERS
file(GLOB THIS_SRCS ${SRC_DIR}/*.cc)
file(GLOB THIS_HEADERS ${INC_DIR}/*.h)

## Append srcs and heads to project srcs and heads
list (APPEND Task_CorrelateSourceCatalogs_SRCS ${THIS_SRCS})
list (APPEND Task_CorrelateSourceCatalogs_SRCS ${THIS_HEADERS})

set (Task_CorrelateSourceCatalogs_SRCS ${Task_CorrelateSourceCatalogs_SRCS} PARENT_SCOPE)

list (APPEND INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
set (INCLUDE_DIRS ${INCLUDE_DIRS} PARENT_SCOPE)
