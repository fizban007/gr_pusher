# - Try to find the CUDPP libraries
# Once done this will define
#
# CUDPP_FOUND - system has CUDPP
# CUDPP_INCLUDE_DIR - the CUDPP include directory
# CUDPP_LIBRARIES - CUDPP library

MESSAGE(STATUS "Looking for CUDPP...")

FIND_PATH(CUDPP_INCLUDE_DIR cudpp.h PATHS 
	/usr/include/ 
	/usr/local/include/
	${CUDA_INCLUDE_DIRS}
	)

# IF(CMAKE_SIZEOF_VOID_P EQUAL 8)
# 	FIND_LIBRARY(CUDPP_LIBRARIES NAMES cudpp64 PATHS 
# 	    /usr/local/cuda/lib64/
#             /usr/local/lib/
# 	    )
# ELSE()
# 	FIND_LIBRARY(CUDPP_LIBRARIES NAMES cudpp PATHS 
# 	    /usr/local/cuda/lib/ 
#             /usr/local/lib/
# 	    )
# ENDIF()
FIND_LIBRARY(CUDPP_LIBRARIES NAMES libcudpp.so cudpp PATHS 
  /usr/local/cuda/lib/ 
  /usr/local/lib/
  )

IF(CUDPP_INCLUDE_DIR AND CUDPP_LIBRARIES)
	SET(CUDPP_FOUND 1)
	IF(NOT CUDPP_FIND_QUIETLY)
		MESSAGE(STATUS "Found CUDPP: libraries = ${CUDPP_LIBRARIES}")
	ENDIF(NOT CUDPP_FIND_QUIETLY)
ELSE()
	SET(CUDPP_FOUND 0 CACHE BOOL "CUDPP not found")
	IF(CUDPP_FIND_REQUIRED)
	    MESSAGE(FATAL_ERROR "Could NOT find CUDPP, error")
	ELSE()
	    MESSAGE(STATUS "Could NOT find CUDPP, disabled")
	ENDIF()
ENDIF()

MARK_AS_ADVANCED(CUDPP_INCLUDE_DIR CUDPP_LIBRARIES)

