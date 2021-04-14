
cmake_minimum_required(VERSION 3.11)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

project("RESILIENT")

find_package(VTK REQUIRED)
    if(NOT VTK_FOUND)
        message(FATAL_ERROR "VTK could not be found. Try setting VTK_DIR.")
    else()
        message(STATUS "Found VTK: " ${VTK_DIR})
    endif()
    if(VTK_VERSION LESS 8.9)
        include(${VTK_USE_FILE})
    endif()

FIND_PATH(TETGEN_DIR tetgen.h HINTS tetgen1.6.0)
FIND_LIBRARY(TETGEN_LIBRARY NAMES tet HINTS ${TETGEN_DIR})
    if(NOT TETGEN_DIR)
         message(FATAL_ERROR "TETGEN could not be found. Try setting TETGEN_DIR.")
    else()
        message(STATUS "Found TETGEN: " ${TETGEN_DIR})
        include_directories(${TETGEN_DIR})
    endif()

FIND_PATH(CORK_INCLUDE_DIR cork.h HINTS cork/include ${CORK_DIR}/include)
FIND_LIBRARY(CORK_LIBRARY NAMES cork HINTS cork/lib ${CORK_DIR}/lib)
SET(CORK_FOUND "NO")
IF (CORK_INCLUDE_DIR AND CORK_LIBRARY)
    SET(CORK_FOUND "YES")
ENDIF (CORK_INCLUDE_DIR AND CORK_LIBRARY)

if(CORK_FOUND)
  message(STATUS "Found CORK: ${CORK_INCLUDE_DIR}")
  include_directories(${CORK_INCLUDE_DIR})
else(CORK_FOUND)
  message(FATAL_ERROR "Could NOT find CORK. Try setting CORK_DIR.")  
endif(CORK_FOUND)

find_library(GMP_LIBRARY gmp)

add_executable(RESILIENT
    src/SetAtrialFiberOrientation.cpp
    src/FiberOrientation.cpp
    src/AveragingOrientation.cpp
    src/Config.cpp
    src/ConvertTetgenio.cpp
    src/DataFormat.cpp
    src/Methods.cpp
    src/Reader.cpp
    src/Writer.cpp
    src/AveragingOrientation.cpp
    src/Config.cpp
    src/ConvertTetgenio.cpp
    src/DataFormat.cpp
    src/Methods.cpp
    src/Reader.cpp
    src/Writer.cpp
    src/AveragingOrientation.cpp
    src/Config.cpp
    src/ConvertTetgenio.cpp
    src/DataFormat.cpp
    src/Methods.cpp
    src/Reader.cpp
    src/Writer.cpp
)

target_link_libraries(RESILIENT
    ${VTK_LIBRARIES}
    ${TETGEN_LIBRARY}    
    ${CORK_LIBRARY}
    ${GMP_LIBRARY}
)

add_executable(TestTetrahedralize
    src/TestTetrahedralize.cpp    
    src/DataFormat.cpp
    src/Reader.cpp
    src/Writer.cpp
    src/Methods.cpp
    src/Config.cpp
    src/ConvertTetgenio.cpp
    src/AveragingOrientation.cpp    
)

target_link_libraries(TestTetrahedralize
    ${VTK_LIBRARIES}
    ${TETGEN_LIBRARY}    
    ${CORK_LIBRARY}
    ${GMP_LIBRARY}
)

add_executable(FindAndMarkSeedPoints
    src/FindAndMarkSeedPoints.cpp    
    src/DataFormat.cpp
    src/Reader.cpp
    src/Writer.cpp
    src/Config.cpp
    src/Methods.cpp
    src/ConvertTetgenio.cpp
    src/AveragingOrientation.cpp    
)

target_link_libraries(FindAndMarkSeedPoints
    ${VTK_LIBRARIES}
    ${TETGEN_LIBRARY}    
    ${CORK_LIBRARY}
    ${GMP_LIBRARY}
)


install (TARGETS RESILIENT TestTetrahedralize FindAndMarkSeedPoints DESTINATION bin/)
