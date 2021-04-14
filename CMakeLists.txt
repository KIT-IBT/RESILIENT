
cmake_minimum_required(VERSION 3.11)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

#Project name
project("RESILIENT")

# Move binaries to vtkToolbox/RESILIENT
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/RESILIENT)


# Try to find the libraries.
# VTK (https://gitlab.kitware.com/vtk/vtk),
# Cork (https://github.com/gilbo/cork), 
# Tetegen (https://github.com/pyvista/tetgen),
# and include path.
#
# Once done this will define:
# VTK_FOUND
# VTK_DIR
#
# CORK_FOUND
# CORK_DIR
# CORK_LIBRARY
#
# TETGEN_FOUND
# TETGEN_DIR

########################
## Set up VTK dependency
########################

# Select only VTK components (libraries to be linked to) that are needed.
# (Removing the COMPONENTS argument from find_package links all of them.)
# See folder structure in https://vtk.org/doc/nightly/html/files.html
# to find out to which component a VTK class belongs.

set(VTK_COMPONENTS
    vtkCommonCore
    vtkCommonDataModel
    vtkCommonExecutionModel
    vtkCommonMath
    #vtkCommonMisc
    #vtkCommonSystem
    vtkFiltersCore
    #vtkFiltersExtraction
    #vtkFiltersFlowPaths
    #vtkFiltersGeneral
    #vtkFiltersGeometry
    #vtkFiltersModeling
    #vtkFiltersParallel
    #vtkFiltersVerdict
    vtkIOCore
    vtkIOGeometry
    vtkIOLegacy
    #vtkIOPLY
    vtkIOXML
)

find_package(VTK REQUIRED COMPONENTS ${VTK_COMPONENTS} NO_MODULE)

    if(NOT VTK_FOUND)
        message(FATAL_ERROR "VTK could not be found. Try setting VTK_DIR.")
    else()
        message(STATUS "Found VTK: " ${VTK_DIR})
        include(${VTK_USE_FILE})
    endif()

    #if(VTK_VERSION LESS 8.9)
    #    include(${VTK_USE_FILE})
    #endif()

########################
## Set up TETGEN dependency
########################


FIND_PATH(TETGEN_DIR tetgen.h)

    if(NOT TETGEN_DIR)
         message(FATAL_ERROR "TETGEN could not be found. Try setting TETGEN_DIR.")
    else()
        message(STATUS "Found TETGEN: " ${TETGEN_DIR})
        include_directories(${TETGEN_DIR})
    endif()

########################
## Set up CORK dependency
########################


FIND_PATH(CORK_INCLUDE_DIR cork.h)

    if(NOT CORK_FOUND)
         message(FATAL_ERROR "CORK could not be found. Try setting CORK_DIR.")
    else()
        message(STATUS "Found CORK: " ${CORK_DIR})
        include_directories(${CORK_DIR}"/include")
    endif()

find_library(CORK_LIBRARY NAMES cork)

    if (NOT CORK_LIBRARY)
	   message(STATUS "Found CORK library: ${CORK_LIBRARY}")
       include_directories(${CORK_LIBRARY})
    else()
       message(FATAL_ERROR "CORK liberay could not be found. Try setting CORK_LIBRARY.")
    endif ()


# Tell compiler where to look for libraries and header files
include_directories( 
    "."
    "./src"
)

add_executable(
    "RESILIENT" "SetAtrialFiberOrientation.cpp" "FiberOrientation.cpp" "AveragingOrientation.cpp" "Config.cpp" "ConvertTetgenio.cpp" "DataFormat.cpp" "Methods.cpp" "Reader.cpp" "Writer.cpp"
    "TestTetrahedralize" "TestTetrahedralize.cpp" "AveragingOrientation.cpp" "Config.cpp" "ConvertTetgenio.cpp" "DataFormat.cpp" "Methods.cpp" "Reader.cpp" "Writer.cpp"
    "FindAndMarkSeedPoints" "FindAndMarkSeedPoints.cpp" "AveragingOrientation.cpp" "Config.cpp" "ConvertTetgenio.cpp" "DataFormat.cpp" "Methods.cpp" "Reader.cpp" "Writer.cpp"
    )
install (TARGETS RESILIENT DESTINATION ${PROJECT_SOURCE_DIR}/bin/)

