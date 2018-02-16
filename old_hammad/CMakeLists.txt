#=============================================================================
# CMake configuration file for SLS Roller
#=============================================================================

cmake_minimum_required(VERSION 2.8)
project(SLS_Roller)
include(cmake/SBELUtils.cmake)
# Policy CMP0046 warning when a dependency does not exist
# According to a dev from kitware The add_dependencies
# command is for inter-target ordering dependencies, not
# file-level dependencies.
# This warning appears in Cmake 3.0
IF(CMAKE_MAJOR_VERSION EQUAL 3)
    IF(COMMAND cmake_policy)
        cmake_policy(SET CMP0046 OLD)
    ENDIF(COMMAND cmake_policy)
ENDIF(CMAKE_MAJOR_VERSION EQUAL 3)
IF(APPLE)
    set(CMAKE_MACOSX_RPATH OFF)
ENDIF()

# Location of additional CMAKE modules
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# ------------------------------------------------------------------------------
# Specify output locations for libraries and executables
# Set the build type depending on the platform
# ------------------------------------------------------------------------------
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# ------------------------------------------------------------------------------
# Find the ChronoEngine library
# ------------------------------------------------------------------------------
add_definitions(-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP)

FIND_PACKAGE(GLM)
FIND_PACKAGE(GLFW)
#FIND_PACKAGE(ZLIB)

find_package(Chrono
             COMPONENTS PARALLEL VEHICLE OPENGL FEA
             CONFIG
)

# ------------------------------------------------------------------------------
# Find CUDA
# ------------------------------------------------------------------------------
enable_cuda_support()
# ------------------------------------------------------------------------------
# Enable OpenMP
# ------------------------------------------------------------------------------
#enable_openmp_support()


#--------------------------------------------------------------
# Generate the configuration file 

#--------------------------------------------------------------
# Set the include paths and libraries

include_directories(
    ${GLFW_INCLUDE_DIR}
    ${CHRONO_INCLUDE_DIRS}
    "${CUDA_TOOLKIT_ROOT_DIR}/include"
    ${ZLIB_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}
    )


SET(LIBRARIES  ${CHRONO_LIBRARIES} ${ZLIB_LIBRARIES} ${GLFW_LIBRARY})


#--------------------------------------------------------------

IF(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${CHRONO_CXX_FLAGS})
ENDIF()

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    SET (CH_LINKERFLAG_EXE  "${CH_LINKERFLAG_EXE} -framework IOKit -framework Cocoa -framework OpenGL")
ENDIF()

#--------------------------------------------------------------
# List of all executables

SET(DEMOS
    sls
    )

MESSAGE(STATUS "Vehicle fording programs...")

FOREACH(PROGRAM ${DEMOS})
    MESSAGE(STATUS "...add ${PROGRAM}")

    ADD_EXECUTABLE(${PROGRAM}  "${PROGRAM}.cpp")
    SOURCE_GROUP(""  FILES "${PROGRAM}.cpp")

    SET_TARGET_PROPERTIES(${PROGRAM} PROPERTIES
        FOLDER demos
        COMPILE_FLAGS "${CH_BUILDFLAGS}"
        LINK_FLAGS "${CH_LINKERFLAG_EXE}"
        )

    TARGET_LINK_LIBRARIES(${PROGRAM} ${LIBRARIES})
    ADD_DEPENDENCIES(${PROGRAM} ${LIBRARIES})

    INSTALL(TARGETS ${PROGRAM} DESTINATION bin)

ENDFOREACH(PROGRAM)


