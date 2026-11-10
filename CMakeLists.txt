cmake_minimum_required(VERSION 3.13...3.31)

project(hemoFlow LANGUAGES C CXX)

# Build with:  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
# Debug with:  gdb --args ./build/hemoFlowd /path/to/config.xml

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Palabos (upstream master after 2.3.0) is written for C++20; its accelerated
# lattice headers, included by palabos3D.h, use C++17 parallel algorithms.
if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 20)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_DEBUG_POSTFIX d)

option(ENABLE_MPI "Build with MPI" ON)
option(BUILD_HDF5 "Build the (parallel) HDF5 writer if parallel HDF5 is found" ON)

include(CheckCXXCompilerFlag)

###############################################################################
# Dependencies
###############################################################################

if(ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
endif()

# cnpy needs zlib to read the compressed .npz geometry files.
find_package(ZLIB REQUIRED)

set(HEMOFLOW_WITH_HDF5 OFF)
if(BUILD_HDF5)
    if(NOT ENABLE_MPI)
        message(WARNING "The HDF5 writer needs MPI; building without HDF5 output.")
    else()
        set(HDF5_PREFER_PARALLEL ON)
        find_package(HDF5 COMPONENTS C)
        if(HDF5_FOUND AND HDF5_IS_PARALLEL)
            message(STATUS "Parallel HDF5 found (${HDF5_VERSION}): HDF5 writer enabled")
            set(HEMOFLOW_WITH_HDF5 ON)
        elseif(HDF5_FOUND)
            message(WARNING "Only a serial HDF5 was found; building without HDF5 output.")
        else()
            message(WARNING "HDF5 was not found; building without HDF5 output.")
        endif()
    endif()
endif()

###############################################################################
# Common compiler settings
###############################################################################

set(HEMOFLOW_WARNING_FLAGS "")
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    list(APPEND HEMOFLOW_WARNING_FLAGS -Wall -Wnon-virtual-dtor)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        list(APPEND HEMOFLOW_WARNING_FLAGS -Wno-unused-result)
    endif()
endif()

# GCC >= 15 rejects ill-formed but never-instantiated templates in the vendored
# Palabos (multiGrid/). Downgrade that error to a warning where supported.
check_cxx_compiler_flag(-Wno-error=template-body HEMOFLOW_HAS_WNO_TEMPLATE_BODY)
if(HEMOFLOW_HAS_WNO_TEMPLATE_BODY)
    list(APPEND HEMOFLOW_WARNING_FLAGS -Wno-error=template-body)
endif()

###############################################################################
# Palabos library (external dependency, vendored in palabos/)
###############################################################################

file(GLOB_RECURSE PALABOS_SRC CONFIGURE_DEPENDS "palabos/src/*.cpp")
file(GLOB_RECURSE PALABOS_EXT_SRC CONFIGURE_DEPENDS "palabos/externalLibraries/tinyxml/*.cpp")

add_library(palabos STATIC ${PALABOS_SRC} ${PALABOS_EXT_SRC})
target_include_directories(palabos SYSTEM PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/palabos/src"
    "${CMAKE_CURRENT_SOURCE_DIR}/palabos/externalLibraries"
    "${CMAKE_CURRENT_SOURCE_DIR}/palabos/externalLibraries/Eigen3")
target_compile_options(palabos PRIVATE ${HEMOFLOW_WARNING_FLAGS})
target_compile_definitions(palabos PUBLIC $<$<CONFIG:Debug>:PLB_DEBUG>)

if(ENABLE_MPI)
    # Skip the deprecated MPI C++ bindings (as upstream Palabos does)
    target_compile_definitions(palabos PUBLIC PLB_MPI_PARALLEL OMPI_SKIP_MPICXX MPICH_SKIP_MPICXX)
    target_link_libraries(palabos PUBLIC MPI::MPI_CXX)
endif()
if(WIN32 OR CYGWIN)
    target_compile_definitions(palabos PUBLIC PLB_WINDOWS)
else()
    target_compile_definitions(palabos PUBLIC PLB_USE_POSIX)
endif()
if(APPLE)
    target_compile_definitions(palabos PUBLIC PLB_MAC_OS_X)
endif()

###############################################################################
# HemoFlow executable
###############################################################################

add_executable(hemoFlow
    hemoFlow.cpp
    geometry.cpp
    opening.cpp
    io.cpp
    external/cnpy/cnpy.cpp)
set_target_properties(hemoFlow PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
target_include_directories(hemoFlow PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(hemoFlow SYSTEM PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/external/cnpy")
target_compile_options(hemoFlow PRIVATE ${HEMOFLOW_WARNING_FLAGS})
target_link_libraries(hemoFlow PRIVATE palabos ZLIB::ZLIB)

if(HEMOFLOW_WITH_HDF5)
    target_compile_definitions(hemoFlow PRIVATE HEMOFLOW_WITH_HDF5 ${HDF5_DEFINITIONS})
    target_include_directories(hemoFlow SYSTEM PRIVATE
        ${HDF5_INCLUDE_DIRS}
        "${CMAKE_CURRENT_SOURCE_DIR}/external/HighFive/include")
    target_link_libraries(hemoFlow PRIVATE ${HDF5_C_LIBRARIES} ${HDF5_LIBRARIES})
endif()

install(TARGETS hemoFlow RUNTIME DESTINATION bin)

message(STATUS "HemoFlow: MPI=${ENABLE_MPI} HDF5 writer=${HEMOFLOW_WITH_HDF5} C++${CMAKE_CXX_STANDARD}")
