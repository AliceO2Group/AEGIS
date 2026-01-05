cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(THijing C CXX Fortran)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS EG)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

link_directories($ENV{HIJING_ROOT}/lib)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/.)
ROOT_GENERATE_DICTIONARY(G__THijing THijing.h LINKDEF THijingLinkDef.h)

#---Create a shared library with geneated dictionary
add_library(THijing SHARED THijing.cxx openDecayTableHijing.F G__THijing.cxx)
target_link_libraries(THijing ${ROOT_LIBRARIES} HIJING MICROCERN)

install(TARGETS THijing
        LIBRARY DESTINATION lib)

if (${ROOT_VERSION} VERSION_GREATER "6.0")
   install(
      FILES
      ${CMAKE_CURRENT_BINARY_DIR}/libTHijing_rdict.pcm
      ${CMAKE_CURRENT_BINARY_DIR}/libTHijing.rootmap
      DESTINATION lib)
endif (${ROOT_VERSION} VERSION_GREATER "6.0")

if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(THijing PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)
