cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(MICROCERN C Fortran)

include_directories(. kerngen/.)

file(GLOB SOURCES "*.F" "*.c")

if(${CMAKE_SYSTEM} MATCHES Darwin)
    add_definitions(-DCERNLIB_LXIA64 -DCERNLIB_BLDLIB -DCERNLIB_CZ -DCERNLIB_PPC)
endif(${CMAKE_SYSTEM} MATCHES Darwin)
if(${CMAKE_SYSTEM} MATCHES Linux)
    add_definitions(-DCERNLIB_LINUX -DCERNLIB_BLDLIB -DCERNLIB_CZ)
endif(${CMAKE_SYSTEM} MATCHES Linux)

# not on Darwin
if(NOT ${CMAKE_SYSTEM} MATCHES Darwin)
    file(GLOB REMOVE "lnblnk.c")
    list(REMOVE_ITEM SOURCES ${REMOVE})
endif(NOT ${CMAKE_SYSTEM} MATCHES Darwin)

add_library(MICROCERN SHARED ${SOURCES})

install(TARGETS MICROCERN
        LIBRARY DESTINATION lib)

if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(MICROCERN PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)
