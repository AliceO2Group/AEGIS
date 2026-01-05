cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(AEGIS)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

add_subdirectory(MICROCERN)

add_subdirectory(TEPEMGEN)

add_subdirectory(GeneratorCosmics)

find_package(nlohmann_json REQUIRED)
add_subdirectory(GeneratorParam)

add_subdirectory(GeneratorSlowNucleons)

add_subdirectory(GeneratorSpectators)

if (DEFINED ENV{HIJING_ROOT})
   add_subdirectory(THijing)
endif (DEFINED ENV{HIJING_ROOT})

if (DEFINED ENV{DPMJET_ROOT})
   add_subdirectory(TDPMjet)
endif (DEFINED ENV{DPMJET_ROOT})

#Install(DIRECTORY Geometry gconfig DESTINATION share/Detectors/)
