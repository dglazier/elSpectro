cmake_minimum_required(VERSION 2.9)
project(elSpectro)

	
set(CMAKE_INSTALL_BINDIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(CMAKE_CXX_FLAGS "-fPIC -O3")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

#include($ENV{HEPMC3}/share/HepMC3/cmake/HepMC3Config.cmake)
#find_package(HepMC3 REQUIRED)

find_package(ROOT REQUIRED MathMore)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}) 
#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/core)
add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/core)
