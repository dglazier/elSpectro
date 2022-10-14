cmake_minimum_required(VERSION 2.9)
project(elSpectro)

#if -DCMAKE_INSTALL_PREFIX=/install/here use that 
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR} CACHE PATH "install to source" FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

IF(NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_BINDIR ${CMAKE_INSTALL_PREFIX}/bin)
  set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
ENDIF(NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

#if no default install dirs, install at source
if(NOT DEFINED CMAKE_INSTALL_BINDIR)
  set(CMAKE_INSTALL_BINDIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
endif(NOT DEFINED CMAKE_INSTALL_BINDIR)

if(NOT DEFINED CMAKE_INSTALL_LIBDIR)
  set(CMAKE_INSTALL_LIBDIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)
endif(NOT DEFINED CMAKE_INSTALL_LIBDIR)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_INSTALL_BINDIR})

message(INFO ${CMAKE_INSTALL_BINDIR} )

set(CMAKE_CXX_FLAGS "-fPIC -O3")

find_package(ROOT REQUIRED MathMore RooFit GenVector EG)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}) 
#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")

 
##INSTALL JPACPHOTO IF EXISTS AS SUBMODULE
function(my_add_subdirectory _dir)
  get_filename_component(_fullpath "${_dir}" REALPATH)
  if (EXISTS "${_fullpath}" AND EXISTS "${_fullpath}/CMakeLists.txt")
    get_property(_included_dirs GLOBAL PROPERTY GlobalAddSubdirectoryOnceIncluded)
    list(FIND _included_dirs "${_fullpath}" _used_index)
    if (${_used_index} EQUAL -1)
      set_property(GLOBAL APPEND PROPERTY GlobalAddSubdirectoryOnceIncluded "${_fullpath}")
      include_directories(${_dir})
      add_subdirectory(${_dir} ${ARGN})
      set(ENV{JPACPHOTO} ${PWD}/${_dir})
      message(INFO " using local jpacPhoto")
    endif()
  else()
    message(WARNING "my_add_subdirectory: Can't find ${_fullpath}/CMakeLists.txt")
  endif()
endfunction(my_add_subdirectory)

my_add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/jpacPhoto)


##INSTALL ELSPECTRO LIBRARY 
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/core)
add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/core)

##INSTALL ELSPECTRO BINARY
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/core/src/)
file(GLOB_RECURSE EXE_FILES "core/src/*.cpp")
foreach( exefile ${EXE_FILES} )
  get_filename_component( exename ${exefile} NAME_WE)
  add_executable( ${exename} ${exefile} )
  target_link_libraries( ${exename} elSpectro)
  target_link_libraries( ${exename} jpacPhoto)
  target_link_libraries( ${exename} ${ROOT_LIBRARIES} -lRooFit -lMathMore -lEG -lGenVector)
  #target_link_libraries( ${exename}  ROOT::Core ROOT::Rint ROOT::RIO ROOT::RooFit ROOT::MathMore ROOT::EG ROOT::GenVector )
endforeach( exefile ${EXE_FILES} )
