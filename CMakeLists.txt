# 1. Use a modern CMake version and define the project
cmake_minimum_required(VERSION 3.15)
project(elSpectro VERSION 1.0 LANGUAGES CXX)

# 2. Use standard CMake modules for installation paths
# This replaces all the manual IF/SET logic for CMAKE_INSTALL_PREFIX
include(GNUInstallDirs)

# Set a default installation prefix if the user hasn't provided one.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    if($ENV{ELSPECTRO})
        # If $ELSPECTRO is set, use it as the default install path.
        set(CMAKE_INSTALL_PREFIX "$ENV{ELSPECTRO}" CACHE PATH "Default install path from $ELSPECTRO env var" FORCE)
        message(STATUS "Install prefix not set by user. Using default from $ELSPECTRO: ${CMAKE_INSTALL_PREFIX}")
    else()
        # If $ELSPECTRO is not set, default to a local 'install' directory.
        set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/" CACHE PATH "Default install path" FORCE)
        message(STATUS "Install prefix not set by user. Using default local path: ${CMAKE_INSTALL_PREFIX}")
    endif()
endif()

# 3. Find dependencies using modern CONFIG mode
# This provides imported targets like ROOT::Core
find_package(ROOT CONFIG REQUIRED COMPONENTS MathMore RooFit GenVector EG)

# 4. Handle the jpacPhoto submodule
# This is simpler than the custom function and adds the directory if it exists.
# This will create a 'JPACPHOTO' library target if the submodule is present.
# This logic checks for an environment variable first, then falls back to the submodule.
if($ENV{JPACPHOTO} AND IS_DIRECTORY $ENV{JPACPHOTO})
    # Store the path in a variable for subdirectories to use
    set(JPACPHOTO_SOURCE_DIR $ENV{JPACPHOTO})
    message(STATUS "Using jpacPhoto from environment variable: ${JPACPHOTO_SOURCE_DIR}")
    add_subdirectory(${JPACPHOTO_SOURCE_DIR})
    
elseif(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/jpacPhoto)
    # Store the path in a variable for subdirectories to use
    set(JPACPHOTO_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/jpacPhoto)
    message(STATUS "Using default jpacPhoto submodule.")
    add_subdirectory(${JPACPHOTO_SOURCE_DIR})
    
else()
    message(WARNING "jpacPhoto dependency not found.")
endif()

# 5. Add the core library subdirectory
# This defines the 'elSpectro' library target.
add_subdirectory(core)

# 6. Add the executables
# It's better practice to have executables in a separate directory (e.g., 'apps')
# instead of globbing all .cpp files. This avoids creating unwanted executables.
# This will define your main application targets.
add_subdirectory(apps)
