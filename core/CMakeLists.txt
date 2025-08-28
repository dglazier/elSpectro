# 1. Automatically find all source (.cpp) and header (.h) files in this directory.
file(GLOB CORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB CORE_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# -- FIX STARTS HERE --
# 2. IMPORTANT: Remove the LinkDef file from the main header list to ensure it's added last.
list(REMOVE_ITEM CORE_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/ElSpectroLinkDef.h")
# -- FIX ENDS HERE --

# 3. Define the library target using the automatically found source files.
add_library(elSpectro SHARED ${CORE_SOURCES})

# 4. Generate the ROOT dictionary.
# The cleaned header list is passed first, and then the LinkDef file is added at the end.
ROOT_GENERATE_DICTIONARY(G__elSpectro
    ${CORE_HEADERS}
    OPTIONS "-I${JPACPHOTO_SOURCE_DIR}"
    LINKDEF ElSpectroLinkDef.h
)
target_sources(elSpectro PRIVATE G__elSpectro.cxx)

# 5. Set Target Properties (Modern way to set C++ standard and flags)
set_target_properties(elSpectro PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    POSITION_INDEPENDENT_CODE ON # Replaces the -fPIC flag
)

# 6. Link dependencies directly to the target
target_link_libraries(elSpectro PUBLIC
    JPACPHOTO
    ROOT::Core
    ROOT::Rint
    ROOT::RIO
    ROOT::RooFit
    ROOT::MathMore
    ROOT::EG
    ROOT::GenVector
)

# 6. Define public include directories for the target
target_include_directories(elSpectro PUBLIC
    ${JPACPHOTO_SOURCE_DIR}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

# 7. Set up installation rules for the library and its headers
install(TARGETS elSpectro
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libelSpectro_rdict.pcm
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)
