# apps/CMakeLists.txt

# Example for creating one executable named 'run_analysis'
add_executable(elspectro elspectro.cpp)

# Link the executable against your libraries.
# It automatically gets the necessary include paths and flags from them.
target_link_libraries(elspectro PRIVATE
    elSpectro
    # JPACPHOTO is already linked via elSpectro, but linking again is fine
    JPACPHOTO
    # You may not need to link ROOT libraries again if elSpectro handles it,
    # but it's safe to be explicit.
    ROOT::Core
    ROOT::Rint
)

# Install the executable
install(TARGETS elspectro DESTINATION ${CMAKE_INSTALL_BINDIR})

# --- Add other executables here following the same pattern ---
# add_executable(another_app another_main.cpp)
# target_link_libraries(another_app PRIVATE elSpectro)
# install(TARGETS another_app DESTINATION ${CMAKE_INSTALL_BINDIR})
