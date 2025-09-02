# elSpectro

`elSpectro` is a C++ framework designed for the simulation of particle physics events, with a focus on electron scattering experiments. It provides a set of core libraries for defining particles, decay models, and production processes, allowing users to build and run custom event generators.


**[View the Interactive Guide & Full Documentation](https://dglazier.github.io/elSpectro/)**

## Quick Build
```bash
git clone --recurse-submodules https://github.com/dglazier/elSpectro
...

***

## Core Concepts

The project is architecturally split into two main components:

* **`core` Library (`libelSpectro`):** This is the heart of the framework. It is a shared library that contains all the fundamental physics classes and logic for defining particles, managing decays, and modeling reactions.
* **`apps` Executables:** These are standalone programs that use the `core` library. The primary executable, `elSpectro`, acts as a custom ROOT environment that simplifies running simulation scripts.

***

## Directory Structure

The repository is organized as follows:
elSpectro/
│

├── apps/              # Contains source code for executables (e.g., the elSpectro runner).

│

├── core/              # Source code for the main elSpectro shared library.

│

├── examples/          # Example ROOT scripts showing how to use the library.

│

├── jpacPhoto/         # A submodule dependency for physics amplitudes.

│

└── CMakeLists.txt     # The main CMake file that orchestrates the build.


***

## Class Structure

This section describes the main classes in the `core` library and how they relate to one another through inheritance (an "is-a" relationship) and composition (a "has-a" relationship).

### Inheritance Hierarchies (`is-a` relationships)

These chains define the different *types* of models and tools available.

* **Decay and Production Models**
    * `class DecayModel`
        * `class Bremsstrahlung : public DecayModel`
        * `class Formation : public DecayModel`
            * `class FormationQ2W : public Formation`
        * `class PhaseSpaceDecay : public DecayModel`
        * `class ProductionModel : public DecayModel`
            * `class TwoBodyProduction : public ProductionModel`
        * `class SDMEDecay : public DecayModel`

* **Final State Kinematics**
    * `class DecayVectors`
        * `class TwoBodyFlat : public DecayVectors`
            * `class TwoBodyEnvelope : public TwoBodyFlat`

* **Top-Level Physics Processes**
    * `class ProductionProcess`
        * `class ElectronScattering : public ProductionProcess`

* **Statistical Distributions**
    * `class Distribution`
        * *(All `Dist...` classes inherit from `Distribution`)*

* **Particle Types**
    * `class Particle`
        * `class DecayingParticle : public Particle`
        * `class CollidingParticle : public Particle`

* **Event Output Writers**
    * `class Writer`
        * `class HepMC3Writer : public Writer`
        * `class LundWriter : public Writer`

#### Standalone Classes

These classes do not inherit from the main hierarchies listed above.

* `class ExcitationSpectra`
* `class JpacTwoBody`
* `class TwoBodytEnvelope`
* ...and other utility classes.

### Compositional Relationships (How Objects Work Together)

This describes how different objects are assembled to define a complete decay process. The workflow is orchestrated by the **`DecayChannel`**.

* **`DecayingParticle`**: This is the top-level object. Each `DecayingParticle` instance **contains one or more `DecayChannel` objects**, representing all the ways it can decay.

* **`DecayChannel`**: This class represents a single, specific decay mode (e.g., $p \rightarrow \pi^0 + n$). It is the central coordinator for a decay and **contains collections** of different models:
    * A `vector` of **`DecayModel`** pointers.
    * A `vector` of **`DecayVectors`** pointers.

During an event, the `DecayChannel` uses these models in sequence:

1.  **`DecayModel` (The "What")**: A `DecayModel` is chosen to **define the final state particles**. Its job is to answer the question, "What particles does the parent decay into?"

2.  **`DecayVectors` (The "How")**: Once the daughter particles are known, a corresponding `DecayVectors` object is used to **generate the decay kinematics**. It answers the question, "What are the momenta and angles of the daughter particles?"

***

## Dependencies

To build and run `elSpectro`, you will need:

* A **C++17** compliant compiler (e.g., GCC, Clang).
* **CMake** (version 3.15 or newer).
* **ROOT** (version 6 or newer): Ensure the `root` executable is in your system's PATH.

***

## Building the Project

The project uses a standard out-of-source CMake build process.

#### 1. Clone the Repository

Clone the repository and its submodules using the `--recursive` flag.
```bash

	git clone --recursive [https://github.com/dglazier/elSpectro.git](https://github.com/dglazier/elSpectro.git)
	cd elSpectro
2. Configure with CMake
It's recommended to install the library locally into a directory within the project. The command below will configure the build system to install everything into a directory named install/.

```bash

cmake -S . -B build -DCMAKE_INSTALL_PREFIX=install
3. Build the Code
Compile the library and all applications using the following command.

```bash

	cmake --build build

4. Install the Library
This final step will copy the compiled library, executables, and headers into the install/ directory you specified.

```bash

	cmake --install build
After this step, you will have a self-contained installation in elSpectro/install/ with bin/, lib/, and include/ subdirectories.

Build Options
The build can be customized with the following options:

-DCMAKE_INSTALL_PREFIX=<path>: The standard CMake flag to override the installation location.

```bash

	# Configure to install at /opt/elSpectro
	cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/opt/elSpectro

$ELSPECTRO Environment Variable: If CMAKE_INSTALL_PREFIX is not set by the user, the build will default to the path specified by this environment variable.

$JPACPHOTO Environment Variable: To use an external build of the jpacPhoto dependency instead of the submodule, set this variable to its source directory.

Using elSpectro
There are two steps to using the software: setting up your environment and running a simulation script.

1. Setting Up Your Environment (One-Time Setup)
To allow your system to find the elSpectro executable and its required libraries from any directory, you need to update your shell's environment variables.

For bash or zsh users:
Run the following commands, replacing /path/to/your/install with the actual path you installed to.

```bash

	# Set a variable for your installation path for convenience
	export ELSPECTRO_INSTALL_DIR="/path/to/your/install"

	# Add the executable directory to your PATH
	export PATH="$ELSPECTRO_INSTALL_DIR/bin:$PATH"

	# Add the library directories to your LD_LIBRARY_PATH
	# You need to add both the elSpectro and jpacPhoto lib directories.
	export LD_LIBRARY_PATH="$ELSPECTRO_INSTALL_DIR/lib:$JPACPHOTO/lib:$LD_LIBRARY_PATH"
For tcsh or csh users:
Run the following commands, replacing /path/to/your/install with your actual installation path.

```tcsh
    # Set a variable for your installation path for convenience
    setenv ELSPECTRO_INSTALL_DIR "/path/to/your/install"

    # Add the executable directory to your path
    set path = ( $ELSPECTRO_INSTALL_DIR/bin $path )

    # Add the library directories to your LD_LIBRARY_PATH
    setenv LD_LIBRARY_PATH "${ELSPECTRO_INSTALL_DIR}/lib:${JPACPHOTO}/lib:${LD_LIBRARY_PATH}"
    
To make this change permanent, add the appropriate commands to your shell's startup file (e.g., ~/.bashrc for bash, ~/.zshrc for zsh, or ~/.cshrc for tcsh).

2. Writing a Simulation Script
Your simulation logic is written in a ROOT script file (e.g., my_simulation.C). The elSpectro executable automatically loads the necessary libraries, so you can directly use its classes.

```C++

	// my_simulation.C

	void my_simulation() {
    	// You can immediately access all the elSpectro tools.
    	   elSpectro::Particle proton("proton", 0.938);

	   std::cout << "elSpectro library is loaded and ready!" << std::endl;
	   std::cout << "Created a particle: " << proton.GetName() << std::endl;
	   std::cout << "Mass: " << proton.GetMass() << " GeV" << std::endl;
	   }
	   
3. Running a Simulation
Once your environment is configured, you can run your script with the elSpectro command from any directory.

```bash

	elSpectro my_simulation.C