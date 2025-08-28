# elSpectro

`elSpectro` is a C++ framework designed for the simulation of particle physics events, with a focus on electron scattering experiments. It provides a set of core libraries for defining particles, decay models, and production processes, allowing users to build and run custom event generators.

***

## Core Concepts

The project is architecturally split into two main components:

* **`core` Library (`libelSpectro`):** This is the heart of the framework. It is a shared library that contains all the fundamental physics classes and logic for defining particles, managing decays, and modeling reactions.
* **`apps` Executables:** These are standalone programs that use the `core` library. The primary executable, `elSpectro`, acts as a custom ROOT environment that simplifies running simulation scripts.

***

## Directory Structure

The repository is organized as follows: