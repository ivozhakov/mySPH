# mySPH

A standalone, self-written solver for the Smoothed Particle Hydrodynamics (SPH) method. SPH is a computational technique used to simulate fluid flows, gases, and other continuum media by representing the system as interacting particles.

## Overview

This project implements a basic SPH solver in C++. It is designed to be lightweight and independent, without relying on external libraries beyond standard C++. The solver can model hydrodynamic phenomena using a mesh-free approach, solving equations like Navier-Stokes through particle interactions.

Key components:
- **SPH Algorithms**: Core logic for particle smoothing, density calculation, pressure, viscosity, and time-stepping.
- **Examples**: Demonstration scenarios for testing and visualizing simulations (e.g., simple fluid flows).

## Features

- Mesh-free particle-based simulation.
- Efficient C++ implementation for performance.
- Modular structure with libraries in `sph_libs/` and testable examples in `examples/`.
- Supports basic fluid dynamics simulations.

## Requirements

- C++ compiler (e.g., g++ or clang++ with C++11 or later).
- Make (for building via Makefile).

No external dependencies required.

## Installation

1. Clone the repository:
   ```
   git clone https://github.com/ivozhakov/mySPH.git
   cd mySPH
   ```

2. Build the project:
   ```
   make
   ```

This will compile the libraries and examples.

## Usage

After building, you can run examples from the `examples/` directory. For instance:

```
./examples/simple_fluid_simulation
```

To integrate into your own project:
- Include headers from `sph_libs/`.
- Link against the compiled SPH library.

### Basic Code Example

```cpp
#include "sph_libs/sph_solver.h"

int main() {
    SPHSolver solver;
    solver.initializeParticles(/* parameters */);
    solver.simulate(/* timesteps */);
    return 0;
}
```

(Adjust parameters based on your simulation needs; refer to source code for details.)

## Directory Structure

- `sph_libs/`: Core SPH library code (classes for particles, kernels, etc.).
- `examples/`: Sample simulations and tests.
- `.clang-format`: Code formatting configuration.
- `.gitignore`: Git ignore file.

## Contributing

Contributions are welcome! Feel free to submit pull requests for improvements, bug fixes, or additional features. Please ensure code adheres to the `.clang-format` style.

## License

Free

## Acknowledgments

- Inspired by standard SPH literature and implementations.
- Author: [ivozhakov](https://github.com/ivozhakov)
