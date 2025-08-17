# J-GenHel: A Julia-based Helicopter Flight Dynamics Simulation

J-GenHel is a comprehensive, nonlinear helicopter flight dynamics simulation model written in the Julia programming language. It is based on the General Helicopter (GenHel) flight dynamics model from NASA, and the specific implementation is representative of a utility helicopter similar to a UH-60 Black Hawk.

This code was developed by Dr. Umberto Saetti at Auburn University and is based on the work presented in the paper "Flight Simulation and Control using the Julia Language" (see [How to Cite](#how-to-cite)).

## Features

*   **High-Fidelity Model:** Implements a 6-DoF rigid-body dynamic model of the fuselage.
*   **Nonlinear Aerodynamics:** Utilizes nonlinear aerodynamic lookup tables for the fuselage, rotor blades, and empennage.
*   **Detailed Rotor Dynamics:** Includes rigid flap and lead-lag rotor blade dynamics.
*   **Advanced Inflow Model:** Features a three-state Pitt-Peters inflow model for the main rotor and a one-state model for the tail rotor.
*   **Complete Analysis Suite:** Provides tools for:
    *   Trimming the aircraft at various flight conditions.
    *   Linearizing the model to obtain state-space representations.
    *   Model order reduction.
    *   Spectral analysis (eigenvalue plots).
    *   Frequency response analysis.
    *   Time-domain simulation (both open-loop and closed-loop with a dynamic inverse controller).

## Getting Started

### Prerequisites

You need to have [Julia](https://julialang.org/downloads/) installed on your system.

The following Julia packages are required:
*   `LinearAlgebra`
*   `SparseArrays`
*   `ControlSystems`
*   `MAT`
*   `Interpolations`
*   `ApproxFun`
*   `Plots`
*   `JLD`

### Installation

1.  Clone this repository.
2.  Open a Julia REPL and navigate to the project directory.
3.  Install the required packages. You can enter the package manager by typing `]` in the REPL and then run:
    ```julia
    add LinearAlgebra SparseArrays ControlSystems MAT Interpolations ApproxFun Plots JLD
    ```

### Running the Simulation

The main script `src/main.jl` runs a complete workflow of trimming, analysis, and simulation. To run it, start a Julia session, navigate to the `src` directory, and execute:

```julia
include("main.jl")
```

The script will generate several plots showing the system's eigenvalues, frequency responses, and time simulation results. It will also save output data to the `output/` directory.

## Project Structure

*   `src/`: Contains all Julia source files (`.jl`).
    *   `main.jl`: The main entry point to run the simulation and analysis.
    *   `GenHel.jl`: The core function defining the helicopter's equations of motion.
    *   `H60_constants.jl`: Defines the physical constants for the UH-60-like model.
    *   Other files contain various helper functions (e.g., for control mixing, trimming, linearization).
*   `data/`: Contains data files (`.mat`, `.jld`) with aerodynamic tables, rotor data, and other model parameters.
*   `output/`: The default directory for saving simulation results and generated data.

## How to Cite

If you use this code in your research, please cite the following publication:

> Saetti, U., and Horn, J. F., "Flight Simulation and Control using the Julia Language", AIAA Scitech Forum, San Diego, CA, Jan 3-7, 2022.
> DOI: [https://arc.aiaa.org/doi/10.2514/6.2022-2354](https://arc.aiaa.org/doi/10.2514/6.2022-2354)

The original mathematical model for the UH-60 is described in:

> Howlett, J. J., “UH-60A Black Hawk Engineering Simulation Program. Volume 1: Mathematical Model,” Tech. rep. NASA-CR-166309, 1980.

## License

This project is licensed under the MIT License. See the [LICENSE](src/LICENSE) file for details.

## Author

*   **Dr. Umberto Saetti**
    *   Assistant Professor, Department of Aerospace Engineering, Auburn University
    *   Email: saetti@auburn.edu
