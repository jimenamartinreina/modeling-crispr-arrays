# Legacy Simulations – CRISPR Spacer Dynamics

This folder contains older versions of C++ scripts used to simulate CRISPR-Cas adaptive immunity in microbial populations. These scripts are preserved for historical reference and do not reflect the latest updates or improvements.

## Contents

- `simulate_spacers1.0.cpp`: The first code for simulating spacers' dynamics. Done by my tutor, Jaime Iranzo, in 2022.
- `simulate_spacers1.1.cpp`: A modified version of the original script that corrected bugs in some parts of the code.

## Purpose

These scripts are retained as legacy versions to document the evolution of the simulation codebase. They are not actively maintained and may lack features or optimizations present in newer versions.

## Usage

While these scripts can still be compiled and executed, users are encouraged to refer to the latest versions for improved performance and accuracy. If needed, the compilation and execution instructions are similar to those provided in the main simulation documentation.

## Notes

- The number of generations, array size, and number of replicates can be adjusted directly in the source code by modifying the following lines:
    ```cpp
    #define MAXGEN 5000    // Number of generations to simulate
    #define NUMSPC 40      // Number of spacers in the array
    #define NUMREPS 5000   // Number of replicates
    ```
- These scripts are primarily intended for archival purposes and may not be suitable for current research needs.


