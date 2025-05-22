# Main Simulations – CRISPR Spacer Dynamics

This folder contains C++ scripts to simulate CRISPR-Cas adaptive immunity in microbial populations. The simulations model how spacer acquisition and virus-host dynamics evolve over time under various parameters.

## Contents

- `simulate_spacers1.2.cpp`: Original version of the simulation.
- `simulate_spacers1.3.cpp`: Modified version that supports extreme conditions like g=60 and Nsp=40, where efficacy values become very high, by making the logarithm.

## Compilation

Use a C++ compiler like `g++`. In my case I compiled the scripts using:

```bash
g++ -O3 simulate_spacers1.2.cpp -o simulate_spacers
g++ -O3 simulate_spacers1.3.cpp -o simulate_spacers
```

## Execution

To execute this script we need three input parameters: 
- `a`: Selection pressure (0.0–1.0).
- `g`: Time where endemic and epidemic incidences cross (0–60).
- `pEndemic`: Initial proportion of endemic spacers (0.0–1.0).

The command to execute is:
```bash
./executable parameter_a parameter_g parameter_pEndemic
```
As an example:
```bash
./simulate_spacers.cpp 0.1 0 0.2 
```

## Different fixed conditions

The number of generations, the array size, and the number of replicates must be changed in the `simulate_spacers1.2.cpp` or `simulate_spacers1.3.cpp` files in the following lines:
#define MAXGEN 5000	    // Number of generations to simulate
#define NUMSPC 40	    // Number of spacers in the array
#define NUMREPS 5000	// Number of replicates
