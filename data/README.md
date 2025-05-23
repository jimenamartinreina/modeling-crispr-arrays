# Data Folder README

This folder contains output files in `.txt` format generated from simulation runs. The data is organized into three main subfolders:

## Folder Structure

- **nsp10/**
    - Contains simulation results for an array size (NSP) of 10 spacers.
    - Only one simulation round was performed, with 10,000 replicates.
    - Data generated using the `simulate_spacers1.2` C++ script.

- **nsp40/**
    - Contains simulation results for an array size (NSP) of 40 spacers.
    - Includes two rounds of simulations (`DIR1` and `DIR2`) with the same parameters.
    - The `MEANS` subfolder contains merged outputs from `DIR1` and `DIR2`.
    - Data generated using the `simulate_spacers1.2` C++ script.

- **nsp40_g60/**
    - Similar structure to `nsp40`, with `DIR1`, `DIR2`, and `MEANS` subfolders.
    - Data generated using the `simulate_spacers1.3` C++ script.

For further details on simulation parameters or file formats, refer to the main project documentation.