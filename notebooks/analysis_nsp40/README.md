# NSP40 Analysis Notebooks

This folder contains Jupyter notebooks for processing and analyzing simulation data related to simulations with an array size of 40 spacers.

## Contents

- **Analysis_Simulations.ipynb**  
    Provides a general overview and summary of the simulation results.

- **Simulations_Data_Beta_1.2.ipynb**  
    Processes simulation output `.txt` files to extract curve properties from beta output and saves the results as `.csv` files in the `ANALYSIS` subfolder.

- **Simulations_Data_Fitness.ipynb**  
    Similar to the above, but focuses on extracting fitness-related properties from simulation outputs and storing them as `.csv` files in `ANALYSIS`.

- **Analysis_Simulations_x.ipynb**  
    Analyze the effect of input parameter (a, g and pEndemic) on output variable x ( age, beta or fitness).

## Folder Structure

```
analysis_nsp10/
├── data/                # Contains processed CSV and TXT files
├── Simulations_Data_Beta_1.2.ipynb
├── Simulations_Data_Fitness.ipynb
├── Analysis_Simulations.ipynb
├── ... (other analysis notebooks)
```

## Usage

1. Run the data processing notebooks (Simulations_Data_x) to generate CSV or TXT files in the `ANALYSIS` folder.
2. Use the analysis notebooks to explore how input parameters affect simulation outcomes.

## Versions

The notebooks ended in 1.1 or 1.2 come from older versions that you can find in the repository folder for older versions (see README from repository).