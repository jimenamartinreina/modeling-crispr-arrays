#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 10:42:01 2024

@author: jimenamartinreina

This is a Python script to create a table of different values for three parameters (alpha, gamma
and pEndemic), based on how many different values you want to study.

It was designed to be used with the model in "simulate_spacers".
"""
import numpy as np
from itertools import product

# In logspace the first number will be base**start (10**-1 = 0.1)
a_numbers = np.logspace(-1, 0, 10)
# g>0
g_numbers = np.linspace(0, 60, 2)
# 0<pEndemic<1
pEndemic_numbers = np.linspace(0, 0.9, 10)

combinations = list(product(a_numbers, g_numbers, pEndemic_numbers))

# We shuffle the list, so it does not influence later that some parameters make the simulationn longer
#combinations = np.random.permutation(combinations)

# Create an output file with these values
with open('parameters.txt', 'w') as f:
    # Write distinct combinations directly
    for line in combinations:
        f.write(" ".join(f"{val:.2f}" for val in line) + "\n")
        
print("See table of parameters in 'parameters.txt'.")
