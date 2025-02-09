# Combined Hansen Solubility Parameters and GCVOL Calculator

This Python tool calculates Hansen Solubility Parameters (HSP) using the GCVOL method for molar volume calculations. It extends and modifies the original Python code for HSP calculations from "Pencil and Paper Estimation of Hansen Solubility Parameters" (Mathieu, ACS Omega, 2018) by integrating the GCVOL method for improved molar volume estimation.

## Overview

The tool provides calculations for:
- Molar volume using the GCVOL method
- Hansen Solubility Parameters (HSP)
  - Dispersive component (δD)
  - Polar component (δP)
  - Hydrogen bonding component (δH)

## Requirements

- Python 3.x
- Required packages:
  - rdkit
  - pandas
  - numpy

## Input Files

The tool requires two Excel files:
1. `GCVOL_Group_Contributions.xlsx`: Contains GCVOL group contribution parameters (part of this repository)
2. `Database.xlsx`: Contains substance database

The tool requires Database.xlsx to contain your compounds data with the following columns:

Column A: Name/identifier of the compound
Column B: SMILES notation of the molecular structure
Column C: Temperature in Kelvin

Important Note:
The identifier in Column A must match exactly with a compound name in the GCVOL table. The tool uses this identifier to look up the pre-calculated group contributions - it does not perform automatic group contribution analysis of the molecule. You need to manually determine the groups present in your molecule and ensure they match with an entry in the GCVOL table.

## Methodology

### Molar Volume Calculation
The molar volume is calculated using the GCVOL method, which uses group contributions and temperature-dependent parameters:

V_mol = Σ(n<sub>i</sub> * (A<sub>i</sub> + B<sub>i</sub> * T * 10<sup>-3</sup> + C<sub>i</sub> * T<sup>2</sup> * 10<sup>-5</sup>))

where:
- n<sub>i</sub>: Number of groups of type i
- A<sub>i</sub>, B<sub>i</sub>, C<sub>i</sub>: GCVOL parameters
- T: Temperature in Kelvin

### Hansen Solubility Parameters
The three HSP components are calculated using:

1. Dispersive (δD): Based on molar refractivity and volume
2. Polar (δP): Using group contribution method
3. Hydrogen bonding (δH): Using group contribution method

## Attribution and Credits
This software is based on and extends two key methodologies:

### HSP Calculation Core:

- Originally described in: Mathieu, D. "Pencil and Paper Estimation of Hansen Solubility Parameters". ACS Omega, 2018, 3 (12)
- Core HSP calculation functions have been adapted from the original implementation
- The original code has been modified under academic fair use for research purposes


### GCVOL Method Integration:

- Based on: Ihmels, E. C.; Gmehling, J. "Extension and Revision of the Group Contribution Method GCVOL for the Prediction of Pure Compound Liquid Densities". Ind. Eng. Chem. Res., 2002, 42 (2)
- Implementation of temperature-dependent molar volume calculations
