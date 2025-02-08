# Combined Hansen Solubility Parameters and GCVOL Calculator

This Python tool calculates Hansen Solubility Parameters (HSP) using the GCVOL method for molar volume calculations. It combines and extends the Python code presented in "Pencil and Paper Estimation of Hansen Solubility Parameters" (Mathieu, ACS Omega, 2018) with additional functionality for molar volume estimation.

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
1. `GCVOL_Gruppebeiträge.xlsx`: Contains GCVOL group contribution parameters
2. `Database2.xlsx`: Contains substance database

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

## Credits

This implementation is based on the methodology described in:
- Mathieu, D. "Pencil and Paper Estimation of Hansen Solubility Parameters". ACS Omega, 2018, 3 (12)
- Ihmels, E. C.; Gmehling, J. "Extension and Revision of the Group Contribution Method GCVOL for the Prediction of Pure Compound Liquid Densities". Ind. Eng. Chem. Res., 2002, 42 (2)
