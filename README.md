# Combined Hansen Solubility Parameters and GCVOL Calculator

This Python tool calculates Hansen Solubility Parameters (HSP) using the GCVOL method for molar volume calculations. It combines and extends the work presented in "Pencil and Paper Estimation of Hansen Solubility Parameters" (Mathieu, ACS Omega, 2018) with additional functionality for molar volume estimation.

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

## Usage

The program accepts input through standard input in CSV format with the following structure:
```
substance_name,SMILES,temperature
```

Example usage:
```bash
echo "ethanol,CCO,298.15" | python hsp_calculator.py
```

Output format:
```
substance_name  δD     δP     δH     Vm
```
where:
- δD: Dispersive component (MPa^0.5)
- δP: Polar component (MPa^0.5)
- δH: Hydrogen bonding component (MPa^0.5)
- Vm: Molar volume (cm³/mol)

## Methodology

### Molar Volume Calculation
The molar volume is calculated using the GCVOL method, which uses group contributions and temperature-dependent parameters:

V_mol = Σ(n_i * (A_i + B_i*T*10^-3 + C_i*T^2*10^-5))

where:
- n_i: Number of groups of type i
- A_i, B_i, C_i: GCVOL parameters
- T: Temperature in Kelvin

### Hansen Solubility Parameters
The three HSP components are calculated using:

1. Dispersive (δD): Based on molar refractivity and volume
2. Polar (δP): Using group contribution method
3. Hydrogen bonding (δH): Using group contribution method

## Credits

This implementation is based on the methodology described in:
- Mathieu, D. "Pencil and Paper Estimation of Hansen Solubility Parameters". ACS Omega, 2018, 3 (12)
- The GCVOL method for molar volume calculations

## License

[Add your chosen license here]

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Contact

[Add your contact information here]
