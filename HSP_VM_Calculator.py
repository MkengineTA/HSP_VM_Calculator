#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Combined HSP and GCVOL calculator
Calculates Hansen Solubility Parameters using GCVOL method for molar volume
"""

import fileinput
from math import sqrt
from collections import Counter, OrderedDict, defaultdict
import pandas as pd
import numpy as np
from rdkit import Chem

# Parameters for molar refractivity calculation
# Format: (code, volume_increment, refractivity_increment)
REFRACTIVITY_PARAMS = [
    ('B30', 7.61, 2.91),
    ('Br10', 27.52, 8.44),
    ('C10', 14.68, 2.32),
    ('C20', 9.89, 4.01),
    ('C21', 22.77, 4.58),
    ('C30', -0.00, 3.15),
    ('C30a', -0.00, 3.48),
    ('C31', 13.23, 4.65),
    ('C31a', 13.23, 4.46),
    ('C32', 27.17, 5.38),
    ('C40', -8.40, 2.60),
    ('C41', 4.46, 3.48),
    ('C42', 16.57, 4.60),
    ('C43', 29.58, 5.74),
    ('Cl10', 24.74, 5.87),
    ('F10', 17.63, 1.06),
    ('N10', 14.09, 1.55),
    ('N20', 7.42, 2.60),
    ('N20a', 7.42, 2.44),
    ('N21', 18.14, 3.55),
    ('N30', -3.08, 2.89),
    ('N30a', -3.08, 2.85),
    ('N31', 7.74, 3.69),
    ('N31a', 7.74, 3.72),
    ('N32', 17.81, 4.60),
    ('O10', 14.89, 1.84),
    ('O20', 6.25, 1.55),
    ('O20a', 6.25, 0.71),
    ('O21', 11.78, 2.51),
]

# Convert to dictionaries for easier access
ri = dict((p[0], p[2]) for p in REFRACTIVITY_PARAMS)

# Parameters for polar component calculation
params_p = {
    'N(1)': 2783,
    'N(2)': 8235,
    'O(0)': 1603,
    'O(1)': 4125,
    'Cl(0)': 1637,
    'C=O': 7492,
    'COOH': -5494,
    'COinAmide': 15972,
    'Carbonate': 19019,
    'Ester': 3653,
    'C#N': 16056,
    'NitroOnC': 13276,
    'O=P': 20310,
}

# Parameters for hydrogen bonding component
H_BOND_PARAMS = {
    'HC': 24.5,
    'HN': -1576,
    'HNamide': 5060,
    'H2N': 5484,
    'HO': 16945,
    'HO_COOH': 7094,
    'N': 3252,
    'O': 1980,
    'X': 412,
}

# SMARTS patterns for polar groups
POLAR_GROUPS = OrderedDict([
    ("COOH", ("[CX3](=O)[OX2H1]", (1, 2))),
    ("NitroOnC", ("[#6][$([NX3](=O)=O),$([NX3+](=O)[O-])]", (1,))),
    ("Carbonate", ("[OX2;R][CX3;R](=O)[OX2]", (0, 1, 3))),
    ("Ester", ("[OX2][CX3]=O", (1,))),
    ("COinAmide", ("[NX3][CX3](=[OX1])[#6]", (1,))),
    ("SO2", ("O=S=O", (1,))),
])


class GCVOLCalculator:
    def __init__(self, gcvol_file, database_file):
        """Initialize GCVOL calculator with input files"""
        self.gcvol_df = pd.read_excel(gcvol_file, header=None)
        self.database_df = pd.read_excel(database_file, header=None)
        self.setup_gcvol_data()

    def setup_gcvol_data(self):
        """Process GCVOL data for calculations"""
        # Extract the coefficient data
        self.gcvol = self.gcvol_df.iloc[3:, 2:].copy()
        self.gcvol = self.gcvol.replace({pd.NA: np.nan})
        self.gcvol = self.gcvol.astype(float)

        # Store coefficients
        self.A = self.gcvol.iloc[:, 0]
        self.B = self.gcvol.iloc[:, 1]
        self.C = self.gcvol.iloc[:, 2]

        # Create substance mapping
        self.substance_columns = {}
        for col in range(5, self.gcvol_df.shape[1]):
            substance = self.gcvol_df.iloc[2, col]
            if isinstance(substance, str):
                self.substance_columns[substance] = col

    def calculate_molar_volume(self, substance_name, temperature):
        """Calculate molar volume using GCVOL method"""
        try:
            if substance_name not in self.substance_columns:
                return None

            col_idx = self.substance_columns[substance_name]
            group_counts = self.gcvol_df.iloc[3:, col_idx]

            v_mol = 0
            for i, count in enumerate(group_counts):
                if count > 0:
                    v_j = count * (
                            self.A[i] +
                            0.001 * temperature * self.B[i] +
                            0.00001 * temperature ** 2 * self.C[i]
                    )
                    v_mol += v_j

            return v_mol * 1e-6  # Convert to m³/mol

        except Exception as e:
            print(f"Error calculating molar volume for {substance_name}: {str(e)}")
            return None

def get_atom_code(atom):
    """ return an atom code consistent with the keys of the 'Vinc' dictionary """
    # be careful to take into account nD = number of deuterium atoms !
    nD = len([x for x in atom.GetNeighbors() if x.GetMass()>2 and x.GetSymbol()=='H'])
    #if nD > 0:
    #    print('nD %u %s' %(nD, arg))
    code = atom.GetSymbol() + str(atom.GetTotalDegree()) + str(atom.GetTotalNumHs()+nD)
    code += 'a'*int(atom.GetIsAromatic())
    return code

def get_ring_descriptors(mol, maxi=6, mini=5):
    """ return dict of ring descriptors for molecule provided as input """
    dic = Counter()
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        size = len(ring)
        if size > maxi:
            label = 'R>'+str(maxi)
        elif size < mini:
            label = 'R<'+str(mini)
        else:
            label = 'R%u' %len(ring)
        dic[label] += 1
        # contribute also +1 aromatic ring ?
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if all(at.GetIsAromatic() for at in atoms):
            dic['aromat'] += 1
    return dic

def GetBChar(bond):
    if bond.GetBondType() == Chem.rdchem.BondType.AROMATIC: return '~'
    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE: return '='
    if bond.GetBondType() == Chem.rdchem.BondType.TRIPLE: return '#'
    return '-'


def get_nH(atom):
    nD = len([x for x in atom.GetNeighbors() if x.GetMass() > 2 and x.GetSymbol() == 'H'])
    nH = nD + atom.GetTotalNumHs()
    return nH


def isNitroN(at, mol):
    if at.GetSymbol() != "N": return False
    if at.GetTotalDegree() != 3: return False
    voisins = at.GetNeighbors()
    Os = [a for a in voisins if a.GetSymbol() == "O"]
    O1s = [a for a in Os if a.GetTotalDegree() == 1]
    return len(O1s) > 1


def isAmideN(at, mol):
    amideSMARTS = "[NX3][CX3](=[OX1])[#6]"
    amidePattern = Chem.MolFromSmarts(amideSMARTS)
    N_indices = [t[0] for t in mol.GetSubstructMatches(amidePattern)]
    return at.GetIdx() in N_indices


def inCOOH(at, mol):
    acSMARTS = "[CX3](=O)[OX2H1]"
    acPattern = Chem.MolFromSmarts(acSMARTS)
    OH_indices = [t[2] for t in mol.GetSubstructMatches(acPattern)]
    return at.GetIdx() in OH_indices


def get_polar_groups(mol):
    global POLARGROUPS
    by_type = defaultdict(list)
    counted_indices = set()
    # first count complex polar groups
    for group_name in POLARGROUPS:
        group, positions = POLARGROUPS[group_name]
        pattern = Chem.MolFromSmarts(group)
        tuples = mol.GetSubstructMatches(pattern)
        for tup in tuples:
            if set(tup) & counted_indices: continue
            counted_indices |= set(tup)
            by_type[group_name].append(1)
    # count insaturated polar bonds
    for bond in mol.GetBonds():
        order = GetBChar(bond)
        if order in ("#", "=", "~"):
            abeg, aend = bond.GetBeginAtom(), bond.GetEndAtom()
            symbols = sorted([abeg.GetSymbol(), aend.GetSymbol()])
            if symbols[0] == symbols[1]: continue  # Skip C=C and C~C bonds
            tup = (abeg.GetIdx(), aend.GetIdx())
            if set(tup) & counted_indices: continue
            bondsymbol = order.join(symbols)
            counted_indices |= set(tup)
            by_type[bondsymbol].append(order)
    # count saturated heteroatoms
    for hetat in mol.GetAtoms():
        idx = hetat.GetIdx()
        if idx in counted_indices: continue
        coo = hetat.GetTotalDegree()
        symbol = hetat.GetSymbol()
        if symbol == "C": continue
        if symbol == "P" and coo > 3: continue
        name = "%s(%u)" % (symbol, get_nH(hetat))
        if name in ("N(0)", "F(0)"): continue
        counted_indices.add(idx)
        by_type[name].append(idx)
    return dict((group_name, len(by_type[group_name])) for group_name in by_type)


def get_HSPp(mol):
    d = get_polar_groups(mol)
    if set(d) - set(params_p): return None
    eP = sum(params_p[k] * d[k] for k in d)
    return eP


# FOR HYDROGEN-BONDING (H) HSP COMPONENT

params_h = {
    'HC': 24.5,
    'HN': -1576,
    'HNamide': 5060,
    'H2N': 5484,
    'HO': 16945,
    'HO_COOH': 7094,
    'N': 3252,
    'O': 1980,
    'X': 412,
}


def get_dic_h(mol):
    dic = defaultdict(int)
    for at in mol.GetAtoms():
        symbol = at.GetSymbol()
        coo = at.GetTotalDegree()
        voisins = [v for v in at.GetNeighbors()]
        if symbol in ("N", "O"):
            if not isNitroN(at, mol):
                dic[symbol] += 1
        if symbol in ("F", "Cl", "Br", "I"):
            dic["X"] += 1
        nH = get_nH(at)
        if nH == 0: continue
        if symbol == "C":
            dic["HC"] += nH
            continue
        if symbol == "N":
            if isAmideN(at, mol):
                dic["HNamide"] += nH
                continue
            if nH == 2:
                dic['H2N'] += 2
            elif nH == 1:
                dic["HN"] += 1
            continue
        if symbol == "O":
            if inCOOH(at, mol):
                dic["HO_COOH"] += nH
            else:
                dic["HO"] += nH
    return dic


def get_HSPh(mol):
    d = get_dic_h(mol)
    if set(d) - set(params_h): return None
    eH = sum(params_h[k] * d[k] for k in d)
    return eH


# UTILITY FOT OUTPUT

def tostr(energy, volume=None):
    if None in (energy, volume): return "  xxxx "
    if volume == 0: return " %5.1f " % energy
    return " %5.1f " % sqrt(max(energy / volume, 0))


def get_Vm_RD(mol, gcvol_calc, substance_name, temperature):
    """Calculate molar volume and refractivity"""
    # Calculate molar refractivity
    atoms = Counter(get_atom_code(atom) for atom in mol.GetAtoms())
    missing_atoms = set(atoms) - set(ri)
    if missing_atoms:
        return None, None

    RD = sum(atoms[k] * ri[k] for k in atoms)

    # Calculate molar volume using GCVOL
    Vm = gcvol_calc.calculate_molar_volume(substance_name, temperature)

    return Vm, RD


def format_hsp_value(energy, volume=None):
    """Format HSP component value for output

    Args:
        energy: Energy value to format
        volume: Volume value for normalization (optional)

    Returns:
        Formatted string with HSP component value
    """
    if None in (energy, volume):
        return "  xxxx "
    if volume == 0:
        return f" {energy:5.1f} "
    return f" {sqrt(max(energy / volume, 0)):5.1f} "


def main():
    # Initialize GCVOL calculator
    gcvol_calc = GCVOLCalculator("GCVOL_Group_Contributions.xlsx", "Database.xlsx")

    for entry in fileinput.input():
        # Parse input line from Database
        data = entry.strip().split(',')
        substance_name = data[0]
        smiles = data[1]
        temperature = float(data[2])

        # Create RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"{substance_name} => None")
            continue

        # Calculate molar volume and refractivity
        Vm, RD = get_Vm_RD(mol, gcvol_calc, substance_name, temperature)

        # Calculate HSP components
        if Vm is None or RD is None:
            hspD = None
        else:
            hspD = sqrt(93.8 + (2016 + 75044. / Vm) * (RD / Vm) ** 2)

        eP = get_HSPp(mol)
        eH = get_HSPh(mol)

        # Format output using format_hsp_value instead of tostr
        hspD_str = format_hsp_value(hspD, 0)
        hspP_str = format_hsp_value(eP, Vm)
        hspH_str = format_hsp_value(eH, Vm)

        if Vm is not None:
            line = f"{substance_name:<12}{hspD_str}{hspP_str}{hspH_str} {Vm:7.1f}"
        else:
            line = f"{substance_name:<12}{hspD_str}{hspP_str}{hspH_str}   None"
        print(line)


if __name__ == "__main__":
    main()