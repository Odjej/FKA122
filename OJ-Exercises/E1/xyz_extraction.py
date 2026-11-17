from ase.build import molecule, bulk
from ase.io import write
import numpy as np
from ase.data.pubchem import pubchem_atoms_search

# 1. Extract a built-in molecule (e.g., CO2, H2O, CH4, C6H6)
co2 = molecule('CO2')
co2.arrays['masses'] = co2.get_masses()

# Write to Extended XYZ format
write('co2.xyz', co2, format='extxyz')


