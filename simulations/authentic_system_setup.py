from pymatgen.io.openmm.generators import OpenMMSolutionGen

from pymatgen.core import Molecule

import numpy as np

EC = "C1COC(=O)O1"  # 1.32
EMC = "CCOC(=O)OC"  # 1.04
EA = "CCOC(C)=O"  # 0.902
PF6 = "F[P-](F)(F)(F)(F)F"
EAf = "CC(=O)OCC(F)(F)F"  # 1.258
FEC = "O=C1OC[C@H](F)O1"  # 1.454
Li = "[Li+]"
fEA = "C(OCC)(=O)C(F)(F)F"  # 1.192
fEAf = "C(OCC(F)(F)F)(=O)C(F)(F)F"  # 1.464

density_dict = {
    "C1COC(=O)O1": 1.32,
    "CCOC(=O)OC": 1.04,
    "CCOC(C)=O": 0.902,
    "CC(=O)OCC(F)(F)F": 1.258,
    "O=C1OC[C@H](F)O1": 1.454,
    "C(OCC)(=O)C(F)(F)F": 1.192,
    "C(OCC(F)(F)F)(=O)C(F)(F)F": 1.464,
}

pf6_charges = np.load('../partial_charges/PF6.npy')
pf6 = Molecule.from_file('../partial_charges/PF6.xyz')
li_charges = np.load('../partial_charges/Li.npy')
li = Molecule.from_file('../partial_charges/Li.xyz')

L_to_A3 = 1e27
A3_to_L = 1e-27

from pymatgen.io.openmm.utils import n_mols_from_volume_ratio, get_box, n_solute_from_molarity

smiles = [EC, EMC]
vol_ratio = [3, 7]
densities = [density_dict[smile] for smile in smiles]
counts = n_mols_from_volume_ratio(600, smiles, vol_ratio, densities)
smile_counts = {smile: count for smile, count in zip(smiles, counts)}
total_density = np.sum(np.array(vol_ratio) * np.array(densities) / np.sum(vol_ratio)) + 0.03
volume = get_box(smile_counts, total_density)[5] ** 3 * A3_to_L
n_solute = n_solute_from_molarity(1.2, volume)

temps = [298, 273, 253, 233]

# gen 2
# for temp in temps:
#     generator = OpenMMSolutionGen(
#         partial_charge_scaling={Li: 0.7, PF6: 0.7},
#         partial_charges=[(pf6, pf6_charges), (li, li_charges)],
#         temperature=temp,
#     )
#     input_set = generator.get_input_set(
#         {EA: 235, FEC: 365, Li: 61, PF6: 61},
#         density=1.154
#     )
#     input_set.write_input(f'real_sims/gen2/{temp}')
#
# # EA 1M
# for temp in temps:
#     generator = OpenMMSolutionGen(
#         partial_charge_scaling={Li: 0.7, PF6: 0.7},
#         partial_charges=[(pf6, pf6_charges), (li, li_charges)],
#         temperature=temp,
#     )
#     input_set = generator.get_input_set(
#         {EA: 522, FEC: 78, Li: 54, PF6: 54},
#         density=0.99
#     )
#     input_set.write_input(f'real_sims/EA/{temp}')

# EA 3M
for temp in temps:
    generator = OpenMMSolutionGen(
        partial_charge_scaling={Li: 0.7, PF6: 0.7},
        partial_charges=[(pf6, pf6_charges), (li, li_charges)],
        temperature=temp,
    )
    input_set = generator.get_input_set(
        {EA: 522, FEC: 78, Li: 165, PF6: 165},
        density=0.99
    )
    input_set.write_input(f'real_sims/EA_3M/{temp}')

# EAf
for temp in temps:
    generator = OpenMMSolutionGen(
        partial_charge_scaling={Li: 0.7, PF6: 0.7},
        partial_charges=[(pf6, pf6_charges), (li, li_charges)],
        temperature=temp,
    )
    input_set = generator.get_input_set(
        {EAf: 512, FEC: 88, Li: 63, PF6: 63},
        density=1.30
    )
    input_set.write_input(f'real_sims/EAf/{temp}')

# fEA
for temp in temps:
    generator = OpenMMSolutionGen(
        partial_charge_scaling={Li: 0.7, PF6: 0.7},
        partial_charges=[(pf6, pf6_charges), (li, li_charges)],
        temperature=temp,
    )
    input_set = generator.get_input_set(
        {fEA: 508, FEC: 92, Li: 66, PF6: 66},
        density=1.24
    )
    input_set.write_input(f'real_sims/fEA/{temp}')

# fEAf
for temp in temps:
    generator = OpenMMSolutionGen(
        partial_charge_scaling={Li: 0.7, PF6: 0.7},
        partial_charges=[(pf6, pf6_charges), (li, li_charges)],
        temperature=temp,
    )
    input_set = generator.get_input_set(
        {fEAf: 498, FEC: 102, Li: 73, PF6: 73},
        density=1.46
    )
    input_set.write_input(f'real_sims/fEAf/{temp}')

