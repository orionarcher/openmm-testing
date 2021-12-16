from pymatgen.io.openmm.generators import OpenMMSolutionGen

from pymatgen.core import Molecule

import numpy as np

EA = "CCOC(C)=O"
PF6 = "F[P-](F)(F)(F)(F)F"
EAf = "CC(=O)OCC(F)(F)F"
FEC = "O=C1OC[C@H](F)O1"
Li = "[Li+]"
fEA = "C(OCC)(=O)C(F)(F)F"
fEAf = "C(OCC(F)(F)F)(=O)C(F)(F)F"

pf6_charges = np.load('../partial_charges/PF6.npy')
pf6 = Molecule.from_file('../partial_charges/PF6.xyz')
li_charges = np.load('../partial_charges/Li.npy')
li = Molecule.from_file('../partial_charges/Li.xyz')


FEC_mixes = {
    "EA300_FEC300": {EA: 300, FEC: 300},
    "EAf300_FEC300": {EAf: 300, FEC: 300},
    "fEA300_FEC300": {fEA: 300, FEC: 300},
    "fEAf300_FEC300": {fEAf: 300, FEC: 300}
}
EA_mixes = {
    "EA300_EAf300": {EA: 300, EAf: 300},
    "EA300_fEA300": {EA: 300, fEA: 300},
    "EA300_fEAf300": {EA: 300, fEAf: 300},
    "all150": {EA: 150, EAf: 150, fEA: 150, fEAf: 150}
}
other_mixes = {
    "EAf300_fEA300": {EAf: 300, fEA: 300},
    "EAf300_fEAf300": {EAf: 300, fEAf: 300},
    "fEA300_fEAf300": {fEA: 300, fEAf: 300},
    "all200": {EAf: 200, fEA: 200, fEAf: 200}
}

generator = OpenMMSolutionGen(
    partial_charge_scaling={Li: 0.7, PF6: 0.7},
    partial_charges=[(pf6, pf6_charges), (li, li_charges)],
)

for amounts in FEC_mixes:
    all_amounts = {**amounts, Li: 55, PF6: 55}
    input_set = generator.get_input_set(
        all_amounts,
        density=1.4
    )
    input_set.write_input(f'mixed_sims/{name}')


temps = [298, 273, 253, 233]

for temp in temps:
    generator = OpenMMSolutionGen(
        partial_charge_scaling={Li: 0.7, PF6: 0.7},
        partial_charges=[(pf6, pf6_charges), (li, li_charges)],
        temperature=temp,
    )
    input_set = generator.get_input_set(
        {EAf: 512, FEC: 88, Li: 54, PF6: 54},
        density=1.36
    )
    input_set.write_input(f'EAf_sims/{temp}')


for temp in temps:
    generator = OpenMMSolutionGen(
        partial_charge_scaling={Li: 0.7, PF6: 0.7},
        partial_charges=[(pf6, pf6_charges), (li, li_charges)],
        temperature=temp,
    )
    input_set = generator.get_input_set(
        {EA: 522, FEC: 78, Li: 54, PF6: 54},
        density=1.06
    )
    input_set.write_input(f'EA_sims/{temp}')

