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

all_mixes = {
    "FEC": FEC_mixes,
    "EA": EA_mixes,
    "other": other_mixes
}

generator = OpenMMSolutionGen(
    partial_charge_scaling={Li: 0.7, PF6: 0.7},
    partial_charges=[(pf6, pf6_charges), (li, li_charges)],
)
for mix_set_name, mix_set in all_mixes.items():
    for run_name, amounts in mix_set.items():
        all_amounts = {**amounts, Li: 55, PF6: 55}
        input_set = generator.get_input_set(
            all_amounts,
            density=1.4
        )
        input_set.write_input(f'mixed_sims/{mix_set_name}/{run_name}')
