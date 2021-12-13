from mpi4py import MPI

from pymatgen.io.openmm.generators import OpenMMSolutionGen
from pymatgen.core import Molecule

import numpy as np
import os

EA = "CCOC(C)=O"
PF6 = "F[P-](F)(F)(F)(F)F"
TFEA = "CC(=O)OCC(F)(F)F"
FEC = "O=C1OC[C@H](F)O1"
Li = "[Li+]"

amounts_dict = {
    "300_EA_300_FEC": {"EA": 300, "TFEA": 0, "FEC": 300},
    "300_EA_300_TFEA": {"EA": 300, "TFEA": 300, "FEC": 0},
    "300_TFEA_300_FEC": {"EA": 0, "TFEA": 300, "FEC": 300},
    "250_EA_250_TFEA_100_FEC": {"EA": 250, "TFEA": 250, "FEC": 100},
}

print(os.getcwd())
pf6_charges = np.load('../partial_charges/PF6.npy')
pf6 = Molecule.from_file('../partial_charges/PF6.xyz')
li_charges = np.load('../partial_charges/Li.npy')
li = Molecule.from_file('../partial_charges/Li.xyz')


input_sets = {}

for name, amounts in amounts_dict.items():
    generator = OpenMMSolutionGen(
        partial_charge_scaling={Li: 0.8, PF6: 0.8},
        partial_charges=[(pf6, pf6_charges), (li, li_charges)],
    )
    input_set = generator.get_input_set(
        {EA: amounts["EA"], TFEA: amounts["TFEA"], FEC: amounts["FEC"], Li: 55, PF6: 55},
        density=1.36
    )
    input_sets[name] = input_set
    print("hi")

tfea_temp_dict = {"TFEA_298": 298, "TFEA_273": 273, "TFEA_253": 253, "TFEA_233": 233}

for name, temp in tfea_temp_dict.items():
    generator = OpenMMSolutionGen(
        partial_charge_scaling={Li: 0.8, PF6: 0.8},
        partial_charges=[(pf6, pf6_charges), (li, li_charges)],
        temperature=temp,
    )
    input_set = generator.get_input_set(
        {TFEA: 512, FEC: 88, Li: 54, PF6: 54},
        density=1.36
    )
    input_sets[name] = input_set

ea_temp_dict = {"EA_298": 298, "EA_273": 273, "EA_253": 253, "EA_233": 233}

for name, temp in ea_temp_dict.items():
    generator = OpenMMSolutionGen(
        partial_charge_scaling={Li: 0.8, PF6: 0.8},
        partial_charges=[(pf6, pf6_charges), (li, li_charges)],
        temperature=temp,
    )
    input_set = generator.get_input_set(
        {EA: 522, FEC: 78, Li: 54, PF6: 54},
        density=1.06
    )
    input_sets[name] = input_set

print(len(input_sets))
