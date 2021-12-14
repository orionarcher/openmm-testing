from mpi4py import MPI

from pymatgen.io.openmm.generators import OpenMMSolutionGen
from pymatgen.io.openmm.simulations import equilibrate_pressure, anneal
from pymatgen.core import Molecule

import numpy as np

import time
import os

from openmm.app import StateDataReporter, PDBReporter, DCDReporter

EA = "CCOC(C)=O"
PF6 = "F[P-](F)(F)(F)(F)F"
TFEA = "CC(=O)OCC(F)(F)F"
FEC = "O=C1OC[C@H](F)O1"
Li = "[Li+]"

pf6_charges = np.load('../partial_charges/PF6.npy')
pf6 = Molecule.from_file('../partial_charges/PF6.xyz')
li_charges = np.load('../partial_charges/Li.npy')
li = Molecule.from_file('../partial_charges/Li.xyz')

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

print("rank: ", rank)

EA_TFEA_FEC = {
	0: {"EA": 300, "TFEA": 0, "FEC": 300},
	1: {"EA": 300, "TFEA": 300, "FEC": 0},
	2: {"EA": 0, "TFEA": 300, "FEC": 300},
	3: {"EA": 250, "TFEA": 250, "FEC": 100},
}
amounts = EA_TFEA_FEC[int(rank)]

generator = OpenMMSolutionGen(
	partial_charge_scaling={Li: 0.8, PF6: 0.8},
	partial_charges=[(pf6, pf6_charges), (li, li_charges)],
)

input_set = generator.get_input_set(
	{EA: amounts["EA"], TFEA: amounts["TFEA"], FEC: amounts["FEC"], Li: 55, PF6: 55},
	density=1.36
)

run_names = { #change
    0: "300_EA_300_FEC",
    1: "300_EA_300_TFEA",
    2: "300_TFEA_300_FEC",
    3: "250_EA_250_TFEA_100_FEC",
}

properties = {"DeviceIndex": f"{rank}"}  # change
sim = input_set.get_simulation(platformProperties={"DeviceIndex": str(rank)})

dir = "mixed_runs"  # change
name = run_names[rank]
os.makedirs(f"{dir}/{name}", exist_ok=True)
sim.reporters.append(
    StateDataReporter(f"{dir}/{name}/state.txt", 1000, step=True, potentialEnergy=True, temperature=True,
                      volume=True, density=True))
sim.reporters.append(DCDReporter(f"{dir}/{name}/trajectory.dcd", 1000))

pdb_reporter = PDBReporter(f"{dir}/{name}/topology.pdb", 1)
pdb_reporter.report(sim, sim.context.getState(getPositions=True))
sim.minimizeEnergy()
start = time.time()
equilibrate_pressure(sim, 1000)
anneal(sim, 400, [1000, 1000, 1000])
sim.step(1000)
print("end time: ", time.time() - start)
