from mpi4py import MPI

from pymatgen.io.openmm.generators import OpenMMSolutionGen
from pymatgen.io.openmm.simulations import equilibrate_pressure, anneal
from pymatgen.core import Molecule
import time

import numpy as np

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

temperature_dict = {0: 298, 1: 273, 2: 253, 3: 233}

temperature = temperature_dict[int(rank)]

generator = OpenMMSolutionGen(
	partial_charge_scaling={Li: 0.8, PF6: 0.8},
	partial_charges=[(pf6, pf6_charges), (li, li_charges)],
	temperature=temperature
	)

input_set = generator.get_input_set(
	{TFEA: 512, FEC: 88, Li: 54, PF6: 54},
	density=1.36
)

run_names = { #change
    0: "TFEA_298",
    1: "TFEA_273",
    2: "TFEA_253",
    3: "TFEA_233",
}

properties = {"DeviceIndex": f"1"}  # change
sim = input_set.get_simulation(platformProperties={"DeviceIndex": str(rank)})

dir = "tfea_runs"  # change
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
print("end time: ", start - time.time())
