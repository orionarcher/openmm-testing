from mpi4py import MPI

from pymatgen.io.openmm.generators import OpenMMSolutionGen
from pymatgen.io.openmm.simulations import equilibrate_pressure, anneal
from pymatgen.core import Molecule

import numpy as np

import time
import os

from openmm.app import StateDataReporter, PDBReporter, DCDReporter
import openmm
from openmm import Platform
from old_setup_sims.old_setup_functions import _smiles_to_simulation


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

sim = _smiles_to_simulation(
    [TFEA, FEC, Li, PF6],
    [512, 88, 62, 62],
    47.8,
    charge_scaling=0.7,
    properties={"DeviceIndex": f"{1}"},
    temperature=298
)

print("made first sim")

temperature_dict = {0: 298, 1: 273, 2: 253, 3: 233}

temperature = temperature_dict[int(rank)]

generator = OpenMMSolutionGen(
    partial_charge_scaling={Li: 0.8, PF6: 0.8},
    partial_charges=[(pf6, pf6_charges), (li, li_charges)],
    temperature=temperature
)

input_set = generator.get_input_set(
    {EA: 522, FEC: 78, Li: 54, PF6: 54},
    density=1.06
)

run_names = {  # change
    0: "EA_298",
    1: "EA_273",
    2: "EA_253",
    3: "EA_233",
}


# platform = Platform.getPlatformByName("OpenCL")
properties = {"DeviceIndex": f"{rank}"}  # change
platform = Platform.getPlatformByName("OpenCL")


system = input_set.inputs['system.xml'].get_system()
topology = input_set.inputs['topology.pdb'].get_topology()
integrator = input_set.inputs['integrator.xml'].get_integrator()
state = input_set.inputs['state.xml'].get_state()

sim = openmm.app.Simulation(
    topology,
    system,
    integrator,
    platform=platform,
    platformProperties=properties
)

sim.context.setState(state)

# sim = input_set.get_simulation(
#     platform=platform,
#     platformProperties=properties,
# )

dir = "ea_runs"  # change
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
mid = time.time()
sim.step(2000)

platform = sim.context.getPlatform()
print("rank: ", rank)
print("platform: ", platform.getName())
print("mid time: ", mid - start)
print("end time: ", time.time() - start)
print("device index: ", platform.getPropertyValue(sim.context, "DeviceIndex"))
