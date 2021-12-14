from old_setup_functions import _smiles_to_simulation
from openmm import MonteCarloBarostat, Platform

from openmm.unit import *

from sys import stdout
from openmm.app import *
import numpy as np

import time

EA = "CCOC(C)=O"
PF6 = "F[P-](F)(F)(F)(F)F"
TFEA = "CC(=O)OCC(F)(F)F"
FEC = "O=C1OC[C@H](F)O1"
Li = "[Li+]"

import os

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

print("rank: ", rank)

properties = {"DeviceIndex": f"{rank}"}

EA_TFEA_FEC = {
    0: {"EA": 300, "TFEA": 1, "FEC": 300},
    1: {"EA": 300, "TFEA": 300, "FEC": 1},
    2: {"EA": 1, "TFEA": 300, "FEC": 300},
    3: {"EA": 250, "TFEA": 250, "FEC": 100},
}
amounts = EA_TFEA_FEC[rank]

run_names = {  # change
    0: "300_EA_300_FEC",
    1: "300_EA_300_TFEA",
    2: "300_TFEA_300_FEC",
    3: "250_EA_250_TFEA_100_FEC",
}
name = run_names[rank]

sim = _smiles_to_simulation(
    [EA, TFEA, FEC, Li, PF6],
    [amounts["EA"], amounts["TFEA"], amounts["FEC"], 55, 55],
    47.8,
    charge_scaling=0.8,
    properties=properties,
)

context = sim.context
platform = context.getPlatform()
system = context.getSystem()
integrator = context.getIntegrator()
positions = context.getState(getPositions=True).getPositions(asNumpy=True)


output_dir = os.mkdir(f"TFEA_output_2_{name}")

sim.reporters.append(
    StateDataReporter(f"TFEA_output_2_{name}/state.txt", 1000, step=True, potentialEnergy=True, temperature=True,
                      volume=True, density=True))
sim.reporters.append(DCDReporter(f"TFEA_output_2_{name}/trajectory.dcd", 1000))

pdb_reporter = PDBReporter(f"TFEA_output_2_{name}/topology.pdb", 1)
pdb_reporter.report(sim, context.getState(getPositions=True))

start = time.time()

# minimize Energy
sim.minimizeEnergy()

# volume equilibration
barostat_force_index = system.addForce(MonteCarloBarostat(1 * atmosphere, 300 * kelvin, 10))
sim.context.reinitialize(preserveState=True)
sim.step(1000000)
system.removeForce(barostat_force_index)
sim.context.reinitialize(preserveState=True)

# annealing
for i in np.arange(0, 100, 1):
    integrator.setTemperature(300 * kelvin + 1 * i * kelvin)
    sim.step(10000)

sim.step(1000000)

for i in np.arange(100, 0, -1):
    integrator.setTemperature(300 * kelvin + 1 * i * kelvin)
    sim.step(10000)

sim.step(5000000)

end = time.time()
print("time: ", end - start)
print(platform.getPropertyNames())
print(platform.getPropertyValue(context, "DeviceIndex"))
