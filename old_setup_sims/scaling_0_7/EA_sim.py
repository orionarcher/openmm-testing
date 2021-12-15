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

temperature_dict = {
    0: 298,
    1: 273,
    2: 253,
    3: 233,
}
temperature = temperature_dict[rank]

print("temp: ", temperature)

sim = _smiles_to_simulation(
    [EA, FEC, Li, PF6],
    [522, 78, 54, 54],
    40.8,
    charge_scaling=0.8,
    properties=properties,
    temperature=temperature,
)


context = sim.context
platform = context.getPlatform()
system = context.getSystem()
integrator = context.getIntegrator()
positions = context.getState(getPositions=True).getPositions(asNumpy=True)

temp_string = str(temperature)

output_dir = os.mkdir(f"EA_output_2_{temp_string}")

sim.reporters.append(StateDataReporter(f"EA_output_2_{temp_string}/state.txt", 1000, step=True, potentialEnergy=True, temperature=True, volume=True, density=True))
sim.reporters.append(DCDReporter(f"EA_output_2_{temp_string}/trajectory.dcd", 1000))

pdb_reporter = PDBReporter(f"EA_output_2_{temp_string}/topology.pdb", 1)
pdb_reporter.report(sim, context.getState(getPositions=True))

start = time.time()

# minimize Energy
sim.minimizeEnergy()

# volume equilibration
barostat_force_index = system.addForce(MonteCarloBarostat(1*atmosphere, 300*kelvin, 10))
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






