from setup.setup_functions import _smiles_to_simulation
from openmm import MonteCarloBarostat, Platform

from openmm.unit import *

from sys import stdout
from openmm.app import *
import numpy as np

import time

EA = "CCOC(C)=O"
PF6 = "F[P-](F)(F)(F)(F)F"
TFEA = "CC(=O)OCC(F)(F)F"
Li = "[Li+]"


from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

print("rank: ", rank)

properties = {"DeviceIndex": f"{rank}"}

scaling_dict = {0: 0.7, 1: 0.8, 2: 0.9, 3: 1.0}
charge_scaling = scaling_dict[rank]

print("charge scaling: ", charge_scaling)

sim = _smiles_to_simulation(
    [EA, TFEA, Li, PF6],
    [600, 100, 50, 50],
    50,
    charge_scaling=0.7,
    properties=properties,
)

sim.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True, volume=True, density=True))


context = sim.context
platform = context.getPlatform()
system = context.getSystem()
integrator = context.getIntegrator()


start = time.time()

# minimize Energy
sim.minimizeEnergy()

# volume equilibration
barostat_force_index = system.addForce(MonteCarloBarostat(1*atmosphere, 300*kelvin, 10))
sim.context.reinitialize(preserveState=True)
sim.step(10000)
system.removeForce(barostat_force_index)
sim.context.reinitialize(preserveState=True)

# annealing
for i in np.arange(0, 100, 1):
    integrator.setTemperature(300 * kelvin + 1 * i * kelvin)
    sim.step(100)

sim.step(1000)

for i in np.arange(100, 0, -1):
    integrator.setTemperature(300 * kelvin + 1 * i * kelvin)
    sim.step(100)

sim.step(1000)

end = time.time()
print("time: ", end - start)
print(platform.getPropertyNames())
print(platform.getPropertyValue(context, "DeviceIndex"))






