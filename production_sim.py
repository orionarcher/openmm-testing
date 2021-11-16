from setup.setup_functions import _smiles_to_simulation
from openmm import MonteCarloBarostat, Platform

from openmm.unit import *

from sys import stdout
from openmm.app import *
import numpy as np

import time

EA = "CCOC(C)=O"
PF6 = "F[P-](F)(F)(F)(F)F"
FEC = "O=C1OC[C@H](F)O1"
TFEA = "CC(=O)OCC(F)(F)F"
Li = "[Li+]"


from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

print("rank: ", rank)

properties = {"DeviceIndex": f"{1}"}

scaling_dict = {0: 0.7, 1: 0.8, 2: 0.9, 3: 1.0}
charge_scaling = scaling_dict[rank]

print("charge scaling: ", charge_scaling)

sim = _smiles_to_simulation(
    [EA, TFEA, FEC, Li, PF6],
    [500, 100, 50, 50, 50],
    50,
    charge_scaling=0.7,
    properties=properties,
)


context = sim.context
platform = context.getPlatform()
system = context.getSystem()
integrator = context.getIntegrator()
positions = context.getState(getPositions=True).getPositions(asNumpy=True)

sim.reporters.append(StateDataReporter("output/state.txt", 1000, step=True, potentialEnergy=True, temperature=True, volume=True, density=True))
sim.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, volume=True, density=True))
sim.reporters.append(DCDReporter("output/test.dcd", 1000))


pdb_reporter = PDBReporter("output/test.pdb", 1000)
pdb_reporter.report(sim, context.getState(getPositions=True))

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
    sim.step(10)

sim.step(1000)

for i in np.arange(100, 0, -1):
    integrator.setTemperature(300 * kelvin + 1 * i * kelvin)
    sim.step(10)

sim.step(100)

end = time.time()
print("time: ", end - start)
print(platform.getPropertyNames())
print(platform.getPropertyValue(context, "DeviceIndex"))






