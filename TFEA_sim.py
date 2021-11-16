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
FEC = "C1C(OC(=O)O1)F"
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
    [TFEA, FEC, Li, PF6],
    [512, 88, 62, 62],
    47.8,
    charge_scaling=charge_scaling,
    properties=properties,
)


context = sim.context
platform = context.getPlatform()
system = context.getSystem()
integrator = context.getIntegrator()
positions = context.getState(getPositions=True).getPositions(asNumpy=True)

sim.reporters.append(StateDataReporter(f"TFEA_output/state_{charge_scaling}.txt", 1000, step=True, potentialEnergy=True, temperature=True, volume=True, density=True))
sim.reporters.append(DCDReporter(f"TFEA_output/traj_{charge_scaling}.dcd", 1000))

pdb_reporter = PDBReporter(f"TFEA_output/top_{charge_scaling}.pdb", 1)
pdb_reporter.report(sim, context.getState(getPositions=True))

start = time.time()

# minimize Energy
sim.minimizeEnergy()

# volume equilibration
barostat_force_index = system.addForce(MonteCarloBarostat(1*atmosphere, 300*kelvin, 10))
sim.context.reinitialize(preserveState=True)
sim.step(200000)
system.removeForce(barostat_force_index)
sim.context.reinitialize(preserveState=True)

# annealing
for i in np.arange(0, 100, 1):
    integrator.setTemperature(300 * kelvin + 1 * i * kelvin)
    sim.step(10000)

sim.step(1000000)

for i in np.arange(100, 0, -1):
    integrator.setTemperature(300 * kelvin + 1 * i * kelvin)
    sim.step(1)

sim.step(5000000)

end = time.time()
print("time: ", end - start)
print(platform.getPropertyNames())
print(platform.getPropertyValue(context, "DeviceIndex"))






