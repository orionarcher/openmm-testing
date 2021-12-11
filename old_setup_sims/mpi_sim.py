from setup.setup_functions import _smiles_to_simulation

from sys import stdout
from openmm.app import *

import time

EA = "CCOC(C)=O"
PF6 = "F[P-](F)(F)(F)(F)F"
TFEA = "CC(=O)OCC(F)(F)F"
Li = "[Li+]"


from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

print(rank)

properties = {"DeviceIndex": f"{rank}"}

sim = _smiles_to_simulation(
    [EA, TFEA, Li, PF6],
    [600, 100, 50, 50],
    50,
    properties=properties,
)

sim.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True, volume=True, density=True))


context = sim.context
platform = context.getPlatform()

sim.minimizeEnergy()
start = time.time()
sim.step(1000)
end = time.time()

print("time:", end - start)
print("name:", platform.getName())
print("device index:", platform.getPropertyValue(context, "DeviceIndex"))
