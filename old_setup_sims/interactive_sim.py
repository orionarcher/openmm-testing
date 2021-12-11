from setup.setup_functions import _smiles_to_simulation

from openmm.unit import *

from openmm.app import PDBReporter
from openmm.app import StateDataReporter
from sys import stdout
from openmm import MonteCarloBarostat, Platform

from openmm import *
from openmm.app import *

import time

EA = "O"
PF6 = "F[P-](F)(F)(F)(F)F"
TFEA = "CCO"
Li = "[Li+]"


properties = {"DeviceIndex": "1"}

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
print(platform.getPropertyNames())


print(platform.getPropertyValue(context, "DeviceIndex"))
start = time.time()
sim.minimizeEnergy()
sim.step(1000)
end = time.time()
print(end - start)
print('hi')
# sim.saveState("state_test")


# 600 EA system
# 50000 steps no energy min
# 1 device OpenCL Perlmutter: 13.53s
# 2 device OpenCL Perlmutter: 17.36
# 4 device OpenCL Perlmutter: 18.49
# 1 node cori haswell: 13 min
# lammps extrapolated: 13.2 min
# my computer 16 threads, extrapolated: 14 min
# my computer OpenCL, extrapolated: 19 min



# 600 EA system
# 1 device OpenCL Perlmutter: 31s
# 2 device OpenCL Perlmutter: 39s
# 4 device OpenCL Perlmutter: 40s
# 1 node cori haswell: 13 min
# lammps extrapolated: 13.2 min
# my computer 16 threads, extrapolated: 14 min
# my computer OpenCL, extrapolated: 19 min
