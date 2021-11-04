from setup.setup_functions import _smiles_to_simulation

from openmm.unit import *

from openmm.app import PDBReporter
from openmm.app import StateDataReporter
from sys import stdout
from openmm import MonteCarloBarostat, Platform

from openmm import *
from openmm.app import *

EA = "CCOC(C)=O"
PF6 = "F[P-](F)(F)(F)(F)F"
TFEA = "CC(=O)OCC(F)(F)F"
Li = "[Li+]"


properties = {"DeviceIndex": "0,1,2,3"}

sim = _smiles_to_simulation(
    [EA],
    [600],
    80,
    platformProperties=properties,
)

sim.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True, volume=True, density=True))


context = sim.context
platform = context.getPlatform()
print(platform.getPropertyNames())
