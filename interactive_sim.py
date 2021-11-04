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

sim = _smiles_to_simulation(
    [EA],
    [600],
    80,
)

context = sim.context
platform = Platform.getPlatformByName("OpenCl")
print(platform.getPropertyValue(context, "DeviceIndex"))
