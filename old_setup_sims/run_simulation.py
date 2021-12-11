from setup.setup_functions import _smiles_to_simulation

EA = "CCOC(C)=O"
PF6 = "F[P-](F)(F)(F)(F)F"
TFEA = "CC(=O)OCC(F)(F)F"
Li = "[Li+]"

sim = _smiles_to_simulation(
    [EA],
    [600],
    40,
)

from openmm.unit import *

from openmm.app import PDBReporter
from openmm.app import StateDataReporter
from sys import stdout
from openmm import MonteCarloBarostat, Platform


context = sim.context
platform = Platform.getPlatformByName("OpenCL")
print(platform.getPropertyValue(context, "DeviceIndex"))

# energy minimization
sim.minimizeEnergy()

# add reporters
# sim.reporters.append(PDBReporter('output/output.pdb', 100))
sim.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True, volume=True, density=True))

# add barostat and run
system = sim.context.getSystem()
integrator = sim.context.getIntegrator()
# barostat_force_index = system.getNumForces()
barostat_force_index = system.addForce(MonteCarloBarostat(1*atmosphere, 300*kelvin, 10))
print(system.usesPeriodicBoundaryConditions())
sim.context.reinitialize(preserveState=True)
sim.step(10000)
system.removeForce(barostat_force_index)
sim.context.reinitialize(preserveState=True)
sim.step(1000)
for i in range(10):
    integrator.setTemperature(300 * kelvin + 10 * i * kelvin)
    sim.step(100)

sim.step(1000)

for i in range(1, 11):
    integrator.setTemperature(400 * kelvin - 10 * i * kelvin)
    sim.step(100)

sim.step(1000)
