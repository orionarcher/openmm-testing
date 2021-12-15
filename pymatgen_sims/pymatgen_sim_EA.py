# from mpi4py import MPI

from pymatgen.io.openmm.generators import OpenMMSolutionGen
from pymatgen.io.openmm.inputs import TopologyInput, StateInput, SystemInput, IntegratorInput
from pymatgen.io.openmm.sets import OpenMMSet
# from pymatgen.io.openmm.simulations import equilibrate_pressure, anneal
from pymatgen.core import Molecule
#
import numpy as np
#
# import time
# import os
#
# from openmm.app import StateDataReporter, PDBReporter, DCDReporter
import openmm
from openmm import Platform, LangevinMiddleIntegrator
from openmm.app import Simulation
from openmm.unit import picosecond, kelvin
from old_setup_functions import _smiles_to_simulation


EA = "CCOC(C)=O"
PF6 = "F[P-](F)(F)(F)(F)F"
TFEA = "CC(=O)OCC(F)(F)F"
FEC = "O=C1OC[C@H](F)O1"
Li = "[Li+]"

pf6_charges = np.load('../partial_charges/PF6.npy')
pf6 = Molecule.from_file('../partial_charges/PF6.xyz')
li_charges = np.load('../partial_charges/Li.npy')
li = Molecule.from_file('../partial_charges/Li.xyz')

# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()

good_sim = _smiles_to_simulation(
    [TFEA, FEC, Li, PF6],
    [512, 88, 62, 62],
    47.8,
    charge_scaling=0.7,
    properties={"DeviceIndex": f"{0}"},
    temperature=298
)

integrator = LangevinMiddleIntegrator(
    300 * kelvin, 1 / picosecond, 0.001 * picosecond
)
good_topology = good_sim.topology
good_system = good_sim.context.getSystem()
good_state = good_sim.context.getState()

good_topology_input = TopologyInput(good_topology)
good_system_input = SystemInput(good_system)
good_state_input = StateInput(good_state)
good_integrator_input = IntegratorInput(integrator)


generator = OpenMMSolutionGen(
    partial_charge_scaling={Li: 0.8, PF6: 0.8},
    partial_charges=[(pf6, pf6_charges), (li, li_charges)],
    temperature=298
)

input_set = generator.get_input_set(
    {EA: 522, FEC: 78, Li: 54, PF6: 54},
    density=1.06
)

properties = {"DeviceIndex": f"{0}"}
opencl = Platform.getPlatformByName("OpenCL")
cpu = Platform.getPlatformByName("CPU")


bad_system = input_set.inputs['system.xml'].get_system()
bad_topology = input_set.inputs['topology.pdb'].get_topology()
bad_integrator = LangevinMiddleIntegrator(
    300 * kelvin, 1 / picosecond, 0.001 * picosecond
)
bad_state = input_set.inputs['state.xml'].get_state()

input_set.inputs['topology.pdb'] = good_topology_input
input_set.inputs['system.xml'] = good_system_input
input_set.inputs['state.xml'] = good_state_input
input_set.inputs['integrator.xml'] = good_integrator_input

sim = Simulation(
    good_topology,
    good_system,
    integrator
)

# sim = input_set.get_simulation(
#     platform=platform,
#     platformProperties=properties,
# )

input_set_2 = OpenMMSet(
    inputs={
        'topology.pdb': good_topology_input,
        'system.xml': good_system_input,
        'state.xml': good_state_input,
        'integrator.xml': good_integrator_input,
    },
    topology_file='topology.pdb',
    system_file='system.xml',
    integrator_file='integrator.xml',
    state_file='state.xml',
)

input_set_2.get_simulation()

input_set_2.get_simulation(platform=cpu)
input_set_2.get_simulation(platform=opencl)

input_set_2.get_simulation(platform=opencl, platformProperties=properties)

# bad_sim = openmm.app.Simulation(
#     topology,
#     system,
#     integrator,
#     platform=platform,
#     # platformProperties=properties
# )
print('complete')

# sim.context.setState(state)

# sim = input_set.get_simulation(
#     platform=platform,
#     platformProperties=properties,
# )
#
# dir = "ea_runs"  # change
# name = run_names[rank]
# os.makedirs(f"{dir}/{name}", exist_ok=True)
# sim.reporters.append(
#     StateDataReporter(f"{dir}/{name}/state.txt", 1000, step=True, potentialEnergy=True, temperature=True,
#                       volume=True, density=True))
# sim.reporters.append(DCDReporter(f"{dir}/{name}/trajectory.dcd", 1000))
#
# pdb_reporter = PDBReporter(f"{dir}/{name}/topology.pdb", 1)
# pdb_reporter.report(sim, sim.context.getState(getPositions=True))
# sim.minimizeEnergy()
# start = time.time()
# equilibrate_pressure(sim, 1000)
# anneal(sim, 400, [1000, 1000, 1000])
# mid = time.time()
# sim.step(2000)
#
# platform = sim.context.getPlatform()
# print("rank: ", rank)
# print("platform: ", platform.getName())
# print("mid time: ", mid - start)
# print("end time: ", time.time() - start)
# print("device index: ", platform.getPropertyValue(sim.context, "DeviceIndex"))
