import os

from simulation_array import input_sets
from mpi4py import MPI
from openmm import Platform
from openmm.app import StateDataReporter, DCDReporter, PDBReporter

from pymatgen.io.openmm.simulations import anneal, equilibrate_pressure

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

run_names = [
    "300_EA_300_FEC",
    "300_EA_300_TFEA",
    "300_TFEA_300_FEC",
    "250_EA_250_TFEA_100_FEC",
]

print("rank: ", rank)

if rank < 4:
    name = run_names[rank]
    input_set = input_sets[name]
    properties = {"DeviceIndex": f"{rank}"}
    platform = Platform.getPlatformByName("OpenCL")
    sim = input_set.get_simulation(
        platformProperties={"DeviceIndex": str(rank)}
    )
    os.mkdir('mixed_runs')
    sim.reporters.append(StateDataReporter(f"mixed_runs/{name}/state.txt", 1000, step=True, potentialEnergy=True, temperature=True, volume=True, density=True))
    sim.reporters.append(DCDReporter(f"mixed_runs/{name}/trajectory.dcd", 1000))

    pdb_reporter = PDBReporter(f"mixed_runs/{name}/topology.pdb", 1)
    pdb_reporter.report(sim, sim.context.getState(getPositions=True))
    sim.minimizeEnergy()
    equilibrate_pressure(sim, 2000000)
    anneal(sim, 400, [1000000, 1000000, 1000000])
    sim.step(5000000)
