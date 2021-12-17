import time
import pathlib
import sys
import os
from mpi4py import MPI

from openmm.app import DCDReporter, StateDataReporter, PDBReporter
from pymatgen.io.openmm.sets import OpenMMSet
from openmm import Platform
from pymatgen.io.openmm.simulations import anneal, equilibrate_pressure


"""
argument 1: the input directory
argument 2: the output directory
"""

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

wd = pathlib.Path.cwd()
input_dir = pathlib.Path(sys.argv[1])
output_dir = pathlib.Path(sys.argv[2])

target_paths = [path for path in input_dir.iterdir()]
assert len(target_paths) == 4, "there must be exactly 4 jobs to run"

target_path = target_paths[rank]

input_set = OpenMMSet.from_directory(target_path)

platform = Platform.getPlatformByName("OpenCL")

sim = input_set.get_simulation(
    platform=platform,
    platformProperties={"DeviceIndex": str(rank)},
)

output_dir.mkdir(parents=True, exist_ok=True)

pdb_reporter = PDBReporter(str(output_dir / "topology.pdb"), 1)
pdb_reporter.report(sim, sim.context.getState(getPositions=True))

sim.reporters.append(DCDReporter(str(output_dir / "trajectory.dcd"), 1000))
sim.reporters.append(
    StateDataReporter(
        str(output_dir / "state.txt"),
        1000,
        step=True,
        potentialEnergy=True,
        temperature=True,
        volume=True,
        density=True,
    )
)

start = time.time()

sim.minimizeEnergy()
equilibrate_pressure(sim, 1000)
anneal(sim, 400, [1000, 1000, 1000])
sim.step(2000)

print("total runtime: ", time.time() - start)
