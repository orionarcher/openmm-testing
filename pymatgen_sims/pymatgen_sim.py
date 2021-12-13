EA = "CCOC(C)=O"
PF6 = "F[P-](F)(F)(F)(F)F"
TFEA = "CC(=O)OCC(F)(F)F"
FEC = "O=C1OC[C@H](F)O1"
Li = "[Li+]"

import os

from mpi4py import MPI

from pymatgen.io.openmm.generators import OpenMMSolutionGen
from pymatgen.io.openmm.simulations import equilibrate_pressure, anneal


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

print("rank: ", rank)

temperature = 1
print("charge scaling: ", charge_scaling)
properties = {"DeviceIndex": f"{rank}"}

platform = Platform.getPlatformByName('OpenCL')
platform.setDefaultPropertyValue("DeviceIndex", str(rank))


generator = OpenMMSolutionGen(
	charge_scaling={Li: 0.8, PF6: 0.8}
	)

input_set = generator.get_input_set(
	{EA: 522, FEC: 78, Li: 54, PF6: 54},
	density = 1.06
)

input_set.get_simulation(platform, {"DeviceIndex": str(rank)})