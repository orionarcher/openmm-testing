from mpi4py import MPI

from pymatgen.io.openmm.generators import OpenMMSolutionGen
from pymatgen.io.openmm.simulations import equilibrate_pressure, anneal

EA = "CCOC(C)=O"
PF6 = "F[P-](F)(F)(F)(F)F"
TFEA = "CC(=O)OCC(F)(F)F"
FEC = "O=C1OC[C@H](F)O1"
Li = "[Li+]"

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

print("rank: ", rank)

temperature_dict = {0: 298, 1: 273, 2: 253, 3: 233}

temperature = temperature_dict[int(rank)]

generator = OpenMMSolutionGen(
	charge_scaling={Li: 0.8, PF6: 0.8},
	temperature=temperature
	)

input_set = generator.get_input_set(
	{TFEA: 512, FEC: 88, Li: 54, PF6: 54},
	density = 1.36
)

properties = {"DeviceIndex": f"{rank}"}
platform = Platform.getPlatformByName('OpenCL')
input_set.get_simulation(platform, {"DeviceIndex": str(rank)})
