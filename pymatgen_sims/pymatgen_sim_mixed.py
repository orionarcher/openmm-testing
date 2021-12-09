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

EA_TFEA_FEC = {
	0: {"EA": 300, "TFEA": 0, "FEC": 300},
	1: {"EA": 300, "TFEA": 300, "FEC": 0},
	2: {"EA": 0, "TFEA": 300, "FEC": 300},
	3: {"EA": 250, "TFEA": 250, "FEC": 100},
}
amounts = EA_TFEA_FEC[int(rank)]

generator = OpenMMSolutionGen(
	partial_charge_scaling={Li: 0.8, PF6: 0.8}
	)

input_set = generator.get_input_set(
	{EA: amounts["EA"], TFEA: amounts["TFEA"], FEC: amounts["FEC"], Li: 55, PF6: 55},
	density = 1.36
)

properties = {"DeviceIndex": f"{rank}"}
platform = Platform.getPlatformByName('OpenCL')
input_set.get_simulation(platform, {"DeviceIndex": str(rank)})