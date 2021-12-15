from pymatgen.io.openmm.sets import OpenMMSet
from openmm import Platform

sim1 = OpenMMSet.from_directory('good_input_set')
sim2 = OpenMMSet.from_directory('bad_input_set')

opencl = Platform.getPlatformByName("OpenCL")
cpu = Platform.getPlatformByName("CPU")

# sim1.get_simulation()
# sim2.get_simulation()
#
# sim1.get_simulation(cpu)
# sim2.get_simulation(cpu)
#
# sim1.get_simulation(opencl)
sim2.get_simulation(opencl)
