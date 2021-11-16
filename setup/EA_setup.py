from openbabel import pybel

EA = 'CCOC(C)=O'
PF6 = 'F[P-](F)(F)(F)(F)F'
TFEA = 'CC(=O)OCC(F)(F)F'
Li = '[Li+]'

# get the correct force field
import openff.toolkit as tk
from openff.toolkit import topology
from openff.toolkit.typing.engines.smirnoff import ForceField
import openmm


from openbabel import pybel
from openmm.app import PDBFile
from openbabel import openbabel
from pymatgen.io.babel import BabelMolAdaptor
import parmed as pmd
import tempfile

from pymatgen.io.lammps.utils import PackmolRunner
from pymatgen.io.xyz import XYZ


def smile_to_mol(smile):
    mol = pybel.readstring("smi", smile)
    mol.addh()
    mol.make3D()
    adaptor = BabelMolAdaptor(mol.OBMol)
    return adaptor.pymatgen_mol


def smile_to_openmm_system(smile):
    molecule = tk.topology.Molecule.from_smiles(smile)
    topology = tk.topology.Topology.from_molecules(molecule)
    forcefield = tk.typing.engines.smirnoff.ForceField('openff_unconstrained-2.0.0.offxml')
    system = forcefield.create_openmm_system(topology)
    return system


def smile_to_pbd_file(smile):
    mol = pybel.readstring("smi", smile)
    mol.addh()
    mol.make3D()
    with tempfile.NamedTemporaryFile() as f:
        mol.write(format='pdb', filename=f.name, overwrite=True)
        pdbfile = PDBFile(f.name)
    return pdbfile


def get_coordinates(mols, numbers, box_list):
    param_list = [{"number": number, "inside box": box_list} for number in numbers]
    controls = {"maxit": 20, "nloop": 600, 'seed': 20}
    with tempfile.NamedTemporaryFile() as f:
        packmol_runner = PackmolRunner(mols, param_list, control_params=controls, auto_box=False, output_file=f.name)
        packed_box = packmol_runner.run()
    coordinates = XYZ(packed_box).as_dataframe()
    return coordinates


def make_boxes(x, y=None, z=None):
    if not (y and z):
        y, z = x, x
    packmol_box = [0, 0, 0, x, y, z]
    openmm_box = [x, y, z, 90, 90, 90]
    return packmol_box, openmm_box

target_smile = EA
EA_packbox, EA_openbox = make_boxes(40.8)
TFEA_packbox, TFEA_openbox = make_boxes(47.8)


from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout


mol = smile_to_mol(target_smile)
openmm_system = smile_to_openmm_system(target_smile)
pbd = smile_to_pbd_file(target_smile)
EA_sys = pmd.openmm.load_topology(pbd.topology, openmm_system)

coordinates = get_coordinates([mol], [600], EA_packbox)
coord_values = coordinates.loc[:, 'x':'z'].values
coord_values -= coord_values.min(axis=0)
EA_600 = EA_sys * 600
EA_600.coordinates = coord_values
EA_600.box = EA_openbox
EA_system = EA_600.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer)
EA_600.save('EA_600.top', overwrite=True)

integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.001*picoseconds)
EA_simulation = Simulation(EA_600.topology, EA_system, integrator)
EA_simulation.context.setPositions(coord_values)
EA_simulation.saveState('../simulations/EA_simulation.xml')
positions = EA_simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(EA_simulation.topology, positions, open('../scratch/output.pdb', 'w'))

print('hi')

