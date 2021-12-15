from openbabel import pybel
import functools
from pymatgen.io.babel import BabelMolAdaptor
import parmed
import tempfile
from pymatgen.io.lammps.utils import PackmolRunner
from pymatgen.io.xyz import XYZ
from openmm.unit import *
from openmm.app import ForceField, Simulation, PME
from openmm import LangevinMiddleIntegrator
from openff.toolkit.typing.engines import smirnoff
import openff.toolkit as tk
import openff
from openff.toolkit import topology
from openmm import Platform
import numpy as np

def smile_to_structure(smile):
    """
    Converts a SMILE to a Parmed structure.

    Parameters
    ----------
    smile: a SMILE.

    Returns
    -------
    Parmed.Structure
    """
    mol = pybel.readstring("smi", smile)
    mol.addh()
    mol.make3D()
    with tempfile.NamedTemporaryFile() as f:
        mol.write(format="pdb", filename=f.name, overwrite=True)
        structure = parmed.load_file(f.name)
    return structure


def _smiles_to_topology(smiles, counts):
    """
    Converts a set of SMILEs and counts to an OpenMM Topology object.

    Parameters
    ----------
    smiles: a list of SMILEs.
    counts: counts for each SMILE.

    Returns
    -------

    """
    assert len(smiles) == len(counts), "smiles and counts must be the same length"
    structures = [smile_to_structure(smile) for smile in smiles]
    combined_structs = parmed.Structure()
    for struct, count in zip(structures, counts):
        combined_structs += struct * count
    return combined_structs.topology


def smile_to_mol(smile):
    """
    Converts a SMILE to a Pymatgen Molecule.

    Parameters
    ----------
    smile: a SMILE.

    Returns
    -------
    Pymatgen.Molecule
    """
    mol = pybel.readstring("smi", smile)
    mol.addh()
    mol.make3D()
    adaptor = BabelMolAdaptor(mol.OBMol)
    return adaptor.pymatgen_mol


def make_box_list(box_size, box_type):
    xyz_to_box_funcs = {
        "packmol": lambda xyz: [0, 0, 0] + xyz,
        "openmm": lambda xyz: xyz + [90, 90, 90],
    }
    assert (
        box_type in xyz_to_box_funcs.keys()
    ), f"allowed box types: {' '.join(xyz_to_box_funcs.keys())}"
    if type(box_size) == int or float:
        xyz = [box_size, box_size, box_size]
        return xyz_to_box_funcs[box_type](xyz)
    elif len(box_size) == 3:
        xyz = box_size
        return xyz_to_box_funcs[box_type](xyz)
    else:
        TypeError(f"box_size must be a number (cube) or a list (rectangular prism), it is f{box_size}")


def _smiles_to_coordinates(smiles, counts, box_size):
    box_size = make_box_list(box_size, 'packmol')
    param_list = [{"number": count, "inside box": box_size} for count in counts]
    controls = {"maxit": 20, "nloop": 600, "seed": 20}
    mols = [smile_to_mol(smile) for smile in smiles]
    with tempfile.NamedTemporaryFile() as f:
        packmol_runner = PackmolRunner(
            mols,
            param_list,
            control_params=controls,
            auto_box=False,
            output_file=f.name,
        )
        packed_box = packmol_runner.run()
    coordinates = XYZ(packed_box).as_dataframe()
    raw_coordinates = coordinates.loc[:, "x":"z"].values
    return raw_coordinates

def _smiles_to_cube_size(smiles, counts, density):
    cm3_to_A3 = 1e24
    NA = 6.02214e23
    mols = [smile_to_mol(smile) for smile in smiles]
    mol_mw = np.array([mol.structure.composition.weight for mol in mols])
    counts = np.array(counts)
    total_weight = sum(mol_mw * counts)
    box_volume = total_weight * cm3_to_A3 / (NA * density)
    side_length = box_volume ** (1 / 3)
    return side_length


def _smiles_to_system(smiles, counts, density=1.5):
    topology = _smiles_to_topology(smiles, counts)
    box_size = _smiles_to_cube_size(smiles, counts, density)
    coordinates = _smiles_to_coordinates(smiles, counts, box_size)
    openff_mols = [
        openff.toolkit.topology.Molecule.from_smiles(smile) for smile in smiles
    ]
    pf6_minus = topology.Molecule.from_smiles("F[P-](F)(F)(F)(F)F")
    pf6_minus.partial_charges = np.asarray([-0.39, 1.34, -0.39, -0.39, -0.39, -0.39, -0.39]) * elementary_charge
    library_charge_type = smirnoff.parameters.LibraryChargeHandler.LibraryChargeType.from_molecule(pf6_minus)
    openff_forcefield = smirnoff.ForceField("openff_unconstrained-2.0.0.offxml")
    openff_forcefield["LibraryCharges"].add_parameter(parameter=library_charge_type)
    openff_topology = openff.toolkit.topology.Topology.from_openmm(
        topology, openff_mols
    )
    system = openff_forcefield.create_openmm_system(openff_topology)
    structure = parmed.openmm.load_topology(topology, system)
    structure.box = make_box_list(box_size, 'openmm')
    system = structure.createSystem(nonbondedMethod=PME, nonbondedCutoff=1 * nanometer)
    return system, topology


def _smiles_to_simulation(
    smiles,
    counts,
    box_size,
    integrator=None,
    properties=None,
    charge_scaling=1,
    temperature=298,
):
    if integrator is None:
        integrator = LangevinMiddleIntegrator(
            300 * kelvin, 1 / picosecond, 0.001 * picoseconds
        )
    # if not platform_properties:
    #     platform_properties = {}
    topology = _smiles_to_topology(smiles, counts)
    coordinates = _smiles_to_coordinates(smiles, counts, box_size)
    openff_mols = [
        openff.toolkit.topology.Molecule.from_smiles(smile) for smile in smiles
    ]
    pf6_minus = openff.toolkit.topology.Molecule.from_smiles("F[P-](F)(F)(F)(F)F")
    pf6_minus.partial_charges = np.asarray([-0.39, 1.34, -0.39, -0.39, -0.39, -0.39, -0.39]) * charge_scaling * elementary_charge
    pf6_charge_type = smirnoff.parameters.LibraryChargeHandler.LibraryChargeType.from_molecule(pf6_minus)
    li_plus = openff.toolkit.topology.Molecule.from_smiles("[Li+]")
    li_plus.partial_charges = np.asarray([1.0]) * charge_scaling * elementary_charge
    li_charge_type = smirnoff.parameters.LibraryChargeHandler.LibraryChargeType.from_molecule(li_plus)
    openff_forcefield = smirnoff.ForceField("openff_unconstrained-2.0.0.offxml")
    openff_forcefield["LibraryCharges"].add_parameter(parameter=pf6_charge_type)
    openff_forcefield["LibraryCharges"].add_parameter(parameter=li_charge_type)
    openff_topology = openff.toolkit.topology.Topology.from_openmm(
        topology, openff_mols
    )
    openff_topology.box_vectors = [box_size, box_size, box_size] * angstrom
    system = openff_forcefield.create_openmm_system(openff_topology, allow_nonintegral_charges=True)
    # structure = parmed.openmm.load_topology(topology, system)
    # structure.box = make_box_list(box_size, 'openmm')
    # system = structure.createSystem(nonbondedMethod=PME, nonbondedCutoff=1 * nanometer)
    # system.addForce(MonteCarloBarostat(1 * atmosphere, 300 * kelvin, 10))
    platform = Platform.getPlatformByName('OpenCL')
    simulation = Simulation(topology, system, integrator, platform=platform, platformProperties=properties)
    simulation.context.setPositions(coordinates)
    # simulation.context.setPeriodicBoxVectors((4, 0, 0), (0, 4, 0), (0, 0, 4))
    return simulation, platform
