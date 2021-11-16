from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.typing.engines.smirnoff.parameters import LibraryChargeHandler

import numpy as np
from openmm import NonbondedForce, unit

from openff.toolkit.topology import Molecule, Topology

from openmm import XmlSerializer
from openmm.app import PDBFile
from openmm.unit import angstrom
from setup_functions import _smiles_to_topology

pf6_minus = Molecule.from_smiles("F[P-](F)(F)(F)(F)F")
li_plus = Molecule.from_smiles("[Li+]")


pf6_minus.partial_charges = np.asarray([-0.39, 1.34, -0.39, -0.39, -0.39, -0.39, -0.39]) * unit.elementary_charge
li_plus.partial_charges = np.asarray([0.99]) * unit.elementary_charge

library_charge_type = LibraryChargeHandler.LibraryChargeType.from_molecule(pf6_minus)
li_charge_type = LibraryChargeHandler.LibraryChargeType.from_molecule(li_plus)

forcefield = ForceField("openff_unconstrained-2.0.0.offxml")
forcefield["LibraryCharges"].add_parameter(parameter=library_charge_type)
forcefield["LibraryCharges"].add_parameter(parameter=li_charge_type)
topology = pf6_minus.to_topology()
li_topology = li_plus.to_topology()

system = forcefield.create_openmm_system(topology)
li_system = forcefield.create_openmm_system(li_topology, allow_nonintegral_charges=True)
# Inspect the openmm.System object to see which charges
# actually ended up on each particle
for force in li_system.getForces():
    if type(force) == NonbondedForce:
        for i in range(force.getNumParticles()):
            print(force.getParticleParameters(i)[0])






print('hi')