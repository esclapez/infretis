from openff.interchange import Interchange
from openff.toolkit import ForceField, Molecule
from openff.units import unit

# load the molecule
mol = Molecule.from_file("../mol.sdf")
topology = mol.to_topology()
topology.box_vectors = unit.Quantity([2.1, 2.1, 2.1], unit.nanometer)
# Load the OpenFF 2.1.0 forcefield called "Sage"
sage = ForceField("openff-2.1.0.offxml")
out = Interchange.from_smirnoff(force_field=sage, topology=topology)
out.to_gro("../mol.gro")
out.to_top("../mol.top")

with open("../gromacs_input/topol.top", "w") as writefile:
    with open("mol.top") as readfile:
        for line in readfile:
            if "[ moleculetype ]" in line:
                writefile.write("; Include tip3p water topology\n")
                writefile.write('#include "amber99.ff/ffnonbonded.itp"\n')
                writefile.write('#include "amber99.ff/tip3p.itp"\n\n')
            writefile.write(line)
