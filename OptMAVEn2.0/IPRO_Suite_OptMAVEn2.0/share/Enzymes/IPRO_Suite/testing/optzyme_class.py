import sys
sys.path.append("/gpfs/scratch/mjg5185/IPRO_Suite/modules/")
from MOLECULES import MoleculeFile
from MOLECULES import DesignGroup
from MOLECULES import OptzymeGroup

mol1 = MoleculeFile("3SQG_A.pdb", \
"/gpfs/home/mjg5185/scratch/IPRO_Suite/IPRO/structures/")[0]
mol2 = MoleculeFile("f430.pdb", \
"/gpfs/home/mjg5185/scratch/IPRO_Suite/IPRO/structures/")[0]
mol3 = MoleculeFile("3SQG_B.pdb", \
"/gpfs/home/mjg5185/scratch/IPRO_Suite/IPRO/structures/")[0]

molecules = [mol1, mol2]
group1 = DesignGroup(1, molecules)
group1.objective = "eliminate"
molecules = [mol1, mol3]
group2 = DesignGroup(2, molecules)
group2.objective = "improve"
for mol in group1:
    print mol.name
print mol1.name, "yes"
print mol2.name, "no"
print mol3.name, "yes"
sys.exit(0)
optzyme = OptzymeGroup(1, group1)
optzyme2 = OptzymeGroup(2, [group1, group2])
optzyme2[1].objective = "improve"
optzyme2[0].objective = "eliminate"
optzyme3 = optzyme2.duplicate()
optzyme2[1].objective = "improve"
optzyme2[0].objective = "eliminate"
