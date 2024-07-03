from Bio.PDB import PDBParser
import sys

pdb_id = sys.argv[1]
input = pdb_id + '.pdb'
output = pdb_id + '.sb2'

def calculate_salt_bridges(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure("pdb", pdb_file)

    charged_residues = {'ARG', 'LYS', 'ASP', 'GLU'}
    charged_atoms = {'ARG': {'NH1': 1, 'NH2': 1, 'NE': 1},
                     'LYS': {'NZ': 1},
                     'ASP': {'OD1': -1, 'OD2': -1},
                     'GLU': {'OE1': -1, 'OE2': -1}}

    intra_salt_bridges = []
    inter_salt_bridges = []

    for model in structure:
        for chain1 in model:
            residues1 = [residue for residue in chain1 if residue.get_resname() in charged_residues]
            for residue2 in residues1:
                for atom1 in residue2:
                    if atom1.get_name() in charged_atoms[residue2.get_resname()]:
                        charge1 = charged_atoms[residue2.get_resname()][atom1.get_name()]
                        for residue3 in residues1:
                            if residue2 != residue3:
                                for atom2 in residue3:
                                    if atom2.get_name() in charged_atoms[residue3.get_resname()]:
                                        charge2 = charged_atoms[residue3.get_resname()][atom2.get_name()]
                                        if charge1 != charge2:
                                            distance = atom1 - atom2
                                            if distance <= 3.2:
                                                salt_bridge = (residue2, residue3, chain1)
                                                intra_salt_bridges.append(salt_bridge)

            for chain2 in model:
                if chain2.id != chain1.id:
                    residues2 = [residue for residue in chain2 if residue.get_resname() in charged_residues]
                    for residue4 in residues1:
                        for atom3 in residue4:
                            if atom3.get_name() in charged_atoms[residue4.get_resname()]:
                                charge3 = charged_atoms[residue4.get_resname()][atom3.get_name()]
                                for residue5 in residues2:
                                    for atom4 in residue5:
                                        if atom4.get_name() in charged_atoms[residue5.get_resname()]:
                                            charge4 = charged_atoms[residue5.get_resname()][atom4.get_name()]
                                            if charge3 != charge4:
                                                distance = atom3 - atom4
                                                if distance <= 3.2:
                                                    salt_bridge = (residue4, residue5, chain1, chain2)
                                                    inter_salt_bridges.append(salt_bridge)

    return intra_salt_bridges, inter_salt_bridges

def print_salt_bridges(salt_bridges, label):
    result=[]
    for bridge in salt_bridges:
        if len(bridge) == 3:
            residue1, residue2, chain = bridge
            pair1 = f"{chain.id}.{residue1.get_resname()}{residue1.get_id()[1]}"
            pair2 = f"{chain.id}.{residue2.get_resname()}{residue2.get_id()[1]}"
            result.append([pair1,pair2])
        elif len(bridge) == 4:
            residue1, residue2, chain1, chain2 = bridge
            pair1 = f"{chain1.id}.{residue1.get_resname()}{residue1.get_id()[1]}"
            pair2 = f"{chain2.id}.{residue2.get_resname()}{residue2.get_id()[1]}"
            result.append([pair1,pair2])
    unique_pair = []
    for pair in result:
        pair.sort()  # Sort the pair to ensure consistent ordering
        if pair not in unique_pair:
            unique_pair.append(pair)

    print(f"{label} salt bridges:",len(unique_pair))
    return unique_pair
# Example usage

intra_salt_bridges, inter_salt_bridges = calculate_salt_bridges(input)
intra = print_salt_bridges(intra_salt_bridges, "Intra")
inter = print_salt_bridges(inter_salt_bridges, "Inter")

file=open(output,'w')
for i in range(len(intra)):
    file.write("%s  %s\n"%(intra[i][0],intra[i][1]))
for i in range(len(inter)):
    file.write("%s  %s\n"%(inter[i][0],inter[i][1]))
file.close()

