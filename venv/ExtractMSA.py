# Extract Multiple Sequence Alignment
# Written by Jim Lennon, based off of earlier work by Arzu Uyar
# Dickson Lab at Michigan State University, May 2018
#
# Given a directory of PDB files, this script will find the residues that fully align between the sequences and then
# output the locations of those atoms to a DCD file for each protein.
#
# To run: python3 <path/to/input/directory> <path/to/output/directory>

import tempfile
import sys
import copy
from os import listdir, path, devnull

import numpy as np
#import mdtraj as md
import matplotlib.pyplot as plt
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.PDB import PDBParser, PDBIO
from Bio.AlignIO import ClustalIO
from Bio.Cluster import pca
from geomm import theobald_qcp, centroid, superimpose, centering
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize


class HiddenPrints:
    '''
    Found at https://stackoverflow.com/questions/8391411/suppress-calls-to-print-python
    '''

    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


class Protein:
    def __init__(self, name, chain, bio_struct=None, carbons=None, pdb_file=None, pca1=None, pca2=None):
        self._name = name
        self._chain = chain
        self._bio_struct = bio_struct
        self._pdb_file = pdb_file
        self._carbons = carbons
        self._pca1 = pca1
        self._pca2 = pca2

    def get_name(self):
        return self._name

    def get_chain(self):
        return self._chain

    def get_bio_struct(self):
        return self._bio_struct

    def get_pdb_file(self):
        return self._pdb_file

    def get_carbons(self):
        return self._carbons

    def get_pca1(self):
        return self._pca1

    def get_pca2(self):
        return self._pca2

    def set_name(self, name):
        self._name = name

    def set_chain(self, chain):
        self._chain = chain

    def set_bio_struct(self, bio_struct):
        self._bio_struct = bio_struct

    def set_pdb_file(self, pdb_file):
        self._pdb_file = pdb_file

    def set_carbons(self, carbons):
        self._carbons = carbons

    def set_pca1(self, pca1):
        self._pca1 = pca1

    def set_pca2(self, pca2):
        self._pca2 = pca2

def get_inputs(args):
    '''
    Get all the PDB files in the directory
    :param args: arguments from the user
    :return: input directory, output directory location of Clustal W executable
                and a list of the input PDB files
    '''
    inputs_dir = args.argv[1]

    # Get all PDB input files and put in a list
    input_files = [(inputs_dir + "/" + f) for f in listdir(inputs_dir) \
                   if path.isfile(path.join(inputs_dir, f)) if not f.startswith(".") if f.endswith("pdb") or f.endswith(".ent")]

    assert len(input_files) > 1, "There must be more than one PDB file in the input directory."

    with open(inputs_dir + '/structures_list.txt', 'r') as f:
        read_data = f.read().splitlines()

    structure_dict = {}

    NAME = 0 # Constant name of the protein
    CHAIN = 1 # Constant name of the chain to use e.g. "A", "B", etc...

    for i in range(2, len(read_data)): # First two lines are column headers
        split_row = read_data[i].split()

        protein_name = split_row[NAME].upper()
        #tmp = split_row[CHAIN]
        structure_dict[protein_name] = Protein(protein_name, split_row[CHAIN])

    return args.argv[1], args.argv[2], args.argv[3], input_files, structure_dict

def get_sequences(input_files, structure_dict):
    '''
    Extract the residues from the PDB files
    :param input_files: list of all PDB files
    :return: temporary file pointer to the FASTA file, list of all residues
    '''

    rec_list = []
    tmpFASTAFP = tempfile.NamedTemporaryFile(mode='w+')

    for seq in input_files:                                 # Open each PDB file
        #with HiddenPrints():
        for record in SeqIO.parse(seq, "pdb-seqres"):       # Get each chain
            # First four letters are the structure (names are in the format [NAME:C], with "C" being the chain ID)
            structure_name = record.name[:4]

            if structure_name == "<unk":
                continue

            if record.annotations["chain"] == structure_dict[structure_name].get_chain():
                tmpFASTAFP.write(record.format("fasta"))    # Write as a FASTA file
                rec_list.append(record)
                break

    tmpFASTAFP.flush()

    get_structures(rec_list, input_files, structure_dict)

    return tmpFASTAFP

def do_msa(tmpFASTAFP, clustal_w_loc):
    '''
    Run ClustalW multiple sequence alignment
    :param tmpFASTAFP: file pointer to the FASTA file
    :return: file pointer to the temporary alignment file
    '''

    #clustal_w_loc = "/Users/jim/DicksonLab/Code/Clustalw/clustalw2"     # Path to ClustalW installation
    tmpPDBFP = tempfile.NamedTemporaryFile(mode='w+')

    clustalw_cline = ClustalwCommandline(clustal_w_loc, infile=tmpFASTAFP.name, outfile=tmpPDBFP.name)
    tmpPDBFP.flush()                                                    # Flush the buffer to make sure
                                                                        # it's been written
    clustalw_cline()

    tmpFASTAFP.close()

    msa = AlignIO.read(tmpPDBFP.name, "clustal")  # Read in the ClustalW file
    tmpPDBFP.close()

    return msa

def get_consensus_sequences(msa, proteins):
    alignment_offset = {}

    for alignment in msa:                                               # Loop thru and get each protein
        for index, letter in enumerate(alignment.seq):                  # Compare each letter to the gap
            if index == 0:                                              # Check first letter
                if letter == "-":                                       # Initial letter is a gap
                    alignment_offset[alignment.id] = [1]
                else:
                    alignment_offset[alignment.id] = [0]
            else:                                                       # Check (not first) letter
                prev = alignment_offset[alignment.id][index - 1]

                if letter == "-":                                       # The letter is a gap, add to offset
                    alignment_offset[alignment.id].append(prev + 1)
                else:                                                   # The letter is not a gap
                    alignment_offset[alignment.id].append(prev)

    consensus = msa.column_annotations["clustal_consensus"]           # Get the consensus line (bottom line in file)
    matches = [i for i, letter in enumerate(consensus) if letter == "*"]    # Get indices of all full alignments

    alignment_indices = {}
    alignment_indices["original"] = matches

    for protein in msa:                                               # Get indexes of fully aligned amino acids
        # Generate a list of the actual indices of full alignment for each protein. The alignment offset accounts for
        # the gaps in the Clustal output
        alignment_indices[protein.id] = \
            [match - alignment_offset[protein.id.replace(":", "_")][match] + 1 for match in matches]

    return alignment_indices

def get_structures(rec_list, input_files, structure_dict):
    parser = PDBParser()

    for protein, file in zip(rec_list, input_files):        # Loop thru each protein and PDB file
        structure = parser.get_structure(protein.id, file)  # Get the structure of each protein

        structure_name = structure.get_id()[:4].upper()
        structure_dict[structure_name].set_bio_struct(structure)

def keep_consensus_residues(structure_dict, aligned_indices):
    for structure in structure_dict.values():
        for chains in structure.get_bio_struct():
            for chain in chains:
                for residue in list(chain):
                    # Remove all residues that aren't fully aligned
                    residue_id = residue.get_id()[1]

                    if residue_id not in aligned_indices[structure.get_bio_struct().id.replace(":", "_")]:     # Don't want to keep hetero residues or atoms:
                        chain.detach_child(residue.get_id())

        #proteins[protein.get_id()] = protein

    #return proteins

def remove_non_hetero_atoms(proteins):
    for protein in proteins.values():
        for chains in protein:
            for chain in chains:
                if chain.get_id() == "A":  # Only keeping 'A' chains
                    for residue in list(chain):
                        residue_flag = residue.get_id()[0]

                        if residue_flag != " ":     # Don't want to keep hetero residues or atoms:
                            chain.detach_child(residue.get_id())
                else:
                    chains.detach_child(chain.get_id())  # Remove all non-'A' chains

        proteins[protein.id] = protein

    return proteins

def correct_indices(structure_dict, alignment_indices):
    # Remove any alignment indices that we don't have residue data for
    # Note: the ClustalW alignment uses the known residues, but we don't necessarily have the coordinates
    #       of all of the atoms comprising each residue. As a result, we remove all of the indices of the atoms
    #       that may be in alignment, but we don't have their coordinates

    starting_indices = []
    ending_indices = []

    for structure in structure_dict.values():
        for chains in structure.get_bio_struct():
            for chain in chains:
                if chain.get_id() == structure.get_chain():
                    starting_indices.append(list(chain)[0].get_id()[1])  # Get the first residue ID
                    ending_indices.append(list(chain)[-1].get_id()[1])  # Get the last residue ID

                    break

    for protein, indices in alignment_indices.items():
        if protein is not "original":
            for index in indices:
                if index < max(starting_indices):
                    alignment_indices[protein].remove(index)

    return alignment_indices

def get_aligned_sequence(tmpFASTAFP, structure_dict, clustal_w_loc):
    '''
    Use the Clustal alignment file to get slices of alignment
    :param tmpPDBFP: file pointer to Clustal alignment output
    :return: Dictionary containing the indexes that have full residue alignment for each protein
    '''

    msa_output = do_msa(tmpFASTAFP, clustal_w_loc)

    alignment_indices = get_consensus_sequences(msa_output, structure_dict)
    #remove_non_hetero_atoms(proteins)

    alignment_indices = correct_indices(structure_dict, alignment_indices)
    keep_consensus_residues(structure_dict, alignment_indices)

    #return proteins

def strip_atoms(structure_dict, missing_carbons):
    for structure in structure_dict.values():
        for chains in structure.get_bio_struct():
            for chain in chains:
                if chain.get_id() == "A":  # Only keeping 'A' chains
                    i = 0

                    for residue in list(chain):
                        if i in missing_carbons:
                            chain.detach_child(residue.get_id())
                        else:
                            for atom in residue:
                                if atom.get_id() not in ("CA", "CB"):
                                    residue.detach_child(atom.get_id())

                        i += 1

def get_carbon_array(structure_dict):
    carbon_found = {}
    carbon_coords = {}

    for structure in structure_dict.values():
        for chains in structure.get_bio_struct():
            for chain in chains:
                if chain.get_id() == structure.get_chain():
                    carbon_found[structure.get_name()] = []
                    carbon_coords[structure.get_name()] = []

                    for residue in chain:
                        if residue.has_id("CA") and residue.has_id("CB"):
                            carbon_found[structure.get_name()].append(True)
                            carbon_coords[structure.get_name()].append((residue["CA"].get_coord(), residue["CB"].get_coord()))
                        else:
                            carbon_found[structure.get_name()].append(False)
                            carbon_coords[structure.get_name()].append((None))

        #print("Coords: ", structure.get_name(), carbon_coords[structure.get_name()])

    # Need to know the number of elements residue elements for indexing below
    max_length = len(carbon_found[list(carbon_found.keys())[0]])

    # Dictionary comprehension source: https://stackoverflow.com/questions/14507591/python-dictionary-comprehension
    coords_list = {key: [] for key in carbon_found}
    key_list = [key for key in carbon_found]

    '''
    carbon_lens = [len(carbon_len) for carbon_len in carbon_coords]
    max_carbon_lens = max(carbon_lens)

    to_delete = []

    for name, values in carbon_coords.items():
        if len(values) != max_carbon_lens:
            to_delete.append(name)

    for name in to_delete:
        del carbon_coords[name]
        del carbon_found[name]
        del coords_list[name]
        del structure_dict[name]
    '''

    missing_carbons = []

    for i in range(max_length):
        # Only use carbon atoms if they are both in all proteins

        print("i =", i)

        if all(i < len(item) and item[i] is True for item in carbon_found.values()):
            for key in key_list:
                coords_list[key].append(carbon_coords[key][i][0])
                coords_list[key].append(carbon_coords[key][i][1])
        else:
            missing_carbons.append(i)
            print("Warning: Carbon atoms coordinates are missing for at least one structure at index ", i)

    strip_atoms(structure_dict, missing_carbons)

    for structure in structure_dict.values():
        structure.set_carbons(coords_list[structure.get_name()])


def superimpose_proteins(structure_dict):
    names = tuple(key for key in structure_dict.keys())

    ref_name = names[0]

    # Move reference structure to the origin
    centered_centroid = centroid.centroid(structure_dict[ref_name].get_carbons())
    centered_ref_coords = [atomic_coords - centered_centroid for atomic_coords in structure_dict[ref_name].get_carbons()]
    structure_dict[ref_name].set_carbons(centered_ref_coords)

    for i in range(1, len(names)):
        moving_carbons = structure_dict[names[i]].get_carbons()

        # Superimpose the structure onto the reference structure and keep only the updated coordinates
        # (do not need rotation matrix or RMSD)
        centered_structure = superimpose.superimpose(np.asarray(centered_ref_coords), np.asarray(moving_carbons))[0]
        structure_dict[names[i]].set_carbons(centered_structure.tolist())

    # TODO update coordinates of all atoms after superimposing (not just carbon alpha and carbon beta atoms)

    '''
    for structure_name, structure in structure_dict.items():
        for chains in structure.get_bio_struct():
            for chain in chains:
                if chain.get_id() == structure.get_chain():
                    i = 0

                    for residue in chain:
                        for atom in residue:
                            if atom.get_id() == "CA":
                                atom.set_coord(structure[i])
                            elif atom.get_id() == "CB":
                                atom.set_coord(carbon_array[protein_name][i + 1])   # CB's are always after the CA
                            else:
                                residue.detach_child(atom.get_id())

                        i += 1

    return proteins
    '''

def get_aligned_structure(structure_dict):
    '''
    Create a dict containing only the residues that have full alignment
    :param rec_list: list of protein records
    :param input_files: list of all input PDB files
    :return: PDB structure containing only fully aligned residues
    '''

    get_carbon_array(structure_dict)
    superimpose_proteins(structure_dict)

# def write_dcd_file(rec_list, input_files, output_dir):
#     '''
#     Write a DCD file containing the locations of the atoms whose residues have full alignment
#     :param rec_list: list of protein records
#     :param input_files: list of all input PDB files
#     :param output_dir: output directory for .dcd files
#     '''
#
#     parser = PDBParser()
#     io = PDBIO()
#
#     for protein, file in zip(rec_list, input_files):                    # Loop thru each protein and PDB file
#         structure = parser.get_structure(protein.id, file)              # Get the structure of each protein
#
#         with tempfile.NamedTemporaryFile(mode='w+') as fp:              # Using temp file to store partial PDB
#             io.save(fp.name)
#             traj = md.load_pdb(fp.name)                                 # Convert PDB to DCD file
#             traj.save_dcd(output_dir + "/" + structure.get_id()[:4] + '.dcd', True)

def get_atomic_coords(structure_dict):
    # Reorganize carbon coords. into lists for each respective atom

    # Get lengths necessary for looping and initializing the 2D array
    structure_len = len(structure_dict) # count of all proteins
    structure_names = [item for item in structure_dict.keys()]

    # count of all atoms in the first protein (all should have the same number of atoms)
    atom_len = len(structure_dict[structure_names[0]].get_carbons())

    coords = np.zeros(shape=(atom_len, structure_len), dtype=(float, 3))

    #print(coords)

    for name, i in zip(structure_names, range(structure_len)):
        for j in range(atom_len):
            coords[j][i] = structure_dict[name].get_carbons()[j]

    #print(coords)

    return coords


def scaler_dot_product(a, b):
    sum = 0

    for a_value, b_value in zip(a, b):
        sum += (a_value[0] * b_value[0] + a_value[1] * b_value[1] + a_value[2] * b_value[2])

    return sum


def get_projection(structure_dict, pcas1, pcas2):
    projected_pcas1 = []
    projected_pcas2 = []

    for structure in structure_dict.values():
        structure.set_pca1(scaler_dot_product(pcas1, structure.get_carbons()))
        structure.set_pca2(scaler_dot_product(pcas2, structure.get_carbons()))

        projected_pcas1.append(structure.get_pca1())
        projected_pcas2.append(structure.get_pca2())

    return projected_pcas1, projected_pcas2

def get_pca(structure_dict):
    coords = get_atomic_coords(structure_dict)

    pca = PCA(n_components=2)
    pcas1 = []
    pcas2 = []

    for coord in coords:
        pca.fit(coord)
        pcas1.append(pca.components_[0])
        pcas2.append(pca.components_[1])

    #print(pcas1, "\n", pcas2)

    projected_pcas1, projected_pcas2 = get_projection(structure_dict, pcas1, pcas2)

    print_pcas(projected_pcas1, projected_pcas2)

    return projected_pcas1, projected_pcas2

#def add_test_chain(structure_dict):
    #for structure in structure_dict:


def print_pcas(pca1, pca2):
    plt.scatter(pca1, pca2, alpha=0.2)
    plt.xlabel("Projected PC 1")
    plt.ylabel("Projected PC 2")

    plt.show()

def write_to_pdb(structure_dict):
    io = PDBIO()

    for structure in structure_dict.values():
        io.set_structure(structure.get_bio_struct())
        io.save(output_dir + "/" + structure.get_name() + '_shifted.pdb')

inputs_dir, output_dir, clustal_w_loc, input_files, structure_dict = get_inputs(sys)       # Get the provided inputs
tmpFASTAFP = get_sequences(input_files, structure_dict)                     # Get the sequences

get_aligned_sequence(tmpFASTAFP, structure_dict, clustal_w_loc)                            # Get the fully aligned residues
get_aligned_structure(structure_dict)                                       # Get the structurally aligned residues

pcas = get_pca(structure_dict)

write_to_pdb(structure_dict)

#write_dcd_file(rec_list, input_files, output_dir)                          # Write the aligned atoms to a .dcd file



