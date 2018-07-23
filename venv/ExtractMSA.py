# Extract Multiple Sequence Alignment
# Written by Jim Lennon, based off of earlier work by Arzu Uyar
# Dickson Lab at Michigan State University, May 2018
#
# Given a directory of PDB files, this script will find the residues that fully align between the sequences and then
# output the locations of those atoms to a DCD file for each protein.
#
# To run: python3 <path/to/input/directory> <path/to/output/directory>

from os import listdir, path
import sys
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.PDB import PDBParser, PDBIO
from Bio.AlignIO import ClustalIO
from Bio.Cluster import pca
import mdtraj as md
import tempfile
import numpy as np
import copy
from geomm import theobald_qcp, centroid, superimpose, centering
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

def getInputs(args):
    '''
    Get all the PDB files in the directory
    :param args: arguments from the user
    :return: input directory, output directory and a list of the input PDB files
    '''
    inputs_dir = args.argv[1]

    # Get all PDB input files and put in a list
    input_files = [(inputs_dir + "/" + f) for f in listdir(inputs_dir) \
                   if path.isfile(path.join(inputs_dir, f)) if not f.startswith(".") if f.endswith("pdb")]

    return args.argv[1], args.argv[2], input_files

def getSequences(input_files):
    '''
    Extract the residues from the PDB files
    :param input_files: list of all PDB files
    :return: temporary file pointer to the FASTA file, list of all residues
    '''

    rec_list = []
    tmpFASTAFP = tempfile.NamedTemporaryFile(mode='w+')

    for seq in input_files:                                 # Open each PDB file
        for record in SeqIO.parse(seq, "pdb-seqres"):       # Get each chain
            if record.annotations["chain"] == "A":          # Only write the "A" Chain
                tmpFASTAFP.write(record.format("fasta"))    # Write as a FASTA file
                rec_list.append(record)
                break

    tmpFASTAFP.flush()

    proteins = getStructures(rec_list, input_files)

    return tmpFASTAFP, proteins

def doMSA(tmpFASTAFP):
    '''
    Run ClustalW multiple sequence alignment
    :param tmpFASTAFP: file pointer to the FASTA file
    :return: file pointer to the temporary alignment file
    '''

    clustalw_exe = "/Users/work/DicksonLab/Code/Clustalw/clustalw2"     # Path to ClustalW installation
    tmpPDBFP = tempfile.NamedTemporaryFile(mode='w+')

    clustalw_cline = ClustalwCommandline(clustalw_exe, infile=tmpFASTAFP.name, outfile=tmpPDBFP.name)
    tmpPDBFP.flush()                                                    # Flush the buffer to make sure
                                                                        # it's been written
    clustalw_cline()

    tmpFASTAFP.close()

    msa = AlignIO.read(tmpPDBFP.name, "clustal")  # Read in the ClustalW file
    tmpPDBFP.close()

    return msa

def getConsensusSequences(msa, proteins):
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

def getStructures(rec_list, input_files):
    parser = PDBParser()
    proteins = {}

    for protein, file in zip(rec_list, input_files):        # Loop thru each protein and PDB file
        structure = parser.get_structure(protein.id, file)  # Get the structure of each protein
        proteins[structure.get_id()] = structure

    return proteins

def keepConsensusResidues(proteins, aligned_indices):
    for protein in proteins.values():
        for chains in protein:
            for chain in chains:
                for residue in list(chain):
                    # Remove all residues that aren't fully aligned
                    residue_id = residue.get_id()[1]

                    if residue_id not in aligned_indices[protein.id.replace(":", "_")]:     # Don't want to keep hetero residues or atoms:
                        chain.detach_child(residue.get_id())

        proteins[protein.get_id()] = protein

    return proteins

def removeNonHeteroAtoms(proteins):
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

def correctIndices(proteins, alignment_indices):
    # Remove any alignment indices that we don't have residue data for
    # Note: the ClustalW alignment uses the known residues, but we don't necessarily have the coordinates
    #       of all of the atoms comprising each residue. As a result, we remove all of the indices of the atoms
    #       that may be in alignment, but we don't have their coordinates

    starting_indices = []
    ending_indices = []

    for protein in proteins.values():
        for chains in protein:
            for chain in chains:
                if chain.get_id() == "A":  # Only keeping 'A' chains
                    starting_indices.append(list(chain)[0].get_id()[1])  # Get the first residue ID
                    ending_indices.append(list(chain)[-1].get_id()[1])  # Get the last residue ID

                    break

    for protein, indices in alignment_indices.items():
        if protein is not "original":
            for index in indices:
                if index < max(starting_indices):
                    alignment_indices[protein].remove(index)

    return alignment_indices

def getAlignedSequence(tmpFASTAFP, proteins):
    '''
    Use the Clustal alignment file to get slices of alignment
    :param tmpPDBFP: file pointer to Clustal alignment output
    :return: Dictionary containing the indexes that have full residue alignment for each protein
    '''

    msa_output = doMSA(tmpFASTAFP)

    alignment_indices = getConsensusSequences(msa_output, proteins)
    proteins = removeNonHeteroAtoms(proteins)

    alignment_indices = correctIndices(proteins, alignment_indices)
    proteins = keepConsensusResidues(proteins, alignment_indices)

    return proteins

def stripAtoms(proteins, missing_carbons):
    for protein in proteins.values():
        for chains in protein:
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

    return proteins

def getCarbonArray(proteins):
    carbon_found = {}
    carbon_coords = {}

    for protein in proteins.values():
        for chains in protein:
            for chain in chains:
                if chain.get_id() == "A":
                    carbon_found[protein.get_id()] = []
                    carbon_coords[protein.get_id()] = []

                    for residue in chain:
                        if residue.has_id("CA") and residue.has_id("CB"):
                            carbon_found[protein.get_id()].append(True)
                            carbon_coords[protein.get_id()].append((residue["CA"].get_coord(), residue["CB"].get_coord()))
                        else:
                            carbon_found[protein.get_id()].append(False)
                            carbon_coords[protein.get_id()].append((None))

    # Need to know the number of elements residue elements for indexing below
    max_length = len(carbon_found[list(carbon_found.keys())[0]])

    # Dictionary comprehension source: https://stackoverflow.com/questions/14507591/python-dictionary-comprehension
    coords_list = {key: [] for key in carbon_found}
    key_list = [key for key in carbon_found]

    missing_carbons = []

    for i in range(max_length):
        # Only use carbon atoms if they are both in all proteins
        if all(item[i] is True for item in carbon_found.values()):
            for key in key_list:
                coords_list[key].append(carbon_coords[key][i][0])
                coords_list[key].append(carbon_coords[key][i][1])
        else:
            missing_carbons.append(i)
            print("Warning: Carbon atoms coordinates are inconsistent at index ", i)

    proteins = stripAtoms(proteins, missing_carbons)

    carbon_array = {key : np.asarray(value, dtype=np.float64) for key, value in coords_list.items()}

    return carbon_array, proteins

def superimposeProteins(proteins, carbon_array):
    names = [key for key in carbon_array]

    ref_name = names[0]
    carbon_array[ref_name] = centering.center(carbon_array[ref_name])

    for i in range(1, len(names)):
        moving_name = names[i]
        carbon_array[moving_name] = superimpose.superimpose(carbon_array[ref_name], carbon_array[moving_name])

    for protein_name, protein in proteins.items():
        for chains in protein:
            for chain in chains:
                if chain.get_id() == "A":
                    i = 0

                    for residue in chain:
                        for atom in residue:
                            if atom.get_id() == "CA":
                                atom.set_coord(carbon_array[protein_name][i])
                            elif atom.get_id() == "CB":
                                atom.set_coord(carbon_array[protein_name][i + 1])   # CB's are always after the CA
                            else:
                                residue.detach_child(atom.get_id())

                        i += 1

    return proteins

def getAlignedStructure(proteins):
    '''
    Create a dict containing only the residues that have full alignment
    :param rec_list: list of protein records
    :param input_files: list of all input PDB files
    :return: PDB structure containing only fully aligned residues
    '''

    carbon_array, proteins = getCarbonArray(proteins)

    '''
    translations = getTranslations(carbon_array)
    rotation_matrices = getRotationMatrix(carbon_array)
    proteins = superimposeStructures(proteins, translations, rotation_matrices)
    '''

    superimposeProteins(proteins, carbon_array)

    return carbon_array

def writeDCDFile(rec_list, input_files, output_dir):
    '''
    Write a DCD file containing the locations of the atoms whose residues have full alignment
    :param rec_list: list of protein records
    :param input_files: list of all input PDB files
    :param output_dir: output directory for .dcd files
    '''

    parser = PDBParser()
    io = PDBIO()

    for protein, file in zip(rec_list, input_files):                    # Loop thru each protein and PDB file
        structure = parser.get_structure(protein.id, file)              # Get the structure of each protein

        with tempfile.NamedTemporaryFile(mode='w+') as fp:              # Using temp file to store partial PDB
            io.save(fp.name)
            traj = md.load_pdb(fp.name)                                 # Convert PDB to DCD file
            traj.save_dcd(output_dir + "/" + structure.get_id()[:4] + '.dcd', True)

def draw_vector(v0, v1, ax=None):
    ax = ax or plt.gca()
    arrowprops=dict(arrowstyle='->',
                    linewidth=2,
                    shrinkA=0, shrinkB=0)
    ax.annotate('', v1, v0, arrowprops=arrowprops)



def getPCA(carbon_arrays):
    pcas = {}
    x = []
    y = []

    for key, coordinates in carbon_arrays.items():
        pcas[key] = PCA(n_components=2)
        pcas[key].fit(coordinates)

        x.append(pcas[key].components_[0][0])
        y.append(pcas[key].components_[0][1])

        #plt.scatter(coordinates[:, 0], coordinates[:, 1], alpha=0.2)

        '''
        for length, vector in zip(pcas[key].explained_variance_, pcas[key].components_):
            v = vector * 3 * np.sqrt(length)
            draw_vector(pcas[key].mean_, pcas[key].mean_ + v)
        plt.axis('equal');
        '''

    plt.scatter(x, y, alpha=0.2)
    plt.show()

    return pcas

def writeToPDB(proteins):
    io = PDBIO()

    for protein in proteins.values():
        io.set_structure(protein)
        io.save(output_dir + "/" + protein.get_id() + '_aligned.pdb')

inputs_dir, output_dir, input_files = getInputs(sys)                    # Get the provided inputs
tmpFASTAFP, proteins = getSequences(input_files)                        # Get the sequences

aligned_proteins = getAlignedSequence(tmpFASTAFP, proteins)             # Get the fully aligned residues
carbon_array = getAlignedStructure(aligned_proteins)                    # Get the structurally aligned residues

pcas = getPCA(carbon_array)

writeToPDB(proteins)

#writeDCDFile(rec_list, input_files, output_dir)                         # Write the aligned atoms to a .dcd file



