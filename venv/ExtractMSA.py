# Extract Multiple Sequence Alignment
# Written by Jim Lennon, based off of earlier work by Arzu Uyar
# Dickson Lab at Michigan State University, May 2018
#
# Given a directory of PDB files, this script will find the residues that fully align between the sequences and then
# output the locations of those atoms to a DCD file for each protein.
#
# To run: python3 <path/to/input/directory> <path/to/output/directory>

from os import listdir
from os.path import isfile, join
import sys
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.PDB import PDBParser, PDBIO, Superimposer
from Bio.AlignIO import ClustalIO
import mdtraj as md
import tempfile
import numpy as np

def getInputs(args):
    '''
    Get all the PDB files in the directory
    :param args: arguments from the user
    :return: input directory, output directory and a list of the input PDB files
    '''
    inputs_dir = args.argv[1]

    # Get all PDB input files and put in a list
    input_files = [(inputs_dir + "/" + f) for f in listdir(inputs_dir) \
                   if isfile(join(inputs_dir, f)) if not f.startswith(".") if f.endswith("pdb")]

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

    return tmpFASTAFP, rec_list

def runClustal(tmpFASTAFP):
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

    return tmpPDBFP

def getConsensusSequences(tmpPDBFP):
    align = AlignIO.read(tmpPDBFP.name, "clustal")  # Read in the ClustalW file
    tmpPDBFP.close()

    alignment_offset = {}

    for protein in align:                                               # Loop thru and get each protein
        alignment_offset[protein.id] = [0]

        for index, letter in enumerate(protein.seq):                    # Compare each letter to the gap
            if index == 0:                                              # Check first letter
                if letter == "-":                                       # Initial letter is a gap
                    alignment_offset[protein.id] = [1]
            else:                                                       # Check (not first) letter
                prev = alignment_offset[protein.id][index - 1]

                if letter == "-":                                       # The letter is a gap, add to offset
                    alignment_offset[protein.id].append(prev + 1)
                else:                                                   # The letter is not a gap
                    alignment_offset[protein.id].append(prev)

    consensus = align.column_annotations["clustal_consensus"]           # Get the consensus line (bottom line in file)
    matches = [i for i, letter in enumerate(consensus) if letter == "*"]    # Get indices of all full alignments

    fully_aligned = {}
    fully_aligned["original"] = matches

    for protein in align:                                               # Get indexes of fully aligned amino acids
        # Generate a list of the actual indices of full alignment for each protein. The alignment offset accounts for
        # the gaps in the Clustal output
        fully_aligned[protein.id] = \
            [match - alignment_offset[protein.id.replace(":", "_")][match] + 1 for match in matches]

    return fully_aligned

def stripHetAtoms(structure):
    for model in structure:
        for chain in model:
            if chain.get_id() == "A":  # Only keeping 'A' chains
                for residue in list(chain):
                    if residue.get_id()[0] != " ":     # Don't want to keep hetero residues or atoms
                        chain.detach_child(residue.get_id())

    return structure

def keepConsensusResidues(rec_list, input_files):
    parser = PDBParser()
    io = PDBIO()

    proteins = {}

    for protein, file in zip(rec_list, input_files):  # Loop thru each protein and PDB file
        structure = parser.get_structure(protein.id, file)  # Get the structure of each protein

        for model in structure:
            for chain in model:
                if chain.get_id() == "A":  # Only keeping 'A' chains
                    for residue in list(chain):
                        # Remove all residues that aren't fully aligned
                        num = residue.get_id()[1]

                        if num not in fullyAligned[protein.id.replace(":", "_")]:
                            chain.detach_child(residue.get_id())
                else:
                    model.detach_child(chain.get_id())  # Remove all non-'A' chains

        io.set_structure(structure)  # Make sure structure is OK
        proteins[structure.get_id()] = structure

    return proteins

def getAlignedSequence(tmpPDBFP):
    '''
    Use the Clustal alignment file to get slices of alignment
    :param tmpPDBFP: file pointer to Clustal alignment output
    :return: Dictionary containing the indexes that have full residue alignment for each protein
    '''

    fully_aligned = getConsensusSequences(tmpPDBFP)

    # Remove any alignment indices that we don't have residue data for
    # Note: the ClustalW alignment uses the known residues, but we don't necessarily have the coordinates
    #       of all of the atoms comprising each residue. As a result, we remove all of the indices of the atoms
    #       that may be in alignment, but we don't have their coordinates

    parser = PDBParser()
    starting_indices = []
    ending_indices = []

    # Find the starting index of known atomic coordinates
    for protein, file in zip(rec_list, input_files):    # Loop thru each protein and PDB file
        structure = parser.get_structure(protein.id, file)  # Get the structure of each protein

        structure = stripHetAtoms(structure)            # Strip off hetero atoms

        for model in structure:
            for chain in model:
                if chain.get_id() == "A":  # Only keeping 'A' chains
                    starting_indices.append(list(chain)[0].get_id()[1])     # Get the first residue ID
                    ending_indices.append(list(chain)[-1].get_id()[1])      # Get the last residue ID

    for protein, indices in fully_aligned.items():
        if protein is not "original":
            for index in indices:
                if index < max(starting_indices):
                    fully_aligned[protein].remove(index)

    return fully_aligned

def getCarbonAtoms(proteins):
    carbon_atoms = {}

    for structure in proteins.values():
        for model in structure:
            for chain in model:
                if chain.get_id() == "A":  # Only keeping 'A' chains
                    carbons = []

                    for residue in list(chain):
                        if residue.has_id("CA"):
                            carbons.append(residue["CA"])
                        elif residue.has_id("CB"):
                            carbons.append(residue["CB"])

            carbon_atoms[structure.get_id()] = carbons

    return carbon_atoms

def getAlignedStructure(rec_list, input_files, fullyAligned):
    '''
    Create a dict containing only the residues that have full alignment
    :param rec_list: list of protein records
    :param input_files: list of all input PDB files
    :return: PDB structure containing only fully aligned residues
    '''

    proteins = keepConsensusResidues(rec_list, input_files)
    carbons = cleanStructures(proteins)

    return proteins

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

inputs_dir, output_dir, input_files = getInputs(sys)                    # Get the provided inputs
tmpFASTAFP, rec_list = getSequences(input_files)                        # Get the sequences
tmpPDBFP = runClustal(tmpFASTAFP)                                       # Run the multiple sequence alignment
fullyAligned = getAlignedSequence(tmpPDBFP)                             # Get the alignment offsets (for gaps)
proteins = getAlignedStructure(rec_list, input_files, fullyAligned)
#writeDCDFile(rec_list, input_files, output_dir)                         # Write the aligned atoms to a .dcd file



