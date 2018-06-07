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
from Bio.PDB import PDBParser, PDBIO
from Bio.AlignIO import ClustalIO
import mdtraj as md
import tempfile

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
                tmpFASTAFP.write(record.format("fasta"))   # Write as a FASTA file

                # Uncomment the line below if you would like to see the protein's code and sequence
                #print(record.format("fasta"))

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

def getAlignment(tmpPDBFP):
    '''
    Use the Clustal alignment file to get slices of alignment
    :param tmpPDBFP: file pointer to Clustal alignment output
    :return: Dictionary containing the indexes that have full residue alignment for each protein
    '''

    alignmentOffset = {}

    align = AlignIO.read(tmpPDBFP.name, "clustal")                      # Read in the ClustalW file
    tmpPDBFP.close()

    for protein in align:                                               # Loop thru and get each protein
        alignmentOffset[protein.id] = [0]

        for index, letter in enumerate(protein.seq):                    # Compare each letter to the gap
            if index == 0:                                              # Check first letter
                if letter == "-":                                       # Initial letter is a gap
                    alignmentOffset[protein.id] = [1]
            else:                                                       # Check (not first) letter
                prev = alignmentOffset[protein.id][index - 1]

                if letter == "-":                                       # The letter is a gap, add to offset
                    alignmentOffset[protein.id].append(prev + 1)
                else:                                                   # The letter is not a gap
                    alignmentOffset[protein.id].append(prev)

    consensus = align.column_annotations["clustal_consensus"]           # Get the consensus line (bottom line in file)
    matches = [i for i, letter in enumerate(consensus) if letter == "*"]    # Get indices of all full alignments

    fullyAligned = {}

    for protein in align:                                               # Get indexes of fully aligned amino acids
        # Generate a list of the actual indices of full alignment for each protein. The alignment offset accounts for
        # the gaps in the Clustal output
        fullyAligned[protein.id] = [match - alignmentOffset[protein.id.replace(":", "_")][match] + 1 for match in matches]

    return fullyAligned

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

        for model in structure:
            for chain in model:
                if chain.get_id() == "A":                               # Only keeping 'A' chains
                    for residue in list(chain):
                        # Remove all residues that aren't fully aligned
                        if residue.get_id()[1] not in fullyAligned[protein.id.replace(":", "_")]:
                            chain.detach_child(residue.get_id())
                else:
                    model.detach_child(chain.get_id())                  # Remove all non-'A' chains

        io.set_structure(structure)                                     # Make sure structure is OK

        with tempfile.NamedTemporaryFile(mode='w+') as fp:              # Using temp file to store partial PDB
            io.save(fp.name)
            traj = md.load_pdb(fp.name)                                 # Convert PDB to DCD file
            traj.save_dcd(output_dir + "/" + structure.get_id()[:4] + '.dcd', True)

inputs_dir, output_dir, input_files = getInputs(sys)                    # Get the provided inputs
tmpFASTAFP, rec_list = getSequences(input_files)                        # Get the sequences
tmpPDBFP = runClustal(tmpFASTAFP)                                       # Run the multiple sequence alignment
fullyAligned = getAlignment(tmpPDBFP)                                   # Get the alignment offsets (for gaps)
writeDCDFile(rec_list, input_files, output_dir)                         # Write the aligned atoms to a .dcd file