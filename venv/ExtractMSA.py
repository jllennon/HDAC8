from os import listdir
from os.path import isfile, join
import sys
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.PDB import PDBParser
from Bio.AlignIO import ClustalIO

def getInputs(args):                                        # Get all the PDB files in the directory
    inputs_dir = args.argv[1]

    # Get all PDB input files and put in a list
    input_files = [(inputs_dir + "/" + f) for f in listdir(inputs_dir) \
                   if isfile(join(inputs_dir, f)) if not f.startswith(".") if f.endswith("pdb")]

    return args.argv[1], args.argv[2], input_files

def getSequences(inputs_dir):
    tmp_file = inputs_dir + "/tmp.fa"                       # TODO: Change into python's temporary file
    rec_list = []

    # Read input files and convert to fasta format
    open(tmp_file, "w+").close()                            # Create or empty the file as necessary

    with open(tmp_file, "a") as handle:
        for seq in input_files:                             # Open each PDB file
            for record in SeqIO.parse(seq, "pdb-seqres"):   # Get each chain
                if record.annotations["chain"] == "A":      # Only write the "A" Chain
                    handle.write(record.format("fasta"))    # Write as a FASTA file
                    rec_list.append(record)
                    #print(dir(record))
                    break

    tmp_out_file = inputs_dir + "/output.aln"               # TODO: Change into python's temporary file

    return tmp_file, tmp_out_file, rec_list

def runClustal(tmp_file, tmp_out_file):                     # Run ClustalW multiple sequence alignment
    clustalw_exe = "/Users/work/DicksonLab/Code/Clustalw/clustalw2"             # Path to ClustalW installation
    clustalw_cline = ClustalwCommandline(clustalw_exe, infile=tmp_file, outfile=tmp_out_file)
    clustalw_cline()

def getAlignmentOffsets(tmp_out_file):
    alignmentOffset = {}

    align = AlignIO.read(tmp_out_file, "clustal")                       # Read in the ClustalW file

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
    matches = [i for i, letter in enumerate(consensus) if letter == "*"]    # Get indexes of all full alignments

    return alignmentOffset, matches

inputs_dir, output_file, input_files = getInputs(sys)                   # Get the provided inputs
tmp_file, tmp_out_file, rec_list = getSequences(inputs_dir)             # Get the sequences
runClustal(tmp_file, tmp_out_file)                                      # Run the multiple sequence alignment
alignmentOffset, matches = getAlignmentOffsets(tmp_out_file)            # Get the alignment offsets (for gaps)




'''
p = PDBParser()
structure = p.get_structure(rec_list[0].id, input_files[0])


for model in structure:
  for chain in model:
    for residue in chain:
      for atom in residue:
        print(residue, atom, atom.get_coord())
'''





# Save coordinates of alignment to file


#use complete_align_list data to extract xyz coordinates from each pdb
#read all pdb files in the input folder
#you need a loop here***************************************
#mdj.load_pdb(pdbid) for pdbid in pdb_paths

#get initial residue number for each pdb


#using initial residue number, find the correct residues corresponding to the complete_align_list


#save the xyz of correct residues into a single file