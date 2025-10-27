#!/usr/bin/env python

import sys
import argparse
import os
from Bio.Seq import reverse_complement
from step1_read_edirect import read_edirect
from step1_read_edirect import align_muscle
from primer_generation import generate_primer_dict, identify_primer_sharing
from step3_locate_primer_indices import convert_target_to_MSA, unique_mismatch_targetseq, target_unique_loci
from evaluate_primer_pairs import check_primer_pairs


def validate_file_path(file_path):
    """
    Custom type function for argparse to validate if a file exists.
    """
    if not os.path.exists(file_path):
        raise argparse.ArgumentTypeError(f"File '{file_path}' does not exist.")
    return file_path


def main():
    # Parse input with argparse
    parser = argparse.ArgumentParser(
        description = f"OptimusPrymer identifies common and unique primers for gene transcripts. Usage is the following:\n"
                      f"./optimus_prymer_main --input <text_file> --catch <NCBI_ID>",
        formatter_class=argparse.RawTextHelpFormatter
    ) # Create parser

    parser.add_argument('-i', '--input', type=validate_file_path, required=True, help='text file with one NCBI accesion ID per line') # Add text file input
    parser.add_argument('-c', '--catch', type=str, required=False, default='All') # Accept catch option
    args = parser.parse_args() # Parse the arguments

    # Print welcome message
    print(f"""

    ███████╗ ██████╗ ████████╗██╗███╗   ███╗██╗   ██╗███████╗
    ██╔═══██╗██╔══██╗╚══██╔══╝██║████╗ ████║██║   ██║██╔════╝
    ██║   ██║██████╔╝   ██║   ██║██╔████╔██║██║   ██║███████╗
    ██║   ██║██╔═══╝    ██║   ██║██║╚██╔╝██║██║   ██║╚════██║
    ╚██████╔╝██║        ██║   ██║██║ ╚═╝ ██║╚██████╔╝███████║
     ╚═════╝ ╚═╝        ╚═╝   ╚═╝╚═╝     ╚═╝ ╚═════╝ ╚══════╝
    
    ██████╗ ██████╗ ██╗   ██╗███╗   ███╗███████╗██████╗ 
    ██╔══██╗██╔══██╗╚██╗ ██╔╝████╗ ████║██╔════╝██╔══██╗
    ██████╔╝██████╔╝ ╚████╔╝ ██╔████╔██║█████╗  ██████╔╝
    ██╔═══╝ ██╔══██╗  ╚██╔╝  ██║╚██╔╝██║██╔══╝  ██╔══██╗
    ██║     ██║  ██║   ██║   ██║ ╚═╝ ██║███████╗██║  ██║
    ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚═╝     ╚═╝╚══════╝╚═╝  ╚═╝ 
    
    """)

    # Check user input
    textfile = args.input 

    # Retrieve genbank records with edirect
    print(f'Retrieving GenBank Records\n')
    seqlist = read_edirect(textfile)

    # Perform alignment with muscle3
    print(f'Aligning Sequences with MUSCLE')
    alignmentlist = align_muscle(seqlist)
    print(f'\n')

    # get order of accession ids from the sequence_files.txt
    sequence_index = {}
    for index, record in enumerate(alignmentlist):
        sequence_index[record.id] = index

    if args.catch != 'All':
        index_to_catch = sequence_index[args.catch]

    # Generate a list of all possible primers
    primer_dict = generate_primer_dict(seqlist, min_len=20, max_len=21)

    # Find shared and unique primers
    shared_primers, unique_primers = identify_primer_sharing(primer_dict)
    print(f"Found {len(shared_primers)} primers shared across all sequences")
    print(f"Found {len(unique_primers)} primers unique to single sequences\n")

     # create an index conversion dictionary
    index_conversion_dict = convert_target_to_MSA(alignmentlist, 2) # just test 2 index for now

    # identify unique loci
    unique_loci_dict = unique_mismatch_targetseq(alignmentlist)

    # screen primers for those present at unique loci
    if args.catch != 'All':
        print(f'Identifying unique primers for {args.catch}\n')
        print(f'Unique Sequence Ranges:')
        validated_unique_primers = target_unique_loci(unique_loci_dict, 2, index_conversion_dict, unique_primers)

    # perform primer pair validation
    if args.catch == 'All':
        valid_primer_pairs = check_primer_pairs(shared_primers)
    else:
        valid_primer_pairs = check_primer_pairs(validated_unique_primers)

    




if __name__ == '__main__':
    main()