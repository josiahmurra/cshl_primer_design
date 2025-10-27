#!/usr/bin/env python

import sys
import argparse
import os
from datetime import datetime
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

    # Get the current datetime
    now = datetime.now()

    default_output = 'optimus_prymer_output_' + now.strftime("%Y%m%d%H%M%S") + '.txt'

    parser.add_argument('-i', '--input', type=validate_file_path, required=True, help='text file with one NCBI accesion ID per line') # Add text file input
    parser.add_argument('-c', '--catch', type=str, required=False, default='All') # Accept catch option
    parser.add_argument('-o', '--output', type=str, required=False, default=default_output)
    parser.add_argument('-minpl', '--minimumPrimerLength', type=int, required=False, default=20)
    parser.add_argument('-maxpl', '--maximumPrimerLength', type=int, required=False, default=22)

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
    print(f'Retrieving GenBank Records from NCBI Entrez Database.')
    seqlist = read_edirect(textfile)
    print('\n', ('-' * 65), '\n')


    # Perform alignment with muscle3
    print(f'Aligning Sequences with:')
    alignmentlist = align_muscle(seqlist)
    print('\n', ('-' * 65), '\n')

    # get order of accession ids from the sequence_files.txt
    sequence_index = {}
    for index, record in enumerate(alignmentlist):
        sequence_index[record.id] = index

    if args.catch != 'All':
        if args.catch not in sequence_index:
            print(f'{args.catch} is not in {textfile}')
        index_to_catch = sequence_index[args.catch]
        print(f"User has selected --catch {args.catch}")
        print(f"Unique primers will be identified for {args.catch}\n")

    # Generate a list of all possible primers
    primer_dict = generate_primer_dict(seqlist, min_len=args.minimumPrimerLength, max_len=args.maximumPrimerLength)

    # Find shared and unique primers
    shared_primers, unique_primers = identify_primer_sharing(primer_dict)
    print(f"Found {len(shared_primers)} primers shared across all sequences")
    print(f"Found {len(unique_primers)} primers unique to single sequences\n")
    print('\n', ('-' * 65), '\n')

    
    # identify unique loci
    unique_loci_dict = unique_mismatch_targetseq(alignmentlist)

    # screen primers for those present at unique loci
    if args.catch != 'All':
        # create an index conversion dictionary
        index_conversion_dict = convert_target_to_MSA(alignmentlist, index_to_catch) 
        # validate unique primers
        print(f'Identifying unique primers for {args.catch}\n')
        print(f'Unique Sequence Ranges:')
        validated_unique_primers = target_unique_loci(unique_loci_dict, index_to_catch, index_conversion_dict, unique_primers)

    # perform primer pair validation
    if args.catch == 'All':
        valid_primer_pairs = check_primer_pairs(shared_primers)
    else:
        valid_primer_pairs = check_primer_pairs(validated_unique_primers)

    # write primer outputs to text file
    print(f"Writing output to {args.output}\n")
    with open(args.output, 'w') as output_file:
        output_file.write(f"Forward_Primer\tReverse_Primer\n")
        for primer_combo in valid_primer_pairs:
            forward_primer = primer_combo[0].seq
            reverse_primer = primer_combo[1].seq
            output_file.write(f"{forward_primer}\t{reverse_primer}\n")


    # print end of program art
    print("""
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⣠⣤⣤⣶⣶⣶⣶⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣶⣶⣶⣶⣤⣤⣤⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠠⣤⣤⣤⣤⣄⣀⣀⣀⣀⠀⠀⣠⣶⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣶⣤⠀⠀⢀⣀⣀⣀⣀⣤⣤⣤⣤⠄
⠀⣿⣿⣿⣿⣿⣿⣿⣿⣿⠀⠀⢻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡟⠀⠀⣿⣿⣿⣿⣿⣿⣿⣿⣿⠀
⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⡇⠀⠸⣿⣿⣿⣿⣿⣿⣿⣯⣍⡉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⢉⣩⣿⣿⣿⣿⣿⣿⣿⣿⠃⠀⢰⣿⣿⣿⣿⣿⣿⣿⣿⡇⠀
⠀⠈⣿⣿⣿⣿⠻⢿⣿⣿⣿⡀⠀⠙⠻⣿⣿⣿⣿⣿⣿⣿⣿⣶⣄⠀⠀⠀⠀⠀⠀⠀⢀⣠⣶⣿⣿⣿⣿⣿⣿⣿⡿⠟⠉⠀⢀⣼⣿⣿⣿⠟⢻⣿⣿⣿⠁⠀
⠀⠀⢹⣿⣿⣿⣄⡀⠉⠻⣿⣿⣷⣤⣀⠀⠙⠻⣿⣿⣿⣿⣿⣿⣿⣿⣦⣄⠀⢀⣠⣶⣿⣿⣿⣿⣿⣿⣿⣿⠟⠋⠀⣀⣤⣾⣿⣿⠟⠋⢀⣠⣿⣿⣿⡏⠀⠀
⠀⠀⠘⣿⣿⣿⣿⣿⣶⣄⡀⠉⠻⣿⣿⣷⣦⣄⠀⠙⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠟⠋⠀⣠⣴⣾⣿⣿⠟⠋⢀⣠⣶⣿⣿⣿⣿⣿⠃⠀⠀
⠀⠀⠀⢻⣿⣿⣧⠈⠛⢿⣿⣶⣄⡀⠙⠻⣿⣿⣿⣶⣄⡈⠙⠿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠟⠋⢀⣠⣴⣿⣿⣿⠟⠋⠀⣠⣶⣿⡿⠛⠁⣸⣿⣿⡟⠀⠀⠀
⠀⠀⠀⠸⣿⣿⣿⣷⣦⣀⠈⠛⢿⣿⣶⣄⡀⠙⣿⣿⣿⣷⠀⠀⡈⠙⠿⣿⣿⣿⣿⠿⠋⢀⠀⠀⣾⣿⣿⣿⠋⠀⣠⣶⣿⡿⠛⠁⢀⣤⣾⣿⣿⣿⠃⠀⠀⠀
⠀⠀⠀⠀⢿⣿⣿⣿⣿⣿⣷⣦⣀⠈⠛⢿⣿⣶⣼⣿⣿⣿⡀⠀⢹⣶⣄⠈⠙⠋⠁⣠⣴⡿⠀⢀⣿⣿⣿⣿⣴⣿⡿⠟⠉⢀⣤⣾⣿⣿⣿⣿⣿⡿⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠙⠻⣿⣿⣿⣿⣿⣿⣷⣦⣀⠈⠛⢿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣦⣴⣿⣿⣿⡇⠀⢸⣿⣿⣿⡿⠛⠉⢀⣤⣾⣿⣿⣿⣿⣿⣿⠿⠋⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⡀⠈⠙⢿⣿⣿⣿⣿⣿⣿⣷⣦⡀⣿⣿⣿⣷⠀⠈⣿⣿⣿⣿⣿⣿⣿⣿⠇⠀⣾⣿⣿⣿⢀⣤⣾⣿⣿⣿⣿⣿⣿⡿⠛⠁⢀⡀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⣿⣷⣤⡀⠉⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡀⠀⣿⣿⣿⣿⣿⣿⣿⣿⠀⢀⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠟⠉⢀⣠⣶⣿⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⣿⣿⣿⣿⣶⠀⠀⠀⠈⠉⠉⠛⠛⠛⠿⠿⢿⡇⠀⢹⣿⣿⣿⣿⣿⣿⣿⠀⢸⡿⠿⠿⠛⠛⠛⠉⠉⠁⠀⠀⠀⣶⣿⣿⣿⣿⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⣿⣿⣿⣿⣿⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⣿⣿⣿⣿⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣿⣿⣿⣿⣿⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⢸⣿⣿⣿⣿⣿⣷⣄⡀⠀⠀⠀⠀⢀⣠⣴⣾⡇⠀⢸⣿⣿⣿⣿⣿⣿⣿⠀⢸⣷⣦⣄⡀⠀⠀⠀⠀⠀⣠⣶⣿⣿⣿⣿⣿⡟⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⡆⠀⢰⣾⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⣿⣿⠀⢸⣿⣿⣿⣿⣷⡆⠀⢰⣿⣿⣿⣿⣿⣿⣿⣿⡇⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⣿⣿⠀⢸⣿⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⡇⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠈⣿⣿⣿⣿⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⣿⣿⠀⢸⣿⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⠇⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⣿⣿⣿⣿⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⣿⣿⠀⢸⣿⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⣿⣿⣿⣿⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⡇⠀⠸⠿⠿⠿⠿⠿⠿⠿⠀⢸⣿⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⢿⣿⣿⣿⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⣇⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣸⣿⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⣿⣿⡿⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠙⠿⣿⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⡿⠋⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣿⣿⣿⡇⠀⢸⣿⣿⣿⣿⣿⡟⠛⠛⠛⠛⠛⠛⠛⠛⠛⠛⢻⣿⣿⣿⣿⣿⡇⠀⢸⣿⣿⣿⠟⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⡇⠀⢸⣿⣿⣿⣿⡿⠀⠀⣶⣶⣶⣶⣶⣶⣶⣶⡀⠀⢻⣿⣿⣿⣿⡇⠀⢸⡿⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⢸⣿⣿⣿⡿⠁⠀⣼⣿⣿⣿⣿⣿⣿⣿⣿⣧⠀⠈⢿⣿⣿⣿⡇⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣿⣿⠃⠀⣼⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣇⠀⠈⣿⣿⠟⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠃⠀⣼⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣆⠀⠈⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀

   
          "Freedom is the right of all sentient beings" 
""")




if __name__ == '__main__':
    main()