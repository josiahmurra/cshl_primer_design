#!/usr/bin/env python

import argparse
import os
from Bio.Seq import reverse_complement
from step1_read_edirect import read_edirect, align_muscle
from primer_check import check_single_primer
from primer_generation import PrimerObject, generate_primer_dict, identify_primer_sharing


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
    seqlist = read_edirect(textfile)

    # Generate a list of all possible primers
    primer_dict = generate_primer_dict(seqlist)

    # Find shared and unique primers
    shared_primers, unique_primers = identify_primer_sharing(primer_dict)
    print(f"Found {len(shared_primers)} primers shared across all sequences")
    print(f"Found {len(unique_primers)} primers unique to single sequences")



if __name__ == '__main__':
    main()