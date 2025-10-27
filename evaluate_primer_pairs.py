#!/usr/bin/env python

import sys
import re
from Bio.Seq import Seq
from Bio import Align
from step1_read_edirect import read_edirect
from primer_generation import generate_primer_dict
from primer_generation import identify_primer_sharing


# return suitable reverse primers for shared forward primers
# start with a list of objects in the primer class

def main():
    textfile = sys.argv[1]
    seqlist = read_edirect(textfile)
    primer_dict = generate_primer_dict(seqlist)

    # Find shared and unique primers
    shared_primers, unique_primers = identify_primer_sharing(primer_dict)
    print(f"Found {len(shared_primers)} primers shared across all sequences")
    print(f"Found {len(unique_primers)} primers unique to single sequences")


    primer_pairs = check_primer_pairs(shared_primers)

dict_primer_pair = {}
def check_primer_pairs(primer_list):
    aligner = Align.PairwiseAligner() # create a pairwise aligner object
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -10 # this prevents gap extension
    

    forward_primers = list()
    reverse_primers = list()
    compatible_primer_list = list()
    incompatable_tm = 0
    incompatable_distance = 0

    if isinstance(primer_list, dict):
        primer_list = list(primer_list.values())

    for primer in primer_list:

        if primer.direction == "forward":
            forward_primers.append(primer)
        elif primer.direction == "reverse":
            reverse_primers.append(primer)

    for for_primer in forward_primers:
        for rev_primer in reverse_primers:
            max_score = 0
            rev_rev_com = Seq(rev_primer.seq).reverse_complement() # reverse complement the reverse primer

            if abs(for_primer.tm - rev_primer.tm) >5:
                incompatable_tm += 1
                continue

            if 50 > for_primer.start - rev_primer.end < 200:
                incompatable_distance += 1
                continue

            alignments = aligner.align(for_primer.seq,  rev_rev_com) # perform alignment

            for align in alignments:
                if align.score > max_score:
                    max_score = align.score
    
            if max_score > 8:
                continue

            compatible_primer_list.append((for_primer, rev_primer))

    print(f'Found {len(compatible_primer_list)} compatible primers') 
    print(f'bad distance: {incompatable_distance}, bad tm: {incompatable_tm}')     

        

if __name__ == '__main__':
    main()