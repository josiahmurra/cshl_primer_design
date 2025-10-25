#!/usr/bin/env python

import sys
import re
from Bio.Seq import Seq
from Bio import Align

debug = False


def calculate_gc_content(sequence):
    """
    Calculate the GC content of a DNA sequence
    """
    if isinstance(sequence, Seq):
        sequence = str(sequence) # create string object if needed
    sequence = sequence.upper()
    calculated_gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)

    if not 0.4 < calculated_gc_content < 0.6:
        if debug:
            print(f"Warning: GC content of primer is {calculated_gc_content}", file = sys.stderr)
    
    return calculated_gc_content



def calculate_tm(sequence):
    """
    Calculate the melting temperature of a DNA sequence
    """
    if isinstance(sequence, Seq):
        sequence = str(sequence) # create string object if needed
    sequence = sequence.upper()

    # Using tm formula from roslaind.bio
    # Assume sequence length > 14 nucelotides, 50mM Na+, and pH 7
    # Tm = 64.9 + 41*(yG+zC-16.4)/(wA+xT+yG+zC)
    gc_param = sequence.count('G') + sequence.count('C') - 16.4
    
    calculated_tm = 64.9 + (41 * gc_param / len(sequence))

    return calculated_tm


def self_complementary_check(sequence):
    """
    Align primer to its reverse complement using Bio.Align
    This returns an alignment score. Gaps >1 are not allowed.
    """

    if not isinstance(sequence, Seq):
        sequence = Seq(sequence.upper()) # create seq object if needed

    rev_com = Seq.reverse_complement(sequence)

    aligner = Align.PairwiseAligner() # create a pairwise aligner object
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -10 # this prevents gap extension
    
    max_score = 0
    alignments = aligner.align(str(sequence), str(rev_com)) # perform alignment

    for align in alignments:
        if align.score > max_score:
            max_score = align.score
    
    if max_score > 8:
        if debug:
            print(f"Warning: Primer has a self complentary score of {max_score}", file = sys.stderr)

    return(max_score) # return the highest self complement alignment score



def  check_gc_clamp(sequence):
    """
    Checks for a gc clamp at the end of a primer sequence.
    Takes a string sequence as input and returns a value of True or False.
    Please provide the primer in 5'->3' oreintation
    """
    if isinstance(sequence, Seq):
        sequence = str(sequence) # create string object if needed
    sequence = sequence.upper()
    seq_end = sequence[-5:] # get the last 5 nucleotides
    seq_clamp = sequence[-2:] # get the last 2 nucleotides

    # Count G or C nucleotides
    gc_end_total = seq_end.count('G') + seq_end.count('C')
    gc_clamp_total = seq_clamp.count('G') + seq_clamp.count('C')

    if gc_end_total > 4 or gc_clamp_total < 1:
        if debug:
            print(f"Warning: GC clamp error: {gc_clamp_total} {gc_clamp_total}", file = sys.stderr)
        return False

    return True



def check_repeats(sequence):
    """
    Check for repeated nucleotides in the sequence.
    Detects single nucleotide repeats (4+ of same nucleotide) and dinucleotide repeats (4+ of same dinucleotide).
    Does not detect trinucleotide or longer repeats.
    Returns True if no problematic repeats are found.
    """
    if isinstance(sequence, Seq):
        sequence = str(sequence) # create string object if needed
    sequence = sequence.upper()
    
    # Find single nucleotide repeats (4+ of same nucleotide)
    # Examples: "AAAA", "TTTT", "CCCC", "GGGG"
    single_nt_repeats = re.findall(r'([ATCG])\1{3,}', sequence)
    
    # Find dinucleotide repeats (4+ of same dinucleotide)
    # Examples: "ATATATAT", "GCGCGCGC", "TATATATA"
    dinucleotide_repeats = re.findall(r'([ATCG]{2})\1{3,}', sequence)
    
    # Combine both types of repeats
    all_repeats = single_nt_repeats + dinucleotide_repeats
    
    if all_repeats:
        if debug:
            print(f'Warning: Found repeated sequences: {all_repeats}', file = sys.stderr)
        return False
    
    return True



def check_single_primer(sequence, min_gc = 0.4, max_gc = 0.6, self_comp_max = 8):
    """
    Uses functions from primer_check.py to check candidate 'PCR primer' quality.
    Checks for: GC content, self-complementation, sequence repeats, and GC clamp.
    If all checks are good, the function returns the melting temp of the primer (Tm).
    """
    primer_tm = calculate_tm(sequence) # get melting temp

    gc_content = min_gc < calculate_gc_content(sequence) < max_gc # check GC content
    self_comp = self_complementary_check(sequence) <= self_comp_max # check self-complementation
    nt_repeats = check_repeats(sequence) # check for repetitive sequences
    gc_clamp = check_gc_clamp(sequence) # check for a GC clamp

    if not (gc_content and self_comp and nt_repeats and gc_clamp):
        return -1
    
    return primer_tm




def main(sequence):

    gc_content = calculate_gc_content(sequence)
    tm = calculate_tm(sequence)
    self_comp_score = self_complementary_check(sequence)
    gc_clamp = check_gc_clamp(sequence)
    repeat_check = check_repeats(sequence)

    print(fr'GC content is {gc_content:.4} and Tm is {tm:.4}.')

    if self_comp_score > 8:
        print(f'Primer has a high self complement score: {self_comp_score}')

    if not gc_clamp:
        print(f"3' end of primer has bad GC conetent")

    if not repeat_check:
        print(f"Primer has a repetetive sequence")


if __name__ == '__main__':
    main(sys.argv[1])
    debug = True
