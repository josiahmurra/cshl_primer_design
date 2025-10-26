#!/usr/bin/env python

import sys
import re
import pprint
from Bio import SeqIO
from step1_read_edirect import read_edirect
from step1_read_edirect import align_muscle


def main():
    textfile = sys.argv[1]
    seqlist = read_edirect(textfile)
    alignmentlist = align_muscle(seqlist)
    # the sequence index according to how the accession ids were read in the sequence_files.txt
    sequence_index = {}
    for index, record in enumerate(alignmentlist):
        sequence_index[record.id] = index
    unique_mismatch_targetseq(alignmentlist)
    # implement a counter for the target sequence
    # it is better to make a set instead of calling it as an ATCG string because it is more flexible and runs faster
    # pos_index refers to the multiple sequence alignment, which is different from the target index
    nts = set("ATCGNatcgn")
    target_index = 0 
    def convert_target_to_MSA(target_index):
        for pos_index in range(len(alignmentlist[0].seq)): # using the first sequence as the target sequence, attribute called "index" in primer class
            if alignmentlist[0].seq[pos_index] in nts:
                target_index += 1
            target_MSA_conversion = pos_index - target_index 
            print (target_MSA_conversion, target_index, pos_index)
        return convert_target_to_MSA
    

def unique_mismatch_targetseq(alignmentlist): 
    # creates a list of empty sets across the length of the sequence in the alignment
    # nt_diff_list is a list of dictionary with the key as the nucleotide and the value contains a list of seq_index. 
    nt_diff_list = [dict() for pos_index in range(len(alignmentlist[0].seq))]
    # there are two ways, range and enumerate, to get the indices.
    # the first forloop assigns a number for each sequence in the group.
    # the second forloop assigns a number for each position throughout the sequence.
    for pos_index in range(len(alignmentlist[0].seq)):
        for seq_index, seq_rec in enumerate(alignmentlist):
            nt = seq_rec.seq[pos_index]
            # print((f"Nucleotide at index {pos_index}: {nt}"))
            # if the nucleotide is not in the specific column of the alignment, there will be no value assigned to the dictionary.
            if nt not in nt_diff_list[pos_index]:
            # this is to initiate a list, we want to append to the list, not iterate over the list when using a forloop.
                nt_diff_list[pos_index][nt] = []
            nt_diff_list[pos_index][nt].append(seq_index)
            if len(nt_diff_list[pos_index][nt]) == 1:
                print(pos_index, nt_diff_list[pos_index][nt], nt)
    # pprint.pprint(nt_diff_list) 





    

    





if __name__ == '__main__':
    main()


