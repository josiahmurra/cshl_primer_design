#!/usr/bin/env python

import sys
import re
import pprint
from Bio import SeqIO
from step1_read_edirect import read_edirect
from step1_read_edirect import align_muscle
from primer_generation import generate_primer_dict, identify_primer_sharing



def main():
    textfile = sys.argv[1]
    seqlist = read_edirect(textfile)
    primer_dict = generate_primer_dict(seqlist)

    # Find shared and unique primers
    shared_primers, unique_primers = identify_primer_sharing(primer_dict)
    print(f"Found {len(shared_primers)} primers shared across all sequences")
    print(f"Found {len(unique_primers)} primers unique to single sequences")

    # Perform alignment
    alignmentlist = align_muscle(seqlist)

    # the sequence index according to how the accession ids were read in the sequence_files.txt
    # sequence_index = {}
    # for index, record in enumerate(alignmentlist):
    #     sequence_index[record.id] = index

    # create an index conversion dictionary
    index_conversion_dict = convert_target_to_MSA(alignmentlist, 2) # just test 2 index for now

    
    # identify unique loci
    unique_loci_dict = unique_mismatch_targetseq(alignmentlist, index_conversion_dict, unique_primers)

    
    
    

def convert_target_to_MSA(alignmentlist, target_index):
    target_index = 0 # implement a counter for the target sequence
    nts = set("ATCGNatcgn") # it is better to make a set
    msa_sequence = alignmentlist[target_index].seq

    conversion_dict = {}

    for msa_index in range(len(msa_sequence)): # using the first sequence as the target sequence, attribute called "index" in primer class
        
        conversion_dict[msa_index] = target_index

        if msa_sequence[msa_index] in nts:
            target_index += 1
            
        print (target_index, msa_index)

    return conversion_dict 
    

def unique_mismatch_targetseq(alignmentlist): 
    # creates a list of empty dictionaries with same lenght as the alignment
    # our goal is to make the key a nucleotide and the value a list of seq_index. 
    nt_diff_list = [dict() for msa_index in range(len(alignmentlist[0].seq))]
    
    unique_loci_dict = {} # this dict is less comprehensive. only unique sites.
    
    # the first forloop assigns a number for each position throught the sequence
    # the second forloop .
    unique_pos = 0

    for msa_index in range(len(alignmentlist[0].seq)):
        
        # enumerate alignment list so that seq_index is 0,1,2,3,4 and seq_rec is a seqObject
        for seq_index, seq_rec in enumerate(alignmentlist):

            if seq_index not in unique_loci_dict:
                unique_loci_dict[seq_index] = {}
            
            nt = seq_rec.seq[msa_index] # get nt at this position in msa
            #print((f"Nucleotide at index {msa_index}: {nt}"))

            # if the nucleotide is not in the specific column of the alignment, 
            # there will be no value assigned to the dictionary.
            # Add a nt to our empty nt_dff_list at the right index
            # only add if we havent already added
            if nt not in nt_diff_list[msa_index]:
                nt_diff_list[msa_index][nt] = [] # initiate a list
            nt_diff_list[msa_index][nt].append(seq_index) # append to the list

        # if length of list is 1, there was only one sequence with this nt at this position
        # save unique index positions in the unique loci list
        for nuc in nt_diff_list[msa_index]:
            #print(nt_diff_list[msa_index][nuc])
            if len(nt_diff_list[msa_index][nuc]) == 1:
                # these are unique
                unique_pos += 1
                index_key = nt_diff_list[msa_index][nuc][0]
                unique_loci_dict[index_key][msa_index] = nuc

    # pprint.pprint(unique_loci_dict)
    # print(f' Found {unique_pos} unique positions') 
    return unique_loci_dict

# def target_unique_loci(unique_loci_dict, catch_sequence_index, unique_primers):


#     unique_pos_list = unique_loci_dict[catch_sequence_index].keys()
#     for ind, pos, in enumerate(unique_pos_list):
#         if [ind+1] - pos < 18




if __name__ == '__main__':
    main()


