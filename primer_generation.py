#!/usr/bin/env python

import sys
import pprint
from Bio.Seq import reverse_complement
from step1_read_edirect import read_edirect
from primer_check import check_single_primer


# define the PrimerObject class
class PrimerObject:
    # set attributes of PrimerObject
    def __init__(self, seq, start, end, tm, direction, seq_id, index = -1):
        self.seq = seq
        self.start = start
        self.end = end
        self.index = index
        self.tm = tm
        self.direction = direction
        self.seq_id = seq_id




def generate_primer_dict(seq_record_list, min_len = 18, max_len = 25):
    """
    Takes as input a list of SeqRecord objects. Returns a dictionary where
    sequence ID (from SeqRecord) is the key and the value is a list of PrimerObjects
    """
    primer_dict = {} # initate primer dictionary

    # using index to keep track of the order of sequence list
    for index, record in enumerate(seq_record_list):
        # add id to primer dictionary and initiate sub-dictionary
        primer_dict[record.id] = {} 
        good_primer_count = 0
        bad_primer_count = 0
        duplicate_primer_count = 0
        # considering primers of length 18-25 by default
        for primer_length in range(min_len , max_len + 1): 

            loops = (len(record.seq) - primer_length) + 1 # total kmers to identify

            for i in range(loops):
                kmer = str(record.seq[i:i + primer_length]) # get current kmer
                kmer_dict = {} # initate dictionary for forward and revcom
                
                kmer_dict['forward'] = kmer
                kmer_dict['reverse'] = reverse_complement(kmer)

                for strand in kmer_dict:
                    kmer_tm = check_single_primer(kmer_dict[strand])

                    candidate_primer = PrimerObject(
                        seq = kmer_dict[strand], 
                        start = i, 
                        end = i + primer_length,
                        tm = kmer_tm,
                        direction = strand,
                        index = index,
                        seq_id = record.id
                    )
                    value = candidate_primer

                    # make sure the primer sequence is unique
                    if candidate_primer.seq in primer_dict[record.id]:
                        duplicate_primer_count += 1
                        value = None

                    # ensure primer Tm is 50-80 and met other check criteria
                    if kmer_tm > 80 or kmer_tm < 50:
                        bad_primer_count += 1
                        value = None
                    else:
                        good_primer_count += 1
                    
                    primer_dict[record.id][candidate_primer.seq] = value

        print(f"Checked {len(primer_dict[record.id])} primers for {record.id}")
        print(f"Good Quality: {good_primer_count}") 
        print(f"Unsuitable:   {bad_primer_count}")
        print(f"Duplicates:   {duplicate_primer_count}")
        print()
    return(primer_dict)

                


def identify_primer_sharing(primer_dict):
    """
    Identify primers that are shared across ALL sequences and primers UNIQUE to single sequences.
    Returns two dictionaries: shared_primers and unique_primers.
    """
    # Dictionary to track which sequences contain each primer
    primer_occurrences = {}
    
    total_sequences = len(primer_dict) # Total number of sequences
    
    # Iterate through each sequence in the primer dictionary
    for seq_id, primers in primer_dict.items():
        # Iterate through each primer in the sequence
        for primer_seq, primer_obj in primers.items():
            # Skip None values (bad primers)
            if primer_obj is not None:
                if primer_seq not in primer_occurrences:
                    primer_occurrences[primer_seq] = []
                primer_occurrences[primer_seq].append(primer_obj)
    
    # Find primers that appear in ALL sequences
    shared_primers = {}
    for primer_seq, primer_obj_list in primer_occurrences.items():
        # Get unique sequence indices (primer_obj has an 'index' attribute)
        unique_indices = set([primer_obj.index for primer_obj in primer_obj_list])
        # Only keep primers that appear in all sequences
        if len(unique_indices) == total_sequences:
            shared_primers[primer_seq] = primer_obj_list[0] # default to zero index seq_id
    
    # Find primers that appear in ONLY ONE sequence
    unique_primers = {}
    for primer_seq, primer_obj_list in primer_occurrences.items():
        # Get unique sequence indices (primer_obj has an 'index' attribute)
        unique_indices = set([primer_obj.index for primer_obj in primer_obj_list])
        # Only keep primers that appear in exactly one sequence
        if len(unique_indices) == 1:
            unique_primers[primer_seq] = primer_obj_list[0]
    
    return shared_primers, unique_primers

    

def main():
    textfile = sys.argv[1]
    seqlist = read_edirect(textfile)
    primer_dict = generate_primer_dict(seqlist)

    # Find shared and unique primers
    shared_primers, unique_primers = identify_primer_sharing(primer_dict)
    print(f"Found {len(shared_primers)} primers shared across all sequences")
    print(f"Found {len(unique_primers)} primers unique to single sequences")

    # Test message (don't need in full script)
    n = 0
    for primer in unique_primers:
        if n == 1:
            primer_obj = unique_primers[primer]
            print(primer_obj.seq)
            print(primer_obj.index)
            print(primer_obj.start, primer_obj.end)
            print(primer_obj.seq_id)
        n += 1
 
if __name__ == "__main__":
    main()

