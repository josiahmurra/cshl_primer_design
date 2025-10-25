#!/usr/bin/env python

import sys
import pprint
from Bio.Seq import reverse_complement
from step1_read_edirect import read_edirect
from primer_check import check_single_primer

# define the PrimerObject class
class PrimerObject:

    # set attributes of PrimerObject
    def __init__(self, seq, start, end, tm, direction, index = -1):
        self.seq = seq
        self.start = start
        self.end = end
        self.index = index
        self.tm = tm
        self.direction = direction




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
                        index = index
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
                    
                    primer_dict[record.id][candidate_primer.seq] = value

        print(f"For {record.id} there were {bad_primer_count} bad primers and {duplicate_primer_count} duplicates")
    return(primer_dict)

                



def main():
    textfile = sys.argv[1]
    seqlist = read_edirect(textfile)
    print(seqlist)
    primer_dict = generate_primer_dict(seqlist)

    print(len(seqlist[0]))
    print(len(primer_dict['NM_001008221.1']))
    #pprint.pprint(primer_dict['NM_001008221.1'])




if __name__ == "__main__":
    main()

