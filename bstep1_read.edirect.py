#!/usr/bin/env python

# efetch will allow you to fetch the sequences using accession numbers.

import sys
import re
from Bio import Entrez
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def main():
    textfile = sys.argv[1]
    seqlist = read_edirect(textfile)
    print(seqlist)

# seqfile_text is the filepath to a text file that contains the accession numbers extracted from NCBI Refseq database. 
def read_edirect(seqfile_textfile):
    with open(seqfile_textfile, "r") as seqfile:
        seqid_list = []
        for line in seqfile:
            line_stripped = line.rstrip()
            seqid_list.append(line_stripped)    
    # Entrez programming utilitis (E-utilities) requires an email address.
    Entrez.email = 'chocoshi@gmail.com'
    # id= this accepts seqid_list, which is generated prior, and returns a list of accession numbers.
    # rettype="gb": Retrieves the record in GenBank format.
    # retmode="text": Returns the data as plain text.
    genbank_fetch_handle = Entrez.efetch(db = "nucleotide", id = seqid_list, rettype = "gb", retmode = "text" )
    seq_object_list = list(SeqIO.parse(genbank_fetch_handle, "genbank"))
    
    # Resources:
    # https://www.geeksforgeeks.org/python/biopython-sequence-alignment/
    # https://biopython.org/docs/1.76/api/Bio.AlignIO.html
    # https://ugene.net/docs/alignment-editor/
    
    seqrecord_list = []
    alignment_file = ()
    for record in seq_object_list:
        record = SeqRecord(Seq(record.seq), id=record.id)
        seqrecord_list.append(record)
    alignment_file = MultipleSeqAlignment(seqrecord_list)

    return alignment_file

if __name__ == '__main__':
    main()