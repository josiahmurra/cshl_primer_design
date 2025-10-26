#!/usr/bin/env python

import sys
import subprocess
from tempfile import NamedTemporaryFile, gettempprefix

# efetch will allow you to fetch the sequences using accession numbers.

from Bio import Entrez
from Bio import SeqIO
from Bio import AlignIO

def main():
    textfile = sys.argv[1]
    seqlist = read_edirect(textfile)
    alignmentlist = align_muscle(seqlist)
    # print(alignmentlist)
    # print(alignmentlist[0].seq)

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
    return seq_object_list

# using the temporary files because MuscleCommandline is deprecated
## Q: What about pymuscle5?
def align_muscle(seq_object_list):
    temp_align = gettempprefix()
    with NamedTemporaryFile(mode = "wt", prefix='muscle', delete_on_close=False) as temp_muscle:
        # writing in fasta format from the list of objects in seq_object_list, and save to temp_muscle
        for seq_object in seq_object_list:
            print(f">{seq_object.id}\n{seq_object.seq}", file = temp_muscle)
            # print(f">{seq_object.id}\n{seq_object.seq}")
        subprocess.run(["muscle", "-in", temp_muscle.name, "-out", temp_align])
    alignment_list = list(AlignIO.read(temp_align, "fasta"))
    return alignment_list

if __name__ == '__main__':
    main()







