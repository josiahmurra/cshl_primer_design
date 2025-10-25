#!/usr/bin/env python

import sys
import subprocess
import tempfile
from tempfile import NamedTemporaryFile, gettempprefix

# efetch will allow you to fetch the sequences using accession numbers.

from Bio import Entrez
from Bio import SeqIO
from Bio import AlignIO

def main():
    textfile = sys.argv[1]
    seqlist = read_edirect(textfile)
    alignment_list = muscle_seqgroup(seqlist)
    for alignment in alignment_list:
        print(alignment)


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
def muscle_seqgroup(seq_object_list):
    temp_align = gettempprefix()
    with NamedTemporaryFile(mode = "wt", prefix='muscle', delete_on_close=False) as temp_muscle:
        # writing in fasta format from the list of objects in seq_object_list, and save to temp_muscle
        for seq_object in seq_object_list:
            print(f">{seq_object.id}\n{seq_object.seq}", file = temp_muscle)
            print(f">{seq_object.id}\n{seq_object.seq}")
        subprocess.run(["muscle", "-in", temp_muscle.name, "-out", temp_align])
    alignments = list(AlignIO.read(temp_align, "fasta"))
    return(alignments)

       




#     muscle_object_align = MuscleCommandline(input=seq_object_list)
#     print(muscle_object_align)
#     stdout, stderr = muscle_object_align()
#     alignment = AlignIO.read(stdout, "fasta")
#     AlignIO.write(alignment, "seqgroup_output_alignment.fasta", "fasta")
#     print("seqgroup_output_alignment.fasta")

# muscle_seqgroup(read_edirect('sequence_files.txt'))







if __name__ == '__main__':
    main()







