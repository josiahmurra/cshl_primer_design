#!/usr/bin/env python

import sys
from step1_read_edirect import read_edirect
from primer_generation import generate_primer_dict
from primer_generation import identify_primer_sharing


# return suitable reverse primers for shared forward primers
# start with a list of objects in the primer class

def main():
    textfile = sys.argv[1]
    seqlist = read_edirect(textfile)
    primer_dict = generate_primer_dict(seqlist)
    shared_primers = identify_primer_sharing(primer_dict)

dict_primer_pair = {}
def search_rvs_primer(shared_primers):
    for primer in shared_primers:
        print(primer)


if __name__ == '__main__':
    main()