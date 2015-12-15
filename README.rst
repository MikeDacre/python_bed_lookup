#################
Python Bed Lookup
#################

A simple data structure to hold a bed file, allows the user to
search by location for a gene name.

For example:

    from bed_lookup import BedFile
    b = BedFile('my_bed.bed')
    gene = b.lookup('chr3', 1000104)

This module requires cython, and should work with recent versions of
python2 and python3, although python3 will be faster.
