#################
Python Bed Lookup
#################

A simple data structure to hold a bed file, allows the user to
search by location for a gene name.

For example::

    from bed_lookup import BedFile
    b = BedFile('my_bed.bed')
    gene = b.lookup('chr3', 1000104)

This module requires cython, and should work with recent versions of
python2 and python3, although python3 will be faster.

It makes use of a cython optimized dictionary lookup for small bed files
and sqlite for larger bed files. Which backend is being used is transparent
to the user, simply use the ``lookup()`` function as demonstrated in the
example above.
